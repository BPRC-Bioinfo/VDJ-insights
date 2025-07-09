"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

from pathlib import Path
import subprocess
import pandas as pd
from typing import Union
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

from .util import make_dir, calculate_available_resources

from .property import log_error
from .logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def write_fasta(sub_df: pd.DataFrame, out_fasta: Path):
    """
    Schrijft een multi-FASTA-bestand voor alle rijen in sub_df.

    Args:
        sub_df (pd.DataFrame): De subset DataFrame met kolommen 'name', 'sequence_clean', 'start', 'stop', 'strand', 'fasta-file', 'tool'.
        out_fasta (Path): Pad naar het FASTA-uitvoerbestand.
    """
    with out_fasta.open("w") as file:
        for _, row in sub_df.iterrows():
            header = (
                f"{row['name']}___{row['start']}___{row['stop']}___"
                f"{row['strand']}___{row['fasta-file']}___{row['tool']}"
            )
            file.write(f">{header}\n{row['sequence_clean']}\n")


@log_error()
def make_blast_db(cwd: Path, library: str, verbose: bool) -> Path:
    """
    Checks if the BLAST database directory exists. If not, creates the directory
    and initializes the database using the provided library file. Uses the `makeblastdb`
    command with the following parameters:

    - `-in`: Specifies the input FASTA file to use as the reference.
    - `-dbtype nucl`: Indicates that the database will be nucleotide-based.
    - `-out`: Specifies the output path and prefix for the BLAST database files.

    Logs if the database is created successfully or if it already exists. Returns the
    path to the BLAST database directory.

    Args:
        cwd (Path): The current working directory of the user.
        library (str): The path to the FASTA reference file used for creating the BLAST database.

    Returns:
        Path: The path to the BLAST database directory.

    Raises:
        subprocess.CalledProcessError: If the `makeblastdb` command fails.
        OSError: If the directory cannot be created or accessed.
    """
    blast_db_path = cwd / "tmp" / "mapping" / "blast_db"
    if not blast_db_path.exists():
        make_dir(blast_db_path)
        cmd = f"makeblastdb -in {library} -dbtype nucl -out {blast_db_path}/blast_db"
        subprocess.run(cmd,
                       shell=True,
                       check=True,
                       stdout=subprocess.PIPE if not verbose else None,
                       stderr=subprocess.PIPE if not verbose else None
                       )
    return blast_db_path


def run_blast(df: pd.DataFrame, db_path: Path, tmp_dir: Path, sample: str, blast_cols: list, threads: int, verbose: bool) -> list[list[str]]:
    df_s = df[df["sample"] == sample]
    df_s_short = df_s[df_s["sequence_clean"].str.len() <= 50]
    df_s_long = df_s[df_s["sequence_clean"].str.len() > 50]

    safe_sample = sample.replace("/", "_").replace(" ", "_").replace("|", "_")
    fasta_short = tmp_dir / f"{safe_sample}_queries_short.fasta"
    fasta_long = tmp_dir / f"{safe_sample}_queries_long.fasta"
    blast_short_out = tmp_dir / f"{safe_sample}_blast_short.tsv"
    blast_long_out = tmp_dir / f"{safe_sample}_blast_long.tsv"

    if not df_s_long.empty:
        write_fasta(df_s_long, fasta_long)
    if not df_s_short.empty:
        write_fasta(df_s_short, fasta_short)

    rows = []
    if not df_s_long.empty:
        cmd_long = (
            f"blastn -task megablast "
            f"-query {fasta_long} "
            f"-db {db_path}/blast_db "
            f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs btop' "
            f"-perc_identity 50 "
            f"-max_target_seqs 5 "
            f"-out {blast_long_out} "
            f"-num_threads {threads}"
        )
        subprocess.run(cmd_long,
                       shell=True,
                       check=True,
                       stdout=subprocess.PIPE if not verbose else None,
                       stderr=subprocess.PIPE if not verbose else None,
                       )

        with open(blast_long_out, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) == len(blast_cols):
                    rows.append(fields)

    if not df_s_short.empty:
        extra = "-word_size 7 -evalue 1000 -penalty -3 -reward 1 -gapopen 5 -gapextend 2 -dust no"
        cmd_short = (
            f"blastn -task megablast "
            f"-query {fasta_short} "
            f"-db {db_path}/blast_db "
            f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs btop' "
            f"-perc_identity 50 "
            f"-max_target_seqs 5 "
            f"{extra} "
            f"-out {blast_short_out} "
            f"-num_threads {threads}"
        )
        subprocess.run(cmd_short,
                       shell=True,
                       check=True,
                       stdout=subprocess.PIPE if not verbose else None,
                       stderr=subprocess.PIPE if not verbose else None,
                       )

        with open(blast_short_out, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) == len(blast_cols):
                    rows.append(fields)
    return rows


@log_error()
def run_blast_per_sample(df: pd.DataFrame, db_path: Path, output_csv: Path, threads: int, verbose: bool) -> None:
    """
       Batch BLAST per sample: splits input by sample and for each sample:
         1) Creates FASTA files for long and short sequences.
         2) Performs BLAST (standard parameters for long sequences, adjusted parameters for short ones).
         3) Collects all results in memory.
       In the end, all results are combined, filtered for 100% query coverage,
       and written to output_csv.

       Args:
           df (pd.DataFrame): Input DataFrame with columns 'reference', 'start', 'stop', 'name', 'strand', 'sequence', 'tool', 'fasta-file'.
           db_path (Path): Path to the BLAST database directory (should contain 'blast_db').
           output_csv (Path): Output path where the combined results will be saved as CSV.
           threads (int): Number of threads to use for BLAST.
           verbose (bool): Whether to show BLAST output in the terminal.
    """
    tmp_dir = db_path.parent / "batch_fasta"
    make_dir(tmp_dir)

    df["sequence_clean"] = df["sequence"].str.replace("-", "", regex=False)
    df["sample"] = df["reference"].str.split("__", n=1).str[0]
    samples = df["sample"].unique()

    blast_cols = ["query", "subject", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "query seq", "subject seq", "query cov", "btop"]
    aggregated_rows = []

    max_jobs = calculate_available_resources(max_cores=threads, threads=4, memory_per_process=12)
    with ThreadPoolExecutor(max_workers=max_jobs) as executor:
        futures = [executor.submit(run_blast, df, db_path, tmp_dir, sample, blast_cols,  threads, verbose) for sample in samples]
        with tqdm(total=len(samples), desc='Blast reevaluation:', unit="sample") as pbar:
            for future in as_completed(futures):
                aggregated_rows += future.result()
                pbar.update(1)

    if aggregated_rows:
        blast_results = pd.DataFrame(aggregated_rows, columns=blast_cols)
    else:
        blast_results = pd.DataFrame(columns=blast_cols)

    blast_results["query cov"] = pd.to_numeric(blast_results["query cov"], errors="coerce")
    blast_results = blast_results.query("`query cov` == 100")
    blast_results.to_csv(output_csv, index=False)


@log_error()
def blast_main(df: pd.DataFrame, blast_file: Union[str, Path], library: str, threads: int, verbose: bool) -> None:
    """
    Main entry point for BLAST processing: create DB, run per-sample BLAST, and save results.

    Args:
        df (pd.DataFrame): DataFrame with sequence metadata and sequences.
        blast_file (str | Path): Path for the output CSV file.
        library (str): Path to the reference FASTA for the BLAST database.
        threads (int): Number of CPU cores for processing.
        verbose (bool): If False, suppress subprocess output.
    """
    cwd = Path.cwd()
    db_path = make_blast_db(cwd, library, verbose)
    run_blast_per_sample(df, db_path, Path(blast_file), threads, verbose)
