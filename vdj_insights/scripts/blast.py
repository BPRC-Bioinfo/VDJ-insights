"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

from pathlib import Path
import subprocess
import pandas as pd
from typing import Union
from tqdm import tqdm

from .util import make_dir
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
            header = row["name"]
            seq = row["sequence_clean"]
            start = row["start"]
            stop = row["stop"]
            strand = row["strand"]
            fasta_fn = row["fasta-file"]
            tool = row["tool"]
            file.write(f">{header}___{start}___{stop}___{strand}___{fasta_fn}___{tool}\n")
            file.write(seq + "\n")


@log_error()
def make_blast_db(cwd: Path, library: str) -> Path:
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
        command = f"makeblastdb -in {library} -dbtype nucl -out {blast_db_path}/blast_db"
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return blast_db_path


@log_error()
def run_blast_per_sample(df: pd.DataFrame, db_path: Path, output_csv: Path, threads: int) -> None:
    """
    Batch-BLAST per sample: split op sample, voor elke sample:
      1) Maak FASTA voor lange en korte sequenties.
      2) Voer BLAST uit (lange met standaard, korte met extra parameters).
      3) Verzamelt alle resultaten in geheugen.
    Uiteindelijk worden alle resultaten samengevoegd, gefilterd op 100% query coverage,
    en weggeschreven naar output_csv.

    Args:
        df (pd.DataFrame): Input DataFrame met kolommen 'reference', 'start', 'stop', 'name', 'strand', 'sequence', 'tool', 'fasta-file'.
        db_path (Path): Pad naar de BLAST-database (map die 'blast_db' bevat).
        output_csv (Path): Pad waar de samengevoegde resultaten als CSV worden weggeschreven.
        threads (int): Aantal threads om BLAST mee te draaien.
    """
    tmp_dir = db_path.parent / "batch_fasta"
    make_dir(tmp_dir)

    df["sequence_clean"] = df["sequence"].str.replace("-", "")
    df["sample"] = df["reference"].str.split("__", n=1).str[0]

    blast_cols = [
        "query", "subject", "% identity", "alignment length",
        "mismatches", "gap opens", "q. start", "q. end",
        "s. start", "s. end", "evalue", "bit score",
        "query seq", "subject seq", "query cov", "btop"
    ]

    aggregated_rows = []

    samples = df["sample"].unique()
    for sample in tqdm(samples, desc="Blast reevaluation", unit="sample"):
        df_s = df[df["sample"] == sample]
        df_s_short = df_s[df_s["sequence_clean"].str.len() <= 50]
        df_s_long  = df_s[df_s["sequence_clean"].str.len() >  50]

        safe_sample = sample.replace("/", "_").replace(" ", "_").replace("|", "_")
        fasta_short = tmp_dir / f"{safe_sample}_queries_short.fasta"
        fasta_long  = tmp_dir / f"{safe_sample}_queries_long.fasta"
        blast_short_out = tmp_dir / f"{safe_sample}_blast_short.tsv"
        blast_long_out  = tmp_dir / f"{safe_sample}_blast_long.tsv"

        if not df_s_long.empty:
            write_fasta(df_s_long, fasta_long)
        if not df_s_short.empty:
            write_fasta(df_s_short, fasta_short)

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
            subprocess.run(cmd_long, shell=True, check=True)

            with open(blast_long_out, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) == len(blast_cols):
                        aggregated_rows.append(fields)

        # 2) Korte sequenties
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
            subprocess.run(cmd_short, shell=True, check=True)

            with open(blast_short_out, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) == len(blast_cols):
                        aggregated_rows.append(fields)
    if aggregated_rows:
        blast_results = pd.DataFrame(aggregated_rows, columns=blast_cols)
    else:
        blast_results = pd.DataFrame(columns=blast_cols)

    blast_results["query cov"] = pd.to_numeric(blast_results["query cov"], errors="coerce")
    blast_results = blast_results.query("`query cov` == 100")
    blast_results.to_csv(output_csv, index=False)


@log_error()
def blast_main(df: pd.DataFrame, blast_file: Union[str, Path], library: str, threads: int) -> None:
    """
    Beheert het volledige BLAST-proces, van database-aanmaak tot weergave van resultaat:

    1. Creëert de BLAST-database indien nog niet aanwezig.
    2. Verdeelt de input-DataFrame op sample en voert voor elke sample BLAST uit (lange en korte sequenties).
    3. Filtert de gecombineerde resultaten op 100% query coverage.
    4. Schrijft het eindresultaat weg als CSV.

    Args:
        df (pd.DataFrame): DataFrame met sequentiegegevens (kolommen 'reference', 'start', 'stop', 'name', 'strand', 'sequence', 'tool', 'fasta-file').
        blast_file (str of Path): Pad voor het uiteindelijke CSV-bestand.
        library (str): Pad naar de FASTA-reference voor de BLAST-database.
        threads (int): Aantal cores om BLAST mee te draaien.
    """
    cwd = Path.cwd()
    db_path = make_blast_db(cwd, library)
    output_csv = Path(blast_file)
    run_blast_per_sample(df, db_path, output_csv, threads)
