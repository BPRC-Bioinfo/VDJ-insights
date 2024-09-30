from concurrent.futures import ThreadPoolExecutor
import subprocess
from tempfile import NamedTemporaryFile as Ntf
from pathlib import Path
import pandas as pd
from tqdm import tqdm

from logger import custom_logger
from util import make_dir
from property import log_error


logger = custom_logger(__name__)


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
        library (str): The library file name to be used for creating the BLAST database.

    Returns:
        Path: The path to the BLAST database directory.

    Raises:
        subprocess.CalledProcessError: If the `makeblastdb` command fails.
        OSError: If the directory cannot be created or accessed.
    """
    blast_db_path = cwd / "mapping" / "blast_db"
    if not blast_db_path.exists():
        reference = cwd / library
        make_dir(blast_db_path)
        command = f"makeblastdb -in {reference} -dbtype nucl -out {blast_db_path}/blast_db"
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info("BLAST database created successfully.")
    else:
        logger.info("BLAST database already exists.")
    return blast_db_path


def construct_blast_command(fasta_file_path: str | Path,
                            database_path: str | Path,
                            identity_cutoff: int,
                            output_file_path: str | Path,
                            length: int,
                            LENGTH_THRESHOLD: int = 15) -> str:
    """
    Constructs a BLAST command string for sequence alignment. This function builds 
    the command using the `blastn` tool with the following parameters:

    - `-task megablast`: Uses the `megablast` algorithm optimized for highly similar sequences.
    - `-query`: Specifies the path to the input FASTA file.
    - `-db`: Specifies the path to the BLAST database.
    - `-outfmt '6 ...'`: Defines custom output columns, including alignment details like 
      query ID, subject ID, percent identity, alignment length, mismatches, gap opens, 
      query and subject start and end positions, e-value, bit score, query and subject 
      sequences, and query coverage.
    - `-perc_identity`: Sets the minimum percentage identity for alignments.
    - `-out`: Specifies the path to the output file for the BLAST results.

    If the sequence length is below the defined `LENGTH_THRESHOLD`, additional parameters 
    are included to fine-tune the alignment:

    - `-word_size 7`: Sets the word size (k-mer length) to 7, enhancing sensitivity for short alignments.
    - `-evalue 1000`: Sets a high e-value threshold to retain more potential matches.
    - `-max_target_seqs 100`: Limits the number of target sequences in the output to 100.
    - `-penalty -3`: Applies a penalty of -3 for mismatches in the alignment.
    - `-reward 1`: Rewards a match with a score of 1.
    - `-gapopen 5`: Sets the gap opening penalty to 5.
    - `-gapextend 2`: Sets the gap extension penalty to 2.
    - `-dust no`: Disables the filtering of low-complexity regions.

    Args:
        fasta_file_path (Path): Path to the temporary FASTA file.
        database_path (Path): Path to the BLAST database.
        identity_cutoff (int): Minimum percentage identity for the BLAST alignment.
        output_file_path (Path): Path to the output BLAST result file.
        length (int): Length of the sequence being aligned.
        LENGTH_THRESHOLD (int, optional): Threshold for applying additional BLAST parameters. Defaults to 15.

    Returns:
        str: The constructed BLAST command string.
    """
    blast_columns = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs"
    extra = "-word_size 7 -evalue 1000 -max_target_seqs 100 -penalty -3 -reward 1 -gapopen 5 -gapextend 2 -dust no"
    command = f"blastn -task megablast -query {fasta_file_path} -db {database_path}/blast_db -outfmt '{
        blast_columns}' -perc_identity {identity_cutoff} -out {output_file_path}"
    if length <= LENGTH_THRESHOLD:
        command += f" {extra}"
    return command


@log_error()
def execute_blast_search(row: pd.Series, database_path: Path, identity_cutoff: int) -> str:
    """
    Executes a BLAST search for a given row from a DataFrame. Creates a temporary FASTA 
    file from the row data and runs the BLAST command against the specified database.

    Extracts the necessary data from the DataFrame row, including the sequence and 
    associated metadata. A temporary FASTA file is generated with a header that includes 
    information like the sequence name, start and stop positions, strand, file name, and haplotype.

    The BLAST search is executed using the constructed command, and the results are 
    saved in a temporary output file. If an error occurs during the execution, it is logged.

    Args:
        row (pd.Series): A single row from the DataFrame containing sequence data.
        database_path (Path): Path to the BLAST database.
        identity_cutoff (int): Minimum percentage identity for the BLAST alignment.

    Returns:
        str: Path to the temporary file containing the BLAST search results.

    Raises:
        subprocess.CalledProcessError: If the BLAST command fails.
    """
    header, sequence, start, stop, fasta_file_name, strand, haplotype = row["name"], row[
        "sequence"], row["start"], row["stop"], row["fasta-file"], row["strand"], row["haplotype"]

    with Ntf(mode='w+', delete=False, suffix='.fasta') as fasta_temp:
        fasta_header = f">{header}:{start}:{stop}:{
            strand}:{fasta_file_name}:{haplotype}\n"
        fasta_temp.write(fasta_header + sequence + "\n")
        fasta_temp.flush()

        blast_result_path = Ntf(mode='w+', delete=False, suffix='_blast.txt').name
        command = construct_blast_command(fasta_temp.name, database_path, identity_cutoff, blast_result_path, len(sequence))
        subprocess.run(command, shell=True, check=True)
        logger.info(f"Executed BLAST search for {header}")
    return blast_result_path



@log_error()
def aggregate_blast_results(dataframe: pd.DataFrame, database_path: Path, CUTOFFS=[100, 75, 50]) -> pd.DataFrame:
    """
    Iterates over a set of identity cutoffs (100%, 75%, 50%) and performs parallel BLAST 
    searches for each row in the input DataFrame. Uses `execute_blast_search()` to obtain 
    the results for each row, aggregating them into a single DataFrame.

    For each cutoff, the function creates temporary result files for storing the BLAST outputs. 
    Once the results are read, the temporary files are deleted. The resulting DataFrame 
    includes an additional 'cutoff' column indicating the identity cutoff used for each search.

    Args:
        dataframe (pd.DataFrame): DataFrame with rows representing sequence mappings.
        database_path (Path): Path to the BLAST database.
        CUTOFFS (list of int, optional): List of identity cutoffs to use. Defaults to [100, 75, 50].

    Returns:
        pd.DataFrame: DataFrame containing the aggregated BLAST search results.

    Raises:
        OSError: If any temporary file cannot be created or deleted.
        subprocess.CalledProcessError: If a BLAST command fails during execution.
    """
    aggregated_results = pd.DataFrame()
    total_searches = len(dataframe) * len(CUTOFFS)
    with tqdm(total=total_searches, desc="Running all BLAST searches", unit="search") as pbar:
        with ThreadPoolExecutor() as executor:
            for cutoff in CUTOFFS:
                futures = {executor.submit(
                    execute_blast_search, row, database_path, cutoff): row for _, row in dataframe.iterrows()}
                for future in futures:
                    result_file_path_str = future.result()
                    blast_columns = [
                        'query', 'subject', '% identity', 'alignment length', 'mismatches', 'gap opens',
                        'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score',
                        'query seq', 'subject seq', 'query cov'
                    ]
                    temp_df = pd.read_csv(
                        result_file_path_str, sep='\t', names=blast_columns)
                    if not temp_df.empty:
                        temp_df['cutoff'] = cutoff
                        aggregated_results = pd.concat(
                            [aggregated_results, temp_df], ignore_index=True)
                    Path(result_file_path_str).unlink()
                    pbar.update(1)

    logger.info("Aggregation of BLAST results completed.")
    return aggregated_results


def run_blast_operations(df: pd.DataFrame, db_path: Path, blast_file_path: Path) -> None:
    """
    Runs BLAST operations on the input DataFrame and saves the results to an Excel file.

    Aggregates BLAST results by performing searches with different identity cutoffs. 
    Filters the results to include only alignments with 100% query coverage. Extracts 
    additional columns from the BLAST output and saves the final results to an Excel file.

    Args:
        df (pd.DataFrame): Input DataFrame with sequence mappings.
        db_path (Path): Path to the BLAST database.
        blast_file_path (Path): Path to save the BLAST results Excel file.

    Returns:
        None

    Raises:
        OSError: If any file operation fails during the process.
    """
    blast_results = aggregate_blast_results(df, db_path)
    blast_results['query cov'] = pd.to_numeric(blast_results['query cov'], errors='coerce')
    blast_results = blast_results.query("`query cov` == 100")
    path_df = blast_results['query'].str.split(':', expand=True)
    blast_results[['start', 'stop']] = path_df[[1, 2]]
    blast_results.to_csv(blast_file_path, index=False)
    
    logger.info("BLAST operations completed and results saved to csv.")


def blast_main(df: pd.DataFrame, blast_file: str | Path, library: str) -> None:
    """
    Manages the entire BLAST process, from database creation to saving results.

    It performs the following steps:
        1. Creates the BLAST database if it does not exist.
        2. Runs the sequence alignment operations using BLAST with different identity cutoffs.
        3. Filters and processes the BLAST results to include only alignments with 100% query coverage.
        4. Saves the processed results in an Excel file.

    Args:
        df (pd.DataFrame): DataFrame containing the sequence data.
        blast_file (str): File path to save the BLAST results.
        library (str): Library file used for creating the BLAST database.

    Returns:
        None

    Raises:
        OSError: If any file or directory operation fails.
        subprocess.CalledProcessError: If a BLAST command or database creation fails.
    """
    cwd = Path.cwd()
    db_path = make_blast_db(cwd, library)
    run_blast_operations(df, db_path, Path(blast_file))
    logger.info("BLAST main process completed.")
