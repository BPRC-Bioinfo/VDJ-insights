from concurrent.futures import ThreadPoolExecutor
import subprocess
from tempfile import NamedTemporaryFile as NTF
from pathlib import Path
import pandas as pd
from logger import custom_logger

"""
Used python packages:
    1. pandas
    2. openpyxl
    
Used CLI packages:
    1. blast
"""

# Method for logging current states of the program.
logger = custom_logger(__name__)


def make_dir(dir):
    """
    Create a directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """
    try:
        Path(dir).mkdir(parents=True, exist_ok=True)
        logger.info(f"Directory created or already exists: {dir}")
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")
    return dir


def make_blast_db(cwd, library):
    """
    Create a blast database (db) directory folder to run blast with.
    It checks if the db exist, if exist return the db. 
    Otherwise create the db.

    Args:
        cwd (Path): Path of the users current location.

    Returns:
        db (Path): Path of the blast db.
    """
    blast_db_path = cwd / "mapping" / "blast_db"
    if not blast_db_path.exists():
        reference = cwd / library
        make_dir(blast_db_path)
        command = f"makeblastdb -in {reference} -dbtype nucl -out {blast_db_path}/blast_db"
        try:
            subprocess.run(command, shell=True, check=True)
            logger.info("BLAST database created successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating BLAST database: {e}")
    else:
        logger.info("BLAST database already exists.")
    return blast_db_path


def construct_blast_command(fasta_file_path, database_path, identity_cutoff, output_file_path, length, LENGTH_THRESHOLD=15):
    """
    Construct the BLAST command string based on the given input and
    output paths and a % identity cutoff. 
    It also uses a custom blast_columns.
    - qseqid: Query Sequence ID
    - sseqid: Subject Sequence ID
    - pident: Percentage of Identical Matches
    - length: Alignment Length
    - mismatch: Number of Mismatches
    - gapopen: Number of Gap Openings
    - qstart: Start of Alignment in Query
    - qend: End of Alignment in Query
    - sstart: Start of Alignment in Subject
    - send: End of Alignment in Subject
    - evalue: Expectation Value
    - bitscore: Bit Score
    - qseq: Aligned Part of Query Sequence
    - sseq: Aligned Part of Subject Sequence
    - qcovs: Query Coverage Per Subject

    Args:
        fasta_file_path (Path): Path of the temporary fasta file.
        database_path (Path): Path of the blast db.
        identity_cutoff (int): Integer to indicate what the minimum
        % Identity must be.
        output_file_path (Path): Path of the temporary blast output file.
        length (int): Length of the segment.

    Returns:
        command (str): Constructed blast command. 
    """
    blast_columns = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs"
    extra = "-word_size 7 -evalue 1000 -max_target_seqs 100 -penalty -3 -reward 1 -gapopen 5 -gapextend 2 -dust no"
    command = f"blastn -task megablast -query {fasta_file_path} -db {database_path}/blast_db -outfmt '{blast_columns}' -perc_identity {identity_cutoff} -out {output_file_path}"
    if length <= LENGTH_THRESHOLD:
        command += f" {extra}"
    logger.debug(f"Constructed BLAST command: {command}")
    return command


def execute_blast_search(row, database_path, identity_cutoff) -> str:
    """
    Executes a BLAST search for a given row from a DataFrame, 
    creating a temporary FASTA file from the row data and running the
    BLAST command against a specified database.
    The output is a temporary file containing BLAST results.


    Args:
        row (Series): Row of the df.
        database_path (Path): Path of the blast db.
        identity_cutoff (int): Integer to indicate what the minimum
        % Identity must be.

    Returns:
        blast_result_path (str): Path to the temporary blast output file.
    """
    header, sequence, start, stop, fasta_file_name, strand, haplotype = row["name"], row[
        "sequence"], row["start"], row["stop"], row["fasta-file"], row["strand"], row["haplotype"]

    # Create a temporary FASTA file
    with NTF(mode='w+', delete=False, suffix='.fasta') as fasta_temp:
        # Including the FASTA file name in the header
        fasta_header = f">{header}:{start}:{stop}:{strand}:{fasta_file_name}:{haplotype}\n"
        fasta_temp.write(fasta_header + sequence + "\n")
        fasta_temp.flush()

        # Create a temporary file for BLAST output
        blast_result_path = NTF(mode='w+', delete=False,
                                suffix='_blast.txt').name
        command = construct_blast_command(
            fasta_temp.name, database_path, identity_cutoff, blast_result_path, len(sequence))
        try:
            subprocess.run(command, shell=True, check=True)
            logger.info(f"Executed BLAST search for {header}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error executing BLAST command for {header}: {e}")
        return blast_result_path


def aggregate_blast_results(dataframe, database_path, CUTOFFS=[100, 75, 50]):
    """
    Iterates over a set of identity cutoffs (100%, 75%, 50%). 
    For every cutoff, it uses parallel searches using 
    "execute_blast_search()" to get the results for every row in the input df. 
    All BLAST results are then saved into a temporary df  
    In the end all temporary df are concatenated to a single df. 
    An additional 'cutoff' column indicating which cutoff is used is also stored. 
    Temporary result files are deleted after their data is read.

    Args:
        dataframe (DataFrame): df with rows representing mappings.
        database_path (Path): Path to the BLAST database.

    Returns:
        pd.DataFrame: df with all the data from the different BLAST searches.
    """
    aggregated_results = pd.DataFrame()
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
                temp_df = pd.read_csv(result_file_path_str,
                                      sep='\t', names=blast_columns)
                if not temp_df.empty:
                    temp_df['cutoff'] = cutoff
                    aggregated_results = pd.concat(
                        [aggregated_results, temp_df], ignore_index=True)
                Path(result_file_path_str).unlink()
    logger.info("Aggregation of BLAST results completed.")
    return aggregated_results


def run_blast_operations(df, db_path, blast_file_path):
    """
    Run BLAST operations and save results to an Excel file.

    Args:
        df (DataFrame): Input DataFrame with mappings.
        db_path (Path): Path to the BLAST database.
        blast_file_path (Path): Path to save the BLAST results Excel file.
    """
    blast_results = aggregate_blast_results(df, db_path)
    blast_results['query cov'] = pd.to_numeric(
        blast_results['query cov'], errors='coerce')
    blast_results = blast_results.query("`query cov` == 100")
    path_df = blast_results['query'].str.split(':', expand=True)
    blast_results[['start', 'stop']] = path_df[[1, 2]]
    blast_results.to_excel(blast_file_path, index=False)
    logger.info("BLAST operations completed and results saved to Excel.")


def blast_main(df, blast_file, library):
    cwd = Path.cwd()
    db_path = make_blast_db(cwd, library)
    run_blast_operations(df, db_path, blast_file)
    logger.info("BLAST main process completed.")
