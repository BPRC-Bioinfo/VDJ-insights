from concurrent.futures import ThreadPoolExecutor
import subprocess
from tempfile import NamedTemporaryFile as NTF
from mapping import mapping_main
from RSS import RSS_main
from write_annotation_report import write_annotation_reports
from pathlib import Path
import pandas as pd


def combine_df():
    """
    Calling the mapping script annotation.py. This script is called 
    with either "bowtie", "bowtie2" or "minimap2" as input value. 
    These all return a df.
    Then all df's are concatinated to form a new complete one. 
    Finally dropping all duplicates in the df based on a
    subset of ["start", "stop"], resetting the index of this df 
    and return this new df.

    Returns:
        unique_combinations (DataFrame): A df containg all the unique
        mapping entries from bowtie(2) and minimap2.
    """
    minimap2_df = mapping_main("minimap2")
    bowtie_df = mapping_main("bowtie")
    bowtie2_df = mapping_main("bowtie2")
    combined = pd.concat([minimap2_df, bowtie_df, bowtie2_df])
    unique_combinations = combined.drop_duplicates(subset=["start", "stop"])
    return unique_combinations.reset_index(drop=True)


def write_report(df, report):
    """
    Creates a excel (xlsx) file of the df and saved as "report.xlsx".

    Args:
        df (DataFrame): df to be saved.
        report (Path): Path of a excel to be saved to.
    """

    df.to_excel(report, index=False)


def get_or_create(annotation_folder):
    """
    Verifies if the report.xlsx is present.
    If present it returns the content of the file as a df. Otherwise
    call the combine_df() function to create the df;
    then save it as report.xlsx and return the df.

    Args:
        annotation_folder (Path): Path of the annotation folder.

    Returns:
        df (DataFrame): df containing information from report.xlsx.
    """
    report = annotation_folder / "report.xlsx"
    if not report.exists():
        df = combine_df()
        write_report(df, report)
        return df
    else:
        return pd.read_excel(report)


def make_dir(dir):
    """
    Create an directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """

    Path(dir).mkdir(parents=True, exist_ok=True)


def make_blast_db(cwd):
    """
    Create a blast database (db) directory folder to run blast with.
    It checks if the db exist, if exist return the db. 
    Otherwise create the db.

    Args:
        cwd (Path): Path of the users current location.

    Returns:
        db (Path): Path of the blast db.
    """
    db = cwd / "blast_db"
    if not db.exists():
        reference = cwd / "library" / "retained.fasta"
        make_dir(db)
        command = f"makeblastdb -in {reference} -dbtype nucl -out {db}/blast_db"
        subprocess.run(command, shell=True)
    return db


def construct_blast_command(fasta_file_path, database_path, identity_cutoff, output_file_path, segment):
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
        segment (str): Segment, either V,D,J. 

    Returns:
        command (str): Constructed blast command. 
    """
    blast_columns = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs"
    command = f"blastn -task megablast -query {fasta_file_path} -db {database_path}/blast_db -outfmt '{blast_columns}' -perc_identity {identity_cutoff} -out {output_file_path}"
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
    header, sequence, start, stop, fasta_file_name, strand = row["name"], row[
        "sequence"], row["start"], row["stop"], row["fasta-file"], row["strand"]

    # Create a temporary FASTA file
    with NTF(mode='w+', delete=False, suffix='.fasta') as fasta_temp:
        # Including the FASTA file name in the header
        fasta_header = f">{header}:{start}:{stop}:{strand}:{fasta_file_name}\n"
        fasta_temp.write(fasta_header + sequence + "\n")
        fasta_temp.flush()

        # Create a temporary file for BLAST output
        blast_result_path = NTF(
            mode='w+', delete=False, suffix='_blast.txt').name
        command = construct_blast_command(
            fasta_temp.name, database_path, identity_cutoff, blast_result_path, row["segment"])
        subprocess.run(command, shell=True)

        return blast_result_path


def aggregate_blast_results(dataframe, database_path):
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
        pd.DataFrame: df with all the datafrom the different BLAST searches.
    """

    aggregated_results = pd.DataFrame()
    with ThreadPoolExecutor() as executor:
        for cutoff in [100, 75, 50]:
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
                temp_df['cutoff'] = cutoff
                aggregated_results = pd.concat(
                    [aggregated_results, temp_df], ignore_index=True)
                Path(result_file_path_str).unlink()

    return aggregated_results


def main():
    """
    Main function for this annotation.py script. 
    It first sets a cwd (current directory the user is in) object 
    based current location. Based on this cwd some 
    input and output directories and files are set. It fetches the 
    needed blast db and the DataFrame (df). The "LOC" are filtered out because
    they gave a distored picture because the lengths are way larger.
    Then it checks if the "blast_results.xlsx" is created. If this is not
    the case the needed df is made. Then the query cov is converted to numeric
    and checked if it is equal to 100%. The start and stop coordinates 
    are filtered out of the query column and also stored in the df.
    This df is saved to "blast_results.xlsx".
    The write_annottion_report() function is called.  
    """
    cwd = Path.cwd()
    annotation_folder = cwd / "annotation"
    make_dir(annotation_folder)
    db = make_blast_db(cwd)
    df = get_or_create(annotation_folder)
    filtered_df = df.query("region != 'LOC'")
    blast_file: Path = annotation_folder / "blast_results.xlsx"
    if not blast_file.exists():
        blast_results = aggregate_blast_results(filtered_df, db)
        blast_results['query cov'] = pd.to_numeric(
            blast_results['query cov'], errors='coerce')
        blast_results = blast_results.query("`query cov` == 100")
        path_df = blast_results['query'].str.split(':', expand=True)
        blast_results[['start', 'stop']] = path_df[[1, 2]]
        blast_results.to_excel(blast_file, index=False)
    else:
        blast_results = pd.read_excel(blast_file)
    write_annotation_reports(annotation_folder)
    RSS_main()


if __name__ == '__main__':
    main()
