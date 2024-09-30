import re
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from Bio import SeqIO

from util import seperate_annotation
from logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)

pd.options.mode.copy_on_write = True


def make_record_dict(fasta):
    """
    Creates a dictionary from a FASTA file, mapping sequence IDs to their corresponding SeqRecord objects.

    Args:
        fasta (str or Path): Path to the input FASTA file.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are SeqRecord objects.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
        IOError: If there is an error reading the FASTA file.
    """
    record_dict = {}
    with open(fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id not in record_dict:
                record_dict[record.id] = record
    return record_dict


def fetch_prefix(name, cell_type):
    """
    Extracts the prefix of a sequence name that starts with the specified cell type.

    Args:
        name (str): The sequence name from which to extract the prefix.
        cell_type (str): The cell type to match against the prefix.

    Returns:
        str: The extracted prefix that starts with the cell type.

    Raises:
        ValueError: If no prefix starting with the cell type is found.
    """
    prefix = [i for i in name.split("_") if i.startswith(cell_type)][0]
    return prefix


def add_region_segment(row, cell_type):
    """
    Determines the region and segment of a sequence based on its name.
    The sequence name is broken into parts using underscores, and the part
    starting with the specified cell type is identified. Numeric values
    are removed from the prefix, and the resulting string is split into the
    region and segment. These values are added as new columns in the DataFrame row.

    Args:
        row (pd.Series): Current row of the DataFrame.
        cell_type (str): The cell type to fetch the prefix.

    Returns:
        pd.Series: The updated row with 'Region', 'Segment', and 'Short name' columns added.
    """
    query = row['Old name-like']
    prefix = fetch_prefix(query, cell_type)
    short_name = prefix
    prefix = re.sub(r"[0-9-]", "", prefix)
    region, segment = prefix[0:3], prefix[3]
    row["Region"], row["Segment"], row["Short name"] = region, segment, short_name
    return row


def filter_group(group):
    """
    Filters out groups of sequences that contain a 100% identity match.

    Args:
        group (pd.DataFrame): A DataFrame group containing sequence alignments.

    Returns:
        bool: True if the group does not contain any 100% identity matches, False otherwise.
    """
    return not (group['% identity'] == 100).any()


def main_df(df):
    """
    Processes the BLAST results DataFrame to filter out entries that contain gaps in the sequences or have 100% identity.
    Converts the '% identity' column to float for numerical operations and filters out 100% identity entries,
    which are considered as non-novel segments. The filtered entries are returned along with a reference DataFrame
    containing only the entries with 100% identity.

    Args:
        df (pd.DataFrame): A DataFrame containing BLAST results.

    Returns:
        tuple: A tuple containing two DataFrames:
            - df (pd.DataFrame): The filtered DataFrame without gaps and 100% identity entries.
            - reference_df (pd.DataFrame): A DataFrame containing only 100% identity entries.
    """
    mask = ~df['query seq'].str.contains('-') & ~df['subject seq'].str.contains('-')
    df = df[mask]
    df['% identity'] = df['% identity'].astype(float)
    reference_df = df.query("`% identity` == 100.000")
    df = df.groupby(['start', 'stop', 'haplotype']).filter(filter_group)
    return df, reference_df


def add_values(df):
    """
    Adds additional columns to the DataFrame, such as the percentage of mismatches relative to the alignment length,
    the lengths of the query and subject sequences, and splits the 'query' column into multiple columns for easier access.

    Args:
        df (pd.DataFrame): A DataFrame containing BLAST results.

    Returns:
        pd.DataFrame: The DataFrame with additional columns for sequence lengths, mismatch percentage, and split query information.
    """
    df['% Mismatches of total alignment'] = (df['mismatches'] / df['alignment length']) * 100
    df['query_seq_length'] = df['query seq'].str.len()
    df['subject_seq_length'] = df['subject seq'].str.len()
    path_df = df['query'].str.split(':', expand=True)
    df[['query', 'start', 'stop', 'strand', 'path', 'haplotype']] = path_df[[0, 1, 2, 3, 4, 5]]
    return df


def add_like_to_df(df):
    """
    Adds a 'Old name-like' column to the DataFrame, which appends '-like' to the reference names.
    This function also sorts the DataFrame by 'Reference' and renames certain columns for consistency.

    Args:
        df (pd.DataFrame): A DataFrame containing BLAST results.

    Returns:
        pd.DataFrame: The DataFrame with the added 'Old name-like' column and renamed columns.
    """
    output_df = df[[
        'subject', 'query',
        'mismatches', '% Mismatches of total alignment',
        'start', 'stop',
        'subject seq', 'query seq',
        'strand', 'path', 'haplotype',
        'query_seq_length', 'subject_seq_length'

    ]]
    output_df.columns = [
        'Reference', 'Old name-like',
        'Mismatches', '% Mismatches of total alignment',
        'Start coord', 'End coord',
        'Reference seq', 'Old name-like seq',
        'Strand', 'Path', 'Haplotype',
        'Reference Length', 'Old name-like Length'
    ]
    output_df = output_df.sort_values(by="Reference")
    output_df['Old name-like'] = output_df['Old name-like'] + '-like'
    return output_df


def orf_function(aa, segment):
    """
    Determines the functional status of a segment based on its amino acid sequence.
    If the sequence ends with a single stop codon ('*') and has no internal stop codons, it is considered 'F/ORF'.
    If the sequence has internal stop codons, it is considered 'P'.
    Additionally, a message is generated to indicate if the stop codon is at a critical position.

    Args:
        aa (str): The amino acid sequence of the segment.
        segment (str): The type of segment (V, D, or J).

    Returns:
        tuple: A tuple containing:
            - message (str): A message indicating the presence of a stop codon.
            - function_type (str): The functional status of the segment ('F/ORF' or 'P').
    """
    end_codon = aa[-1] == "*" and aa.count("*") == 1
    function_dict = {
        "V": "F/ORF" if end_codon or not "*" in aa else "P",
        "D": "PF/ORF",
        "J": "PF/ORF"
    }
    function_type = function_dict.get(segment, "Unknown")
    message = ""
    if end_codon:
        if segment not in function_dict:
            message = "Segment not recognized."
        else:
            message = "STOP-CODON at position 116 (last 3' codon of germline CDR3-IMGT) may disappear during rearrangements"
    return message, function_type


def trim_sequence(sequence, strand):
    """
    Trims a nucleotide sequence to ensure it is a multiple of three bases long,
    and translates it into its corresponding amino acid sequence.
    If the strand is '-', the sequence is reverse complemented before trimming and translation.

    Args:
        sequence (str): The nucleotide sequence to be trimmed and translated.
        strand (str): The strand orientation ('+' or '-').

    Returns:
        tuple: A tuple containing:
            - aa (str): The translated amino acid sequence.
            - sequence (str): The trimmed nucleotide sequence.
    """
    sequence = Seq(sequence)
    seq_length = len(sequence)
    excess_aa = seq_length % 3
    if strand == "-":
        sequence = sequence.reverse_complement()
    if excess_aa != 0:
        sequence = sequence[:-excess_aa]
    aa = sequence.translate()
    return aa, sequence


def add_orf(row):
    """
    Adds a 'Function' column to the DataFrame row, indicating the functional status of the segment.
    The status is determined by translating the nucleotide sequence into amino acids and checking for stop codons.
    If the sequence contains a single stop codon at the end, it is marked as 'F/ORF'. Otherwise, it is marked as 'P'.
    Additionally, a 'Message' column is added to provide information about the stop codon position.

    Args:
        row (pd.Series): The current row of the DataFrame.

    Returns:
        pd.Series: The updated row with 'Function' and 'Message' columns added.
    """
    sequence, strand, segment = row[['Old name-like seq', 'Strand', "Segment"]]
    aa, sequence = trim_sequence(sequence, strand)
    message, function_type = orf_function(aa, segment)
    row["Function"] = function_type
    row["Message"] = message
    return row


def filter_df(row, cell_type):
    """
    Filters a group of sequences to identify the best reference sequence and creates a list of similar references.
    The best reference is determined based on the presence of the 'Specific Part' (a combination of region and segment)
    in the 'Reference' column, sorting by 'Mismatches' and 'Reference', and taking the first hit. 
    If no specific hit is found, the first row of the group is used.
    The remaining similar references are stored in the 'Similar references' column, and the 'Old name-like' is updated.

    Args:
        row (pd.DataFrame): The current group of sequences.
        cell_type (str): The cell type used to fetch the prefix.

    Returns:
        pd.Series: The best reference sequence with an additional column for similar references.
    """
    row['Specific Part'] = row["Region"] + row["Segment"]
    specific_part_in_reference = row.apply(
        lambda x: x['Specific Part'] in x['Reference'], axis=1)
    query_subject_length_equal = row['Reference seq'].str.len(
    ) == row['Old name-like seq'].str.len()
    filtered_rows = row[specific_part_in_reference & query_subject_length_equal]
    best_row = filtered_rows.sort_values(by=['Mismatches', 'Reference']).head(1)

    if best_row.empty:
        best_row = row.head(1)
    all_references = ', '.join(
        set(row['Reference']) - set(best_row['Reference']))
    best_row['Similar references'] = all_references
    best_row["Old name-like"] = best_row["Reference"]
    best_row["Short name"] = best_row["Reference"].apply(
        lambda ref: fetch_prefix(ref, cell_type))
    if int(best_row.Mismatches.iloc[0]) != 0:
        best_row["Old name-like"] = best_row["Reference"] + "-like"
        best_row["Short name"] = best_row["Reference"].apply(
            lambda ref: fetch_prefix(ref, cell_type)) + "-like"
    return best_row.squeeze()


def add_reference_length(row, record):
    """
    Adds the length of the reference sequence to the DataFrame row.
    The length is determined by looking up the sequence in the record dictionary.

    Args:
        row (pd.Series): The current row of the DataFrame.
        record (dict): A dictionary mapping reference names to SeqRecord objects.

    Returns:
        pd.Series: The updated row with the 'Library Length' column added.
    """
    reference = row["Reference"]
    row["Library Length"] = len(record[reference].seq)
    return row


def extract_sample(path):
    filename = path.split("/")[-1]
    sample_pattern = re.compile(r'(GCA|GCF|DRR|ERR)_?\d{6,9}(\.\d+)?')
    match = sample_pattern.search(filename)
    if match:
        return match.group(0)
    else:
        return filename.split("_")[0]


def run_like_and_length(df, record, cell_type):
    """
    Adds 'Old name-like', 'Region', 'Segment', and 'Library Length' columns to the DataFrame.
    Filters the DataFrame to retain only rows where the reference length, old name-like length, and library length are equal.
    Adds a 'Sample' column, which is extracted from the path.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        record (dict): A dictionary mapping reference names to SeqRecord objects.
        cell_type (str): The cell type used to fetch the prefix.

    Returns:
        pd.DataFrame: The processed DataFrame with additional columns and filtered rows.
    """
    df = add_like_to_df(df)
    df = df.apply(add_region_segment, axis=1, cell_type=cell_type)
    df = df.apply(add_reference_length, axis=1, record=record)
    length_mask = df[["Reference Length", "Old name-like Length", "Library Length"]].apply(lambda x: x.nunique() == 1, axis=1)
    df = df[length_mask]
    df["Sample"] = df["Path"].apply(extract_sample)
    return df


def group_similar(df, cell_type):
    """
    Groups sequences by their start coordinate, end coordinate, and haplotype, and filters each group to identify the best reference sequence.
    The best reference sequence is selected based on the number of mismatches and the reference name, and similar references are stored.

    Args:
        df (pd.DataFrame): The DataFrame to group and filter.
        cell_type (str): The cell type used to fetch the prefix.

    Returns:
        pd.DataFrame: The grouped and filtered DataFrame.
    """
    df = df.groupby(['Start coord', 'End coord', 'Haplotype']).apply(lambda group: filter_df(group, cell_type))
    df = df.reset_index(drop=True)
    return df


def annotation_long(df, annotation_folder):
    """
    Generates a condensed version of the annotation report, saving it as 'annotation_report_long.xlsx'.
    The report includes key columns such as reference names, coordinates, functions, paths, and regions.

    Args:
        df (pd.DataFrame): The DataFrame containing BLAST results.
        annotation_folder (Path): The directory where the report will be saved.

    Raises:
        OSError: If the file cannot be created or written to.
    """
    console_log.info("Generating annotation_report_long.xlsx!")
    df = df[['Reference', 'Old name-like', 'Mismatches',
             '% Mismatches of total alignment', 'Start coord',
             'End coord', 'Function', 'Path', 'Region', 'Segment',
             'Haplotype', 'Sample', 'Short name', 'Message']]
    df.to_excel(annotation_folder / 'annotation_report_long.xlsx', index=False)


def annotation(df: pd.DataFrame, annotation_folder, file_name, no_split):
    """
    Generates a full annotation report and saves it as the specified file name.
    The report includes key columns such as reference names, coordinates, functions, similar references, paths, and regions.

    Args:
        df (pd.DataFrame): The DataFrame containing BLAST results.
        annotation_folder (Path): The directory where the report will be saved.
        file_name (str): The name of the output Excel file.

    Raises:
        OSError: If the file cannot be created or written to.
    """
    console_log.info(f"Generating {file_name}!")
    df = df[['Reference', 'Old name-like', 'Mismatches',
             '% Mismatches of total alignment', 'Start coord',
             'End coord', 'Function', 'Similar references', 'Path',
             'Strand', 'Region', 'Segment', 'Haplotype', 'Sample',
             'Short name', 'Message', 'Old name-like seq', 'Reference seq',]]
    df["Status"] = "Known" if "known" in file_name else "Novel"

    full_annotation_path = annotation_folder / file_name
    df.to_excel(full_annotation_path, index=False)

    if not no_split:
        console_log.info("Creating individual sample excel files...")
        df.groupby("Sample").apply(lambda group: seperate_annotation(group, annotation_folder, file_name))


def report_main(annotation_folder: str | Path, blast_file: str | Path, cell_type: str, library: str | Path, no_split: bool):
    """
    Main function to process and generate the annotation reports from the BLAST results.
    It performs the following steps:
        1. Loads the configuration and reference sequence records.
        2. Processes the BLAST results DataFrame to add necessary columns and filter the results.
        3. Generates the condensed and full annotation reports.

    Args:
        annotation_folder (Path): The directory where the reports will be saved.
        blast_file (Path or str): Path to the input BLAST results file (Excel format).
        cell_type (str): The cell type used to fetch prefixes.
        library (Path or str): Path to the reference sequence library in FASTA format.

    Raises:
        Exception: If any step fails, logs the error and raises an exception.
    """
    cwd = Path.cwd()
    record = make_record_dict(cwd / library)
    df = pd.read_csv(blast_file)
    df = add_values(df)
    df, ref_df = main_df(df)
    df = run_like_and_length(df, record, cell_type)
    ref_df = run_like_and_length(ref_df, record, cell_type)
    df, ref_df = df.apply(add_orf, axis=1), ref_df.apply(add_orf, axis=1)
    annotation_long(df, annotation_folder)
    df, ref_df = group_similar(df, cell_type), group_similar(ref_df, cell_type)

    annotation(df, annotation_folder, 'annotation_report_novel.xlsx', no_split)
    annotation(ref_df, annotation_folder,'annotation_report_known.xlsx', no_split)
