"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

import re
from typing import Union

import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from Bio import SeqIO

from .util import seperate_annotation, log_error
from .logger import console_logger, file_logger

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
    df['% identity'] = df['% identity'].astype(float)
    reference_df = df[df['% identity'] == 100.0]

    df = df.groupby(['start', 'stop', 'path']).filter(filter_group)
    return df, reference_df


def parse_btop(btop):
    snps = 0
    insertions = 0
    deletions = 0
    i = 0
    while i < len(btop):
        if btop[i].isdigit():
            while i < len(btop) and btop[i].isdigit():
                i += 1
        elif i + 1 < len(btop) and btop[i].isalpha() and btop[i + 1].isalpha():
            snps += 1
            i += 2
        elif btop[i] == "-":
            if i > 0 and btop[i - 1].isalpha():
                deletions += 1
            elif i + 1 < len(btop) and btop[i + 1].isalpha():
                insertions += 1
            i += 1
        else:
            i += 1
    return snps, insertions, deletions


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
    df[['SNPs', 'Insertions', 'Deletions']] = df['btop'].apply(lambda x: pd.Series(parse_btop(x)))

    split_query_df = df['query'].str.split('#', expand=True)
    #split_query_df = df['query'].str.split(r'[#:]', expand=True)

    df[['query', 'start', 'stop', 'strand', 'path', 'haplotype', 'tool', 'mapping_accuracy']] = split_query_df[[0, 1, 2, 3, 4, 5, 6, 7]]
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
        'tool', 'mapping_accuracy', '% identity', 'strand', 'path', 'haplotype',
        'query_seq_length', 'subject_seq_length', 'btop', 'SNPs', 'Insertions', 'Deletions'

    ]]
    output_df.columns = [
        'Reference', 'Old name-like',
        'Mismatches', '% Mismatches of total alignment',
        'Start coord', 'End coord',
        'Reference seq', 'Old name-like seq',
        'tool', 'mapping_accuracy', '% identity','Strand', 'Path', 'Haplotype',
        'Reference Length', 'Old name-like Length', 'BTOP', 'SNPs', 'Insertions', 'Deletions'
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
        else: #ALLEEN VOOR V SEGMENTEN
            message = "the STOP-CODON at the 3' end of the V-REGION can be deleted by rearrangement"
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
    excess_bases = seq_length % 3
    if strand == "-":
        sequence = sequence.reverse_complement()
    if excess_bases != 0:
        sequence = sequence[:-excess_bases]
    amino_acid_sequence = sequence.translate()
    return amino_acid_sequence, sequence


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
    sequence = sequence.replace("-", "")
    amino_acid_sequence, sequence = trim_sequence(sequence, strand)
    message, function_type = orf_function(amino_acid_sequence, segment)
    row["Function"] = function_type
    row["Message"] = message
    return row


def filter_df(group_df, cell_type):
    """
    Filters a group of sequences to identify the best reference sequence and creates a list of similar references.
    The best reference is determined based on the presence of the 'Specific Part' (a combination of region and segment)
    in the 'Reference' column, sorting by 'Mismatches' and 'Reference', and taking the first hit.
    If no specific hit is found, the first row of the group is used.
    The remaining similar references are stored in the 'Similar references' column, and the 'Old name-like' is updated.

    Args:
        group_df (pd.DataFrame): The current group of sequences.
        cell_type (str): The cell type used to fetch the prefix.

    Returns:
        pd.Series: The best reference sequence with an additional column for similar references.
    """
    group_df['Specific Part'] = group_df["Region"] + group_df["Segment"]
    specific_part_in_reference = group_df.apply(
        lambda x: x['Specific Part'] in x['Reference'], axis=1)
    query_subject_length_equal = group_df['Reference seq'].str.len(
    ) == group_df['Old name-like seq'].str.len()
    filtered_rows = group_df[specific_part_in_reference & query_subject_length_equal]
    best_row = filtered_rows.sort_values(by=['Mismatches', 'Reference']).head(1)

    if best_row.empty:
        best_row = group_df.head(1)
    all_references = ', '.join(
        set(group_df['Reference']) - set(best_row['Reference']))
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

    #length_mask = df[["Reference Length", "Old name-like Length", "Library Length"]].apply(lambda x: x.nunique() == 1, axis=1)
    #df = df[length_mask]
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
    df = df.groupby(['Start coord', 'End coord', 'Sample']).apply(lambda group: filter_df(group, cell_type))
    df = df.reset_index(drop=True)

    cols_to_convert = ["Start coord", "End coord", "% identity", "mapping_accuracy", "Mismatches"]
    df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors='coerce')

    filterd_df = (
        df
        .sort_values(
            by=['% identity', 'mapping_accuracy', 'Mismatches'],
            ascending=[False, False, True]
        )
        .groupby(['Start coord', 'End coord'], as_index=False)
        .first()
    )
    filterd_df["Alignment_length"] = filterd_df["End coord"] - filterd_df["Start coord"]
    longest_sequences = (
        filterd_df
        .sort_values(["% identity", "Alignment_length"], ascending=[False, False])
        .groupby("Start coord")
        .head(1)
    )

    processed_groups = []
    grouped = longest_sequences.groupby(["Sample", "Region", "Segment"], as_index=False)
    for name, group in grouped:
        group = group.sort_values(by="Start coord").reset_index(drop=True)

        merged_intervals = []
        current_interval = group.iloc[0]

        for idx in range(1, len(group)):
            next_interval = group.iloc[idx]
            if next_interval["Start coord"] <= current_interval["End coord"]:
                group.loc[current_interval.name, "End coord"] = max(current_interval["End coord"],
                                                                    next_interval["End coord"])
            else:
                merged_intervals.append(current_interval)
                current_interval = next_interval

        merged_intervals.append(current_interval)

        merged_group = pd.DataFrame(merged_intervals)
        processed_groups.append(merged_group)

    result_df = pd.concat(processed_groups, ignore_index=True)
    return result_df


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
    file_log.info("Generating annotation_report_long.xlsx!")
    df = df[['Reference', 'Old name-like', 'Mismatches',
             '% Mismatches of total alignment', 'Start coord',
             'End coord', 'Function', 'Path', 'Region', 'Segment',
             'Haplotype', 'Sample', 'Short name', 'Message', 'tool','mapping_accuracy', '% identity']]
    df.to_excel(annotation_folder / 'annotation_report_long.xlsx', index=False)


def annotation(df: pd.DataFrame, annotation_folder, file_name, no_split, metadata_folder):
    """
    Generates a full annotation report and saves it as the specified file name.
    The report includes key columns such as reference names, coordinates, functions, similar references, paths, and regions.

    Args:
        df (pd.DataFrame): The DataFrame containing BLAST results.
        annotation_folder (Path): The directory where the report will be saved.
        file_name (str): The name of the output Excel file.
        no_split (bool): A flag indicating whether to split the report by sample.
        metadata_folder (Path): The path to the metadata file.

    Raises:
        OSError: If the file cannot be created or written to.
    """
    file_log.info(f"Generating {file_name}!")

    df["Status"] = "Known" if "known" in file_name else "Novel"
    df = df[["Sample", "Haplotype", "Region", "Segment", "Start coord", "End coord", "Strand", "Reference", "Old name-like", "Short name", "Similar references", "Old name-like seq", "Reference seq", "Mismatches", "% Mismatches of total alignment", "% identity", "BTOP", "SNPs", "Insertions", "Deletions", "mapping_accuracy", "tool", "Function", "Status", "Message", "Path"]]

    if metadata_folder:
        metadata_df = pd.read_excel(metadata_folder)
        df = df.merge(metadata_df[['Accession', 'Population']], left_on='Sample', right_on='Accession', how='left')
        df = df.drop(columns=["Accession"])

    full_annotation_path = annotation_folder / file_name
    df.to_excel(full_annotation_path, index=False)

    if not no_split:
        file_log.info("Creating individual sample excel files...")
        df.groupby("Sample").apply(lambda group: seperate_annotation(group, annotation_folder, file_name))


@log_error()
def report_main(annotation_folder: Union[str, Path], blast_file: Union[str, Path], cell_type: str, library: Union[str, Path], no_split: bool, metadata_folder: Union[str, Path]):
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
        no_split (bool): A flag indicating whether to split the report by sample.
        metadata_folder (Path or str): Path to the metadata file.

    Raises:
        Exception: If any step fails, logs the error and raises an exception.
    """
    cwd = Path.cwd()
    segments_library = make_record_dict(library)

    df = pd.read_csv(blast_file, low_memory=False)
    df = add_values(df)

    novel_df, known_df = main_df(df)
    if not novel_df.empty:
        novel_df = run_like_and_length(novel_df, segments_library, cell_type)
        novel_df = novel_df.apply(add_orf, axis=1)
        annotation_long(novel_df, annotation_folder)
        novel_df = group_similar(novel_df, cell_type)
        annotation(novel_df, annotation_folder, 'annotation_report_novel.xlsx', no_split, metadata_folder)
    if not known_df.empty:
        known_df = run_like_and_length(known_df, segments_library, cell_type)
        known_df = known_df.apply(add_orf, axis=1)
        known_df = group_similar(known_df, cell_type)
        annotation(known_df, annotation_folder,'annotation_report_known.xlsx', no_split, metadata_folder)

