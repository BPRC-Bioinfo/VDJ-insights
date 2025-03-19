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

from .util import log_error, make_dir
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
    query = row['Target name']
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


def split_df(df):
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
                insertions += 1
            elif i + 1 < len(btop) and btop[i + 1].isalpha():
                deletions += 1
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
    df[['query', 'start', 'stop', 'strand', 'path', 'tool']] = split_query_df[[0, 1, 2, 3, 4, 5]]
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
        'subject', 'query','mismatches', '% Mismatches of total alignment',
        'start', 'stop', 'subject seq', 'query seq', 'tool', '% identity', 'strand', 'path',
        'query_seq_length', 'subject_seq_length', 'btop', 'SNPs', 'Insertions', 'Deletions'

    ]]
    output_df.columns = [
        'Library name', 'Target name', 'Mismatches', '% Mismatches of total alignment',
        'Start coord', 'End coord', 'Library sequence', 'Target sequence', 'Mapping tool', '% identity','Strand', 'Path',
        'Target length', 'subject library length', 'BTOP', 'SNPs', 'Insertions', 'Deletions'
    ]

    output_df = output_df.sort_values(by="Library name")
    return output_df



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


def filter_df(group_df, cell_type):
    """
    Filters a group of sequences to identify the best reference sequence and creates a list of similar references.
    The best reference is determined based on the presence of the 'Specific Part' (a combination of region and segment)
    in the 'Reference' column, sorting by 'Mismatches' and 'Reference', and taking the first hit.
    If no specific hit is found, the first row of the group is used.
    The remaining similar references are stored in the 'Similar references' column, and the 'Target name' is updated.

    Args:
        group_df (pd.DataFrame): The current group of sequences.
        cell_type (str): The cell type used to fetch the prefix.

    Returns:
        pd.Series: The best reference sequence with an additional column for similar references.
    """
    group_df['Specific Part'] = group_df["Region"] + group_df["Segment"]
    specific_part_in_reference = group_df.apply(
        lambda x: x['Specific Part'] in x['Library name'], axis=1)
    query_subject_length_equal = group_df['Library sequence'].str.len(
    ) == group_df['Target sequence'].str.len()
    filtered_rows = group_df[specific_part_in_reference & query_subject_length_equal]
    best_row = filtered_rows.sort_values(by=['Mismatches', 'Library name']).head(1)

    if best_row.empty:
        best_row = group_df.head(1)
    all_references = ', '.join(
        set(group_df['Library name']) - set(best_row['Library name']))
    best_row['Similar references'] = all_references
    best_row["Target name"] = best_row["Library name"]
    best_row["Short name"] = best_row["Library name"].apply(
        lambda ref: fetch_prefix(ref, cell_type))
    if float(best_row["% identity"].iloc[0]) < 100.0:
        best_row["Target name"] = best_row["Library name"] + "-like"
        best_row["Short name"] = best_row["Library name"].apply(
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
    reference = row["Library name"]
    row["Library length"] = len(record[reference].seq)
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
    Adds 'Target name', 'Region', 'Segment', and 'Library Length' columns to the DataFrame.
    Filters the DataFrame to retain only rows where the reference length, Target name length, and library length are equal.
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
    df["Sample"] = df["Path"].apply(extract_sample)
    return df


def group_similar(df, cell_type):
    """
    Groups sequences by their start coordinate, end coordinate, and filters each group to identify the best reference sequence.
    The best reference sequence is selected based on the number of mismatches and the reference name, and similar references are stored.

    Args:
        df (pd.DataFrame): The DataFrame to group and filter.
        cell_type (str): The cell type used to fetch the prefix.

    Returns:
        pd.DataFrame: The grouped and filtered DataFrame.
    """
    df = df.groupby(['Sample', 'Start coord', 'End coord']).apply(lambda group: filter_df(group, cell_type))
    df = df.reset_index(drop=True)

    cols_to_convert = ["Start coord", "End coord", "% identity", "Mismatches"]
    df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors='coerce')


    filterd_df = (
        df
        .sort_values(by=['Sample', '% identity', 'Mismatches'], ascending=[True, False, True])
        .groupby(['Start coord', 'End coord'], as_index=False)
        .first()
    )

    filterd_df["Alignment_length"] = filterd_df["End coord"] - filterd_df["Start coord"]
    longest_sequences = (
        filterd_df
        .sort_values(by=["Sample", "% identity", "Alignment_length"], ascending=[True, False, False])
        .groupby(["Sample", "Start coord"], as_index=False)
        .first()
    )
    processed_groups = []
    grouped = longest_sequences.groupby(["Sample", "Region"], as_index=False)
    for name, group in grouped:
        group = group.sort_values(by=["Start coord", "% identity"], ascending=[True, False]).reset_index(drop=True)
        merged_intervals = []
        current_interval = group.iloc[0]

        for idx in range(1, len(group)):
            next_interval = group.iloc[idx]
            if next_interval["Start coord"] <= current_interval["End coord"]:
                if next_interval["% identity"] > current_interval["% identity"]:
                    current_interval = next_interval
                else:
                    group.loc[current_interval.name, "End coord"] = max(current_interval["End coord"], next_interval["End coord"])
            else:
                merged_intervals.append(current_interval)
                current_interval = next_interval

        merged_intervals.append(current_interval)

        merged_group = pd.DataFrame(merged_intervals)
        processed_groups.append(merged_group)

    result_df = pd.concat(processed_groups, ignore_index=True)
    return result_df


def annotation(df: pd.DataFrame, annotation_folder, file_name, metadata_folder):
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
    df = df[["Sample", "Region", "Segment", "Start coord", "End coord", "Strand", "Target name", "Library name", "Short name", "Similar references", "Target sequence", "Library sequence", "Mismatches", "% Mismatches of total alignment", "% identity", "BTOP", "SNPs", "Insertions", "Deletions", "Mapping tool", "Status", "Path"]]

    if metadata_folder:
        metadata_df = pd.read_excel(metadata_folder)
        merge_cols = ["Accession"]
        if "Population" in metadata_df.columns:
            merge_cols.append("Population")
        if "Haplotype" in metadata_df.columns:
            merge_cols.append("Haplotype")
        if len(merge_cols) > 1:
            df = df.merge(metadata_df[merge_cols], left_on="Sample", right_on="Accession", how="left")
            df = df.drop(columns=["Accession"])

    full_annotation_path = annotation_folder / file_name
    df.to_excel(full_annotation_path, index=False)

    for sample, sample_df in df.groupby("Sample"):
        path = annotation_folder / "individual" / sample
        make_dir(path)
        new_file = path / file_name
        sample_df.to_excel(new_file, index=False)


def add_suffix_to_short_name(group):
    """
    Adds a suffix to duplicate 'Short name' entries to ensure uniqueness.

    Args:
        group (pd.DataFrame): Grouped DataFrame containing segment data.

    Returns:
        pd.DataFrame: Updated DataFrame with unique 'Short name' values.
    """
    unique_sequences = group['Target sequence'].unique()
    if len(unique_sequences) > 1:
        seq_to_suffix = {seq: f"_{i + 1}" for i, seq in enumerate(unique_sequences)}
        group['Short name'] = group['Short name'] + group['Target sequence'].map(seq_to_suffix)
        group['Target name'] = group['Target name'] + group['Target sequence'].map(seq_to_suffix)
    return group


@log_error()
def report_main(annotation_folder: Union[str, Path], blast_file: Union[str, Path], cell_type: str, library: Union[str, Path], metadata_folder: Union[str, Path]):
    """
    Main function to process and generate the annotation reports from the BLAST results.
    """
    segments_library = make_record_dict(library)

    df = pd.read_csv(blast_file, low_memory=False)
    df = add_values(df)

    novel_df, known_df = split_df(df)

    if not novel_df.empty:
        novel_df = run_like_and_length(novel_df, segments_library, cell_type)
        novel_df["Status"] = "Novel"
    if not known_df.empty:
        known_df = run_like_and_length(known_df, segments_library, cell_type)
        known_df["Status"] = "Known"

    combined_df = pd.concat([novel_df, known_df], ignore_index=True)

    grouped_df = group_similar(combined_df, cell_type)

    known_df = grouped_df[grouped_df["Status"] == "Known"]
    novel_df = grouped_df[grouped_df["Status"] == "Novel"]

    if not novel_df.empty:
        novel_df = novel_df.groupby('Library name', group_keys=False).apply(add_suffix_to_short_name)
        annotation(novel_df, annotation_folder, 'annotation_report_novel.xlsx', metadata_folder)

    if not known_df.empty:
        annotation(known_df, annotation_folder, 'annotation_report_known.xlsx', metadata_folder)

    combined_df = pd.concat([novel_df, known_df], ignore_index=True)
    if not combined_df.empty:
        annotation(combined_df, annotation_folder, 'annotation_report_all.xlsx', metadata_folder)

    console_log.info(f"Known segments detected: {known_df.shape[0]}")
    console_log.info(f"Novel segments detected: {novel_df.shape[0]}")


    #pivot_table_sheet1 = combined_df.pivot_table(index='Sample', columns=['Region', 'Status', 'Function', 'Segment'], aggfunc='size', fill_value=0)
    pivot_table_sheet2 = combined_df.pivot_table(index='Sample', columns=['Region', 'Segment'], aggfunc='size', fill_value=0)
    #pivot_table_sheet3 = combined_df.pivot_table(index='Sample', columns=['Region', 'Function', 'Segment'], aggfunc='size', fill_value=0)
    pivot_table_sheet4 = combined_df.pivot_table(index='Sample', columns=['Region', 'Status', 'Segment'], aggfunc='size', fill_value=0)
    pivot_table_sheet5 = combined_df.pivot_table(index='Sample', columns=['Region'], aggfunc='size', fill_value=0)
    pivot_table_sheet6 = combined_df.pivot_table(index='Sample', columns=['Region', 'Status'], aggfunc='size', fill_value=0)

    with pd.ExcelWriter(f"{annotation_folder}/pivot_tables.xlsx") as writer:
        #pivot_table_sheet1.to_excel(writer, sheet_name="Total")
        pivot_table_sheet2.to_excel(writer, sheet_name="Segment")
        #pivot_table_sheet3.to_excel(writer, sheet_name="Function")
        pivot_table_sheet4.to_excel(writer, sheet_name="Status")
        pivot_table_sheet5.to_excel(writer, sheet_name="Region_total")
        pivot_table_sheet6.to_excel(writer, sheet_name="Function_total")
