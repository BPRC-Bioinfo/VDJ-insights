import re
import yaml
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from Bio import SeqIO
import sys
from logger import custom_logger

# Method for logger current states of the program.
logger = custom_logger(__name__)

CONFIG = None
pd.options.mode.copy_on_write = True


def load_config(cwd):
    global CONFIG
    config_file = Path(cwd / 'config' / 'config.yaml')
    if config_file.exists():
        with open(config_file, 'r') as file:
            CONFIG = yaml.safe_load(file)
    else:
        logger.error("No configuration file provided, closing application!")
        sys.exit()


def make_record_dict(fasta):
    with open(fasta, 'r') as fasta_file:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        return record_dict


def fetch_prefix(name):
    cell_type = CONFIG.get("SPECIES", {}).get("cell", "TR")
    prefix = [i for i in name.split("_") if i.startswith(cell_type)][0]
    return prefix


def add_region_segment(row):
    """
    Determines based on the potential name of the segment what the 
    segment and the region is. It first brakes the name based on the "_" 
    in to a list and check if it finds a a part that starts with a value 
    specified in the option variable. When a part if found all the 
    numeric values are removed with regex.
    The two values are appointed to the columns "Region" and "Segment" 
    and the new row is returned.

    Args:
        row (Series): Current row of the df.

    Returns:
        row (series): Current row with the two extra columns 
        "Region" and "Segment".
    """
    query = row['Old name-like']
    prefix = fetch_prefix(query)
    short_name = prefix
    prefix = re.sub(r"[0-9-]", "", prefix)
    region, segment = prefix[0:3], prefix[3]
    row["Region"], row["Segment"], row["Short name"] = region, segment, short_name
    return row


def filter_group(group):
    """
    Check if group contains a 100% identity. If this is the case,
    the group is discarded.

    Args:
        group (Series): Dataframe grouped by ["start", "stop"].

    Returns:
        bool: Boolean that indicates if a group contains a 100% identity
        (True otherwise False).
    """
    return not (group['% identity'] == 100).any()


def main_df(df):
    """
    Processes the BLAST results df and filtering out entries.
    It filters out entries where either 'query seq' or 'subject seq'
    contains a "-", as these represent gaps in the alignment.
    Also converts the '% identity' column to float for easier numerical 
    operations. Additionally, filter out entries with 100% identity, 
    as these do not contain new segments. 
    These entries are filtered and saved in the reference_df.

    Args:
        df (pd.DataFrame): A df containing BLAST results.
    Returns:
        df (DataFrame): The filtered DataFrame with entries having no 
        dashes in sequences and excluding 100% identity entries.
        reference_df (DataFrame): A reference df containing only
        the entries with 100% identity.
    """
    mask = ~df['query seq'].str.contains(
        '-') & ~df['subject seq'].str.contains('-')
    df = df[mask]
    df['% identity'] = df['% identity'].astype(float)
    reference_df = df.query("`% identity` == 100.000")
    df = df.groupby(['start', 'stop', 'haplotype']).filter(filter_group)
    return df, reference_df


def add_values(df):
    """
    Add different columns to the df such as "% Mismatches of total alignment". 
    This is an indicator of the percentage of the amount mismatches in a sequence.
    Also add the length of the query (query_seq_length) and
    subject (subject_seq_length) to the df. Finally split (":") the 
    query in to 5 columns ['query', 'start', 'stop', 'strand', 'path'] and
    add them to to df.


    Args:
        df (DataFrame): An df containing BLAST results.

    Returns:
        df (DataFrame): BLAST result df with the extra columns added.
    """
    df['% Mismatches of total alignment'] = (
        df['mismatches'] / df['alignment length']) * 100
    df['query_seq_length'] = df['query seq'].str.len()
    df['subject_seq_length'] = df['subject seq'].str.len()
    path_df = df['query'].str.split(':', expand=True)
    df[['query', 'start', 'stop', 'strand', 'path',
        'haplotype']] = path_df[[0, 1, 2, 3, 4, 5]]
    return df


def add_like_to_df(df):
    """
    Sorts the BLAST results df by 'reference'. It also creates a new 
    column 'Old name-like' which contains the reference with -like 
    behind it. It also extract the important columns and renames them. 

    Args:
        df (DataFrame): An df containing BLAST results.

    Returns:
        output_df (DataFrame): df with added "Old name-like" column.
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
    Check if the sequence contains a * when translated to amino acids. 
    When the strand is "-" the reverse_complete() is taken from the 
    sequence and also checked. A column named "Function" is created. 
    When the sequence contains a "*" a "P" is stored in the 
    "Function" otherwise a "F/ORF".

    Args:
        row (Series): An row from the BLAST df.

    Returns:
        row (Series): The row with the extra "Function" column.
    """
    sequence, strand, segment = row[['Old name-like seq', 'Strand', "Segment"]]
    aa, sequence = trim_sequence(sequence, strand)
    message, function_type = orf_function(aa, segment)
    row["Function"] = function_type
    row["Message"] = message
    return row


def filter_df(row):
    """
    Filters the BLAST row/section to identify the best reference and creates 
    a list of leftover similar references. First the Specific Part,
    which begins with the chosen cell type. It is determined from the 
    Old name-like column. Then selects the best reference based on the presence of 
    this 'Specific Part' in the 'Reference' column, sorting by 
    'Mismatches' and 'Reference' and taking the first hit. 
    If no specific hit is found, the first row of the df is 
    used. Remaining similar references are stored in the 'Similar 
    references' column. In the end the Old name-like naming is 
    constructed again.

    Args:
        row (pd.Series): The current section that is being filtered.

    Returns:
        best_row: the BLAST df with now the best reference chosen and 
        a extra column "Similar references" containing the list with similar 
        references.
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
        fetch_prefix)
    if int(best_row.Mismatches.iloc[0]) != 0:
        best_row["Old name-like"] = best_row["Reference"] + "-like"
        best_row["Short name"] = best_row["Reference"].apply(
            fetch_prefix) + "-like"
    return best_row.squeeze()


def add_reference_length(row, record):
    reference = row["Reference"]
    row["Library Length"] = len(record[reference].seq)
    return row


def run_like_and_length(df, record):
    df = add_like_to_df(df)
    df = df.apply(add_region_segment, axis=1)
    df = df.apply(add_reference_length, axis=1, record=record)
    length_mask = df[["Reference Length",	"Old name-like Length",
                      "Library Length"]].apply(lambda x: x.nunique() == 1, axis=1)
    df = df[length_mask]
    df["Sample"] = df["Path"].str.split("/").str[-1].str.split("_").str[0]
    return df


def group_similar(df):
    df = df.groupby(['Start coord', 'End coord', 'Haplotype']).apply(
        filter_df)
    df = df.reset_index(drop=True)
    return df


def annotation_long(df, annotation_folder):
    """
    Generates a condensed df for the 
    "annotation_report_long.xlsx" file, 
    selecting only the specific columns needed for
    "annotation_report_long.xlsx".

    Args:
        df (DataFrame): An df containing BLAST results.
        annotation_folder (Path): Path to the annotation_report_long.xlsx file 
    """
    logger.info("Generating annotation_report_long.xlsx!")
    df = df[['Reference', 'Old name-like', 'Mismatches',
             '% Mismatches of total alignment', 'Start coord',
             'End coord', 'Function', 'Path', 'Region', 'Segment',
             'Haplotype', 'Sample', 'Short name', 'Message']]
    df.to_excel(annotation_folder / 'annotation_report_long.xlsx', index=False)


def annotation(df, annotation_folder, file_name):
    """
    Generates a condensed df for the 
    "annotation_report.xlsx" file, 
    selecting only the specific columns needed for
    "annotation_report.xlsx".

    Args:
        df (DataFrame): An df containing BLAST results.
        annotation_folder (Path): Path to the annotation_report.xlsx file.
        file_name (str): Name of the excel file the df is written to. 
    """
    # Loose columns 'Reference seq', 'Old name-like seq', 'Reference Length',
    # 'Old name-like Length', 'Library Length',

    logger.info(f"Generating {file_name}!")
    df = df[['Reference', 'Old name-like', 'Mismatches',
             '% Mismatches of total alignment', 'Start coord',
             'End coord', 'Function', 'Similar references', 'Path',
             'Strand', 'Region', 'Segment', 'Haplotype', 'Sample',
             'Short name', 'Message']]
    df["Status"] = "Known" if "100%" in file_name else "Novel"
    df.to_excel(annotation_folder / file_name, index=False)


def write_annotation_reports(annotation_folder):
    """
    Main function for generating reports from BLAST results. 
    It reads the initial BLAST DataFrame (df), processes it through various filtering 
    and data manipulation steps, and generates multiple Excel reports.

    Reads the initial BLAST results from 'blast_results.xlsx'. 
    Modify columns in the df. Filters the df to retain only 
    the best rows based on specific criteria. In the end get 
    different subsets of the df and saves them into various Excel files,
    which include annotated long report, a standard annotation report, 
    and an RSS report.


    Args:
        annotation_folder (Path): Path of the initial annotation folder.
    """
    cwd = Path.cwd()
    load_config(cwd)
    record = make_record_dict(cwd / "library" / "library.fasta")
    df = pd.read_excel(annotation_folder / "blast_results.xlsx")
    df = add_values(df)
    df, ref_df = main_df(df)
    df, ref_df = run_like_and_length(
        df, record), run_like_and_length(ref_df, record)
    df, ref_df = df.apply(add_orf, axis=1), ref_df.apply(add_orf, axis=1)
    annotation_long(df, annotation_folder)
    df, ref_df = group_similar(df), group_similar(ref_df)
    annotation(df, annotation_folder, 'annotation_report.xlsx')
    annotation(ref_df, annotation_folder, 'annotation_report_100%.xlsx')


if __name__ == '__main__':
    cwd = Path.cwd()
    annotation_folder = cwd / "annotation"
    write_annotation_reports(annotation_folder)
