import re
from typing import Union

import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from Bio import SeqIO

from .util import log_error
from .logger import console_logger, file_logger

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


def filter_group(group):
    """
    Filters out groups of sequences that contain a 100% identity match.

    Args:
        group (pd.DataFrame): A DataFrame group containing sequence alignments.

    Returns:
        bool: True if the group does not contain any 100% identity matches, False otherwise.
    """
    return not (group['% identity'] == 100).any()


def select_best_regions(df):
    """
    Select the best regions when overlaps exist, including nested overlaps.

    Args:
        df (pd.DataFrame): DataFrame containing columns: Start, Stop, %Identity, and Mismatches.

    Returns:
        pd.DataFrame: Filtered DataFrame with only the best regions selected.
    """
    cols_to_convert = ['start', 'stop', '% identity', 'cutoff', 'mapping_accuracy', 'mismatches', 'evalue']
    df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors='coerce')

    filterd_df = (
        df
        .sort_values(
            by=['% identity', 'cutoff', 'mapping_accuracy', 'mismatches'],
            ascending=[False, False, False, True]
        )
        .groupby(['start', 'stop'], as_index=False)
        .first()
    )

    longest_sequences = (
        filterd_df
        .sort_values([ "% identity", "alignment length"], ascending=[False, False])
        .groupby("start")
        .head(1)
    )
    return longest_sequences


def extract_sample(path):
    filename = path.split("/")[-1]
    sample_pattern = re.compile(r'(GCA|GCF|DRR|ERR)_?\d{6,9}(\.\d+)?')
    match = sample_pattern.search(filename)
    if match:
        return match.group(0)
    else:
        return filename.split("_")[0]


def report_main(annotation_folder: Union[str, Path], blast_file: Union[str, Path], cell_type: str, library: Union[str, Path], no_split: bool, metadata_folder: Union[str, Path]):
    df = pd.read_csv(blast_file, low_memory=False)

    df['% Mismatches of total alignment'] = (df['mismatches'] / df['alignment length']) * 100

    df[['SNPs', 'Insertions', 'Deletions']] = df['btop'].apply(lambda x: pd.Series(parse_btop(x)))
    split_query_df = df['query'].str.split('#', expand=True)
    df[['query', 'start', 'stop', 'strand', 'path', 'haplotype', 'tool', 'mapping_accuracy']] = split_query_df[[0, 1, 2, 3, 4, 5, 6, 7]]

    df["Sample"] = df["path"].apply(extract_sample)

    df = select_best_regions(df)

    known = df[df['% identity'] == 100.0]
    novel = df.groupby(['start', 'stop', 'path']).filter(filter_group)

    known.to_excel(f"{annotation_folder}/annotation_report_known.xlsx", index=False)
    novel.to_excel(f"{annotation_folder}/annotation_report_novel.xlsx", index=False)

