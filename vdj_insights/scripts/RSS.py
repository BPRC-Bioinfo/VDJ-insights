"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

import os
import subprocess

import pandas as pd
from Bio import SeqIO, motifs
from pathlib import Path
from Bio.Seq import Seq
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

from .util import load_config, calculate_available_resources, log_error


def open_files(cwd: Path) -> pd.DataFrame:
    """
    Opens and concatenates known and novel annotation reports.

    Args:
        cwd (Path): Current working directory.

    Returns:
        tuple:
            - pd.DataFrame: DataFrame with non-VDJ segments.
            - pd.DataFrameGroupBy: Grouped DataFrame by 'Region' and 'Segment' for VDJ segments.
    """
    data_known = pd.read_excel(cwd / "annotation" / "annotation_report_known.xlsx")
    data_novel = pd.read_excel(cwd / "annotation" / "annotation_report_novel.xlsx")
    data = pd.concat([data_known, data_novel])

    data_c = data[~data["Segment"].isin(["V", "D", "J"])]
    data_vdj = data[data["Segment"].isin(["V", "D", "J"])]

    vdj_grouped = data_vdj.groupby(["Region", "Segment"])
    return data_c, vdj_grouped


#@log_error()
def run_meme(locus_fasta_file_name: Path, meme_output: Path, rss_length: int, sum_lenght_seq: int) -> None:
    """
    Runs the MEME command to find motifs in the given FASTA file.

    Args:
        locus_fasta_file_name (Path): Path to the input FASTA file.
        meme_output (Path): Path to the output directory for MEME results.
        rss_length (int): Length of the RSS motif.
        sum_lenght_seq (int): Sum of the lengths of sequences in the FASTA file.
    """
    meme_command = f"meme {locus_fasta_file_name} -o {meme_output} -dna -mod zoops -nmotifs 1 -minw {rss_length} -maxw {rss_length} -maxsize {sum_lenght_seq}"
    if not meme_output.exists():
        subprocess.run(meme_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


@log_error()
def run_fimo(fimo_output: Path, meme_output: Path, locus_fasta_file_name: Path) -> None:
    """
    Runs the FIMO command to find motif occurrences in the given FASTA file.

    Args:
        fimo_output (Path): Path to the output directory for FIMO results.
        meme_output (Path): Path to the MEME output directory.
        locus_fasta_file_name (Path): Path to the input FASTA file.
    """
    fimo_command = f"fimo --thresh 0.0001 --o {fimo_output} {meme_output}/meme.txt {locus_fasta_file_name}"
    if not fimo_output.exists():
        subprocess.run(fimo_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def get_fimo_output(fimo_intput: Path) -> pd.DataFrame:
    """
    Reads the FIMO output file and processes it into a DataFrame.

    Args:
        fimo_intput (Path): Path to the FIMO output directory.

    Returns:
        pd.DataFrame: Processed DataFrame with FIMO output, including sequence index, scores, and p/q values.
    """
    fimo_output = fimo_intput / "fimo.txt"
    df_fimo = pd.read_csv(fimo_output, sep='\t')
    df_fimo["index_group_df"] = df_fimo["sequence name"].str.split("_").str[-1]
    df_fimo['index_group_df'] = df_fimo['index_group_df'].astype(int)
    return df_fimo


def process_group_locus(group_locus, df_fimo, rss_layout, rss_length):
    """
    Processes each group locus and updates it with FIMO results.

    Args:
        group_locus (pd.DataFrame): Grouped locus DataFrame.
        df_fimo (pd.DataFrame): FIMO output DataFrame.
        rss_layout (str): RSS layout type (e.g., "end_plus" or "start_minus").
        rss_length (int): Length of the RSS motif.

    Returns:
        pd.DataFrame: Updated group locus DataFrame with added FIMO results.
    """
    for index_segment, row in group_locus.iterrows():
        is_present = index_segment in df_fimo['index_group_df'].values and df_fimo.loc[df_fimo['index_group_df'] == index_segment, 'score'].values[0] > 0

        if rss_layout == "end_plus":
            direction = 3
        elif rss_layout == "start_minus":
            direction = 5

        group_locus.loc[index_segment, f"{direction}'-RSS"] = str(is_present)
        if is_present:
            group_locus.loc[index_segment, f"{direction}'-p-value"] = df_fimo.loc[df_fimo['index_group_df'] == index_segment, 'p-value'].values[0]
            group_locus.loc[index_segment, f"{direction}'-q-value"] = df_fimo.loc[df_fimo['index_group_df'] == index_segment, 'q-value'].values[0]
            group_locus.loc[index_segment, f"{direction}'-score"] = df_fimo.loc[df_fimo['index_group_df'] == index_segment, 'score'].values[0]
            group_locus.loc[index_segment, f"{direction}'-RSS seq"] = df_fimo.loc[df_fimo['index_group_df'] == index_segment, 'matched sequence'].values[0].upper()

    return group_locus


def process_variant(locus_gene_type, group_locus, config, output_base, cwd):
    """
    Processes each variant of a locus gene type by generating FASTA files and running MEME/FIMO.

    Args:
        locus_gene_type (tuple): Tuple representing the locus gene type (e.g., ("IGH", "V")).
        group_locus (pd.DataFrame): Grouped locus DataFrame.
        config (dict): Configuration dictionary.
        output_base (Path): Base output directory for results.
        cwd (Path): Current working directory.

    Returns:
        pd.DataFrame: Combined DataFrame with updated FIMO results for the group locus.
    """
    combined_results = pd.DataFrame()

    rss_variants = config['RSS_LAYOUT'].get("".join(locus_gene_type), {})
    for rss_variant in rss_variants.keys():
        rss_length = config['RSS_LENGTH'][rss_variant]
        mer1, mer2 = config['RSS_MERS'].get(str(rss_variant), [])

        locus_type = "".join(locus_gene_type)
        locus_fasta_file_name = cwd / output_base / "fasta"
        locus_fasta_file_name.mkdir(exist_ok=True, parents=True)
        locus_fasta_file_name = locus_fasta_file_name / f"{locus_type}_{rss_variant}.fasta"

        locus_with_gaps_fasta_file_name = cwd / output_base / "fasta_gaps"
        locus_with_gaps_fasta_file_name.mkdir(exist_ok=True, parents=True)
        locus_with_gaps_fasta_file_name = locus_with_gaps_fasta_file_name / f"{locus_type}_{rss_variant}.fasta"

        sum_lenght_seq = 0
        sum_seq = 0
        with open(locus_fasta_file_name, 'w') as locus_fasta_file, open(locus_with_gaps_fasta_file_name, 'w') as locus_fasta_file2:
            for index_segment, row in group_locus.iterrows():
                start_coord_segment = row["Start coord"]
                end_coord_segment = row["End coord"]
                strand = row["Strand"]
                fasta_path = row["Path"]

                with open(fasta_path, 'r') as fasta_file:
                    record = SeqIO.read(fasta_file, 'fasta')

                rss_layout = rss_variants[rss_variant][strand]
                if rss_layout == "end_plus":
                    start_rss, end_rss = [end_coord_segment, (end_coord_segment + rss_length)]
                elif rss_layout == "start_minus":
                    start_rss, end_rss = [(start_coord_segment - rss_length), start_coord_segment]

                rss = record.seq[start_rss:end_rss].replace("-", "N")
                if rss:
                    if strand == '-':
                        rss = str(Seq(str(rss)).reverse_complement())

                    seq_l, spacer, seq_r = rss[0:mer1], rss[mer1:-mer2], rss[-mer2:]
                    locus_fasta_file.write(f">{row['Old name-like']}_{index_segment}\n{seq_l}{spacer}{seq_r}\n")

                    if row["Function"] != "P" and row["Status"] == "Known":
                        spacer = len(spacer) * "N"
                        locus_fasta_file2.write(f">{row['Old name-like']}_{index_segment}\n{seq_l}{spacer}{seq_r}\n")
                        sum_seq = sum_seq + 1

                    sum_lenght_seq += len(rss)

        if sum_seq > 1:
            meme_output = cwd / output_base / "meme_output" / f"{locus_type}_{rss_variant}"
            os.makedirs(cwd / output_base / "meme_output", exist_ok=True)
            run_meme(locus_fasta_file_name=locus_with_gaps_fasta_file_name, meme_output=meme_output, rss_length=rss_length, sum_lenght_seq=sum_lenght_seq)

            fimo_output = cwd / output_base / "fimo_output" / f"{locus_type}_{rss_variant}"
            os.makedirs(cwd / output_base / "fimo_output", exist_ok=True)
            run_fimo(fimo_output=fimo_output, meme_output=meme_output, locus_fasta_file_name=locus_fasta_file_name)

            df_fimo = get_fimo_output(fimo_intput=fimo_output)
            group_locus = process_group_locus(group_locus=group_locus, df_fimo=df_fimo, rss_layout=rss_layout, rss_length=rss_length)
    combined_results = pd.concat([combined_results, group_locus])
    return combined_results


def add_suffix_to_short_name(group):
    """
    Adds a suffix to duplicate 'Short name' entries to ensure uniqueness.

    Args:
        group (pd.DataFrame): Grouped DataFrame containing segment data.

    Returns:
        pd.DataFrame: Updated DataFrame with unique 'Short name' values.
    """
    unique_sequences = group['Old name-like seq'].unique()
    if len(unique_sequences) > 1:
        seq_to_suffix = {seq: f"_{i + 1}" for i, seq in enumerate(unique_sequences)}
        group['Short name'] = group['Short name'] + group['Old name-like seq'].map(seq_to_suffix)
    return group


def main_rss(threads: int = 8) -> None:
    """
    Main function to process RSS annotations in parallel using multiple threads.

    Args:
        threads (int, optional): Number of threads to use for parallel processing. Defaults to 8.
    """
    cwd = Path.cwd()
    config = load_config(cwd / "config" / "config.yaml")

    output_base = cwd / "RSS"
    data_c, vdj_grouped = open_files(cwd)

    combined_df = data_c

    max_jobs = calculate_available_resources(max_cores=threads, threads=1, memory_per_process=2)
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = [
            executor.submit(process_variant, locus_gene_type, group_locus, config, output_base, cwd)
            for locus_gene_type, group_locus in vdj_grouped
        ]
        with tqdm(total=vdj_grouped.ngroups, desc="Processing RSS", unit='task') as pbar:
            for future in as_completed(futures):
                future.result()
                combined_df = pd.concat([combined_df, future.result()])
                pbar.update(1)

    known = combined_df[combined_df["Status"] == "Known"]
    novel = combined_df[combined_df["Status"] == "Novel"]

    novel = novel.groupby('Short name', group_keys=False).apply(add_suffix_to_short_name)

    known = known.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
    novel = novel.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])

    known.to_excel(cwd / "annotation" / "annotation_report_known_rss.xlsx", index=False)
    novel.to_excel(cwd / "annotation" / "annotation_report_novel_rss.xlsx", index=False)


if __name__ == '__main__':
    main_rss(threads=12)
