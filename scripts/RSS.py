import multiprocessing
import re
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq

from concurrent.futures import ProcessPoolExecutor, as_completed
from util import make_dir, load_config, seperate_annotation, calculate_available_resources
from logger import console_logger, file_logger
from property import log_error

from overlap import remove_overlapping_segments
from tqdm import tqdm

pd.set_option('display.max_rows', None)

console_log = console_logger(__name__)
file_log = file_logger(__name__)


@log_error()
def write_fasta_file(dictionary, folder: str):
    """
    Writes multiple FASTA files for each region and segment in the provided dictionary.
    Each file is named "{key}{region}.fasta" and stored in the specified folder.
    Headers in the FASTA files are made unique by appending an incrementing number.
    """
    make_dir(folder)
    for key, value in dictionary.items():
        for region, segments in value.items():
            ffile = folder / f"{key}{region}.fasta"
            with open(ffile, 'w') as f:
                for header, segment in segments.items():
                    for x, rss in enumerate(segment):
                        f.write(f">{header}-{x+1}\n{rss}\n")


def calculate_position(position, variant_length, operation):
    """
    Calculates the adjusted start or end position for cutting the RSS sequence based on the given operation.
    It either returns an end coordinate ("end_plus") or a start coordinate ("start_minus").
    """
    if operation.endswith("plus"):
        return str(position + variant_length)
    elif operation.endswith("minus"):
        return str(position - variant_length)
    else:
        file_log.error(f"Invalid operation provided: {operation}")
        raise ValueError(f"Invalid operation: {operation}")


def rss_type(start, end, combi, rss_variant, strand, config):
    """
    Determines the RSS coordinates based on the combination of segment and region, and whether the sequence is on the forward or reverse strand.
    It uses the start or end coordinate to calculate the appropriate RSS coordinates.
    """
    lengths = config['RSS_LENGTH']
    variants = config['RSS_LAYOUT'].get(combi, {})
    variant_length = lengths.get(str(rss_variant), 0)
    layout = variants.get(str(rss_variant), {}).get(strand, "")

    if layout == "end_plus":
        return [str(end), calculate_position(end, variant_length, layout)]
    elif layout == "start_minus":
        return [calculate_position(start, variant_length, layout), str(start)]
    else:
        file_log.error(f"Invalid layout or combination: {combi}, {strand}")
        raise KeyError(f"Invalid combination or strand: {combi}, {strand}")


@log_error()
def fetch_sequence(row, segment, rss_variant, config):
    """
    Extracts the RSS sequence from the specified row of the DataFrame based on the segment and RSS variant.
    The sequence is cut from the appropriate FASTA file based on the coordinates.
    If the sequence is on the reverse strand, it is reverse complemented.
    """
    fasta, start, end, strand, region = row[['Path', 'Start coord', 'End coord', 'Strand', 'Region']]
    with open(fasta, 'r') as fasta_file:
        record = SeqIO.read(fasta_file, 'fasta')
        rss_start, rss_stop = rss_type(start, end, f"{region}{segment}", rss_variant, strand, config)
        rss = record.seq[int(rss_start):int(rss_stop)]
        if strand == '-':
            rss = str(Seq(str(rss)).reverse_complement())
    return rss


@log_error()
def add_to_dict(query, dictionary, rss):
    """
    Adds an RSS sequence to the specified dictionary under the given query.
    If the query does not exist, it creates a new entry.
    """
    dictionary.setdefault(query, []).append(rss)


@log_error()
def get_mers(rss, rss_variant, config):
    """
    Extracts the heptamer and nonamer sequences from the RSS based on the RSS variant.
    """
    mers = config["RSS_MERS"].get(str(rss_variant), [])
    mer1, mer2 = mers
    return [rss[0:mer1], rss[-mer2:]]


def add_to_row(row, mer1, mer2, rss_variant):
    """
    Adds the heptamer and nonamer sequences to the specified row of the DataFrame.
    Depending on the length of mer1, it determines which sequence represents the heptamer or nonamer.
    """
    heptamer, nonamer = (''.join(mer1), ''.join(mer2)) if len(
        mer1) == 7 else (''.join(mer2), ''.join(mer1))
    row[f"{rss_variant}_heptamer"], row[f"{rss_variant}_nonamer"] = heptamer, nonamer
    return row


@log_error()
def add_segment(query, segment, region, separated_segments, rss_sequence, function):
    """
    Adds the RSS sequence to the separated_segments dictionary based on the segment, region, and function.
    Only adds the RSS sequence if the function is not "P" (pseudogene).
    """
    if function != "P":
        add_to_dict(query, separated_segments.setdefault(region, {}).setdefault(segment, {}), str(rss_sequence))


@log_error()
def add_base_rss_parts(row, config):
    """
    Adds base RSS components (heptamer and nonamer) to the row based on the segment type and RSS variant.
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    row["12_heptamer"], row["12_nonamer"], row["23_heptamer"], row["23_nonamer"], row["24_heptamer"], row["24_nonamer"], row["13_heptamer"], row["13_nonamer"] = "", "", "", "", "", "", "", ""
    rss_variants = config['RSS_LAYOUT'].get(full, {}).keys()

    for rss_variant in rss_variants:
        rss_sequence = fetch_sequence(row, segment, rss_variant, config)
        mer1, mer2 = get_mers(rss_sequence, rss_variant, config)
        add_to_row(row, mer1, mer2, rss_variant)

    return row


@log_error()
def run_meme(out, rss_file, rss_variant, THREADS=8):
    """
    Executes the MEME suite to identify motifs in the RSS sequences.
    Chooses between a single-command or multi-command MEME run based on the number of sequences in the input file.
    """
    total_nucleotides = sum(len(record.seq) for record in SeqIO.parse(rss_file, "fasta"))

    multi_command = f"meme {rss_file} -o {out} -dna -mod zoops -nmotifs 1 -minw {rss_variant} -maxsize {total_nucleotides}"
    single_command = f"meme {rss_file} -o {out} -dna -mod anr -nmotifs 1 -minw {rss_variant}"
    amount = subprocess.run(f"cat {rss_file} | egrep '^>' | wc -l", shell=True, capture_output=True)

    if int(amount.stdout) > 1:
        subprocess.run(multi_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        subprocess.run(single_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


@log_error()
def get_reference_mers(regex_string, rss_variant, config):
    """
    Extracts the heptamer and nonamer sequences from a MEME-generated regular expression string.
    """
    rss = re.findall(r'\[[^\]]*\]|.', regex_string)
    mers = config["RSS_MERS"].get(str(rss_variant), [])
    mer1, mer2 = mers
    return rss[0:mer1], rss[-mer2:]


@log_error()
def make_ref_dict(segment, ref_rss_dict, mer1, mer2):
    """
    Updates the reference RSS dictionary with the heptamer and nonamer sequences for a specific segment.
    """
    if len(mer1) == 7:
        ref_rss_dict.setdefault(segment, {}).setdefault(
            "heptamer", ''.join(mer1))
        ref_rss_dict.setdefault(segment, {}).setdefault(
            "nonamer", ''.join(mer2))
    else:
        ref_rss_dict.setdefault(segment, {}).setdefault(
            "heptamer", ''.join(mer2))
        ref_rss_dict.setdefault(segment, {}).setdefault(
            "nonamer", ''.join(mer1))
    return ref_rss_dict


def run_meme_on_file(meme_directory: Path, rss_file: Path, config: dict) -> None:
    """
    Helper function to run the MEME suite on a single RSS file.
    """
    stem = rss_file.stem
    rss_variant = stem.split("_")[-1]
    RSS_convert = config.get("RSS_LENGTH", {})
    out = meme_directory / stem
    meme = out / "meme.txt"

    if not meme.exists():
        run_meme(out, rss_file, RSS_convert.get(rss_variant))


@log_error()
def create_meme_directory(meme_tasks):
    """
    Processes all rss_files, running MEME to generate motifs, with a single progress bar.
    """
    total_tasks = len(meme_tasks)
    max_jobs = calculate_available_resources(max_cores=24, threads=8, memory_per_process=2)

    with tqdm(total=total_tasks, desc="Processing Files") as pbar:
        with ProcessPoolExecutor(max_workers=max_jobs) as executor:
            futures = {
                executor.submit(run_meme_on_file, task['meme_directory'], task['rss_file'], task['config']): task
                for task in meme_tasks
            }
            for future in as_completed(futures):
                task = futures[future]
                try:
                    future.result()
                except Exception as e:
                    file_log.error(f"Error processing {task['rss_file']}: {e}")
                pbar.update(1)


@log_error()
def make_reference_rss(ref_meme_directory, config):
    """
    Generates a reference dictionary of RSS heptamers and nonamers from MEME results.
    """
    ref_rss_dict = {}

    for meme in ref_meme_directory.iterdir():
        meme_text = meme / "meme.txt"
        command = f'cat {meme_text} | egrep -A2 "regular expression"'
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        hits = result.stdout.replace("-", "").replace("\t", "").strip().split("\n")
        hits = [hit for hit in hits if hit]
        if hits:
            split_stem = meme.stem.split("_")
            rss_variant = int(split_stem[1])
            val1, val2 = get_reference_mers(hits[1], rss_variant, config)
            make_ref_dict(meme.stem, ref_rss_dict, val1, val2)
    return ref_rss_dict


def seperate_pattern(pattern):
    """
    Separates a regular expression pattern into individual components, preserving elements in brackets.
    """
    seperated = list()
    for part in pattern.split('['):
        if ']' in part:
            reg, rest = part.split(']', 1)
            seperated.append(f'[{reg}]')
            seperated.extend(rest)
        else:
            seperated.extend(part)
    return seperated


def rss_mismatch(rss, pattern):
    """
    Compares an RSS sequence to a regular expression pattern and counts the number of mismatches.
    """
    pattern, rss = seperate_pattern(pattern), [*rss]
    return [bool(re.match(x, i)) for x, i in zip(pattern, rss)].count(False)


def check_ref_rss(row, ref_rss_dict, rss_variant):
    """
    Compares the heptamer and nonamer sequences of a segment with the reference sequences.
    Updates the DataFrame row with match status and reference sequences.
    """
    region_segment_key = f"{row['Region']}{row['Segment']}_{rss_variant}"
    ref_rss = ref_rss_dict.get(region_segment_key, "")

    for motif in ["heptamer", "nonamer"]:
        ref_seq = ref_rss.get(motif, "") if ref_rss else ""
        query_seq = row.get(f"{rss_variant}_{motif}", "")

        if ref_seq:
            matches = rss_mismatch(query_seq, ref_seq)
            row[f"{rss_variant}_ref_{motif}"] = ref_seq
            row[f"{rss_variant}_{motif}_matched"] = matches <= 1
        else:
            row[f"{rss_variant}_{motif}_matched"] = True

    return row


@log_error()
def combine_df(original_df, new_df):
    """
    Merges the original DataFrame with the new DataFrame containing RSS heptamers and nonamers.
    Ensures that all relevant columns are included in the merged DataFrame.
    """
    columns_to_merge = [
        'Reference', 'Old name-like', 'Start coord', 'End coord',
        '12_heptamer', '12_ref_heptamer', '12_heptamer_matched',
        '12_nonamer', '12_ref_nonamer', '12_nonamer_matched',
        '13_heptamer', '13_ref_heptamer', '13_heptamer_matched',
        '13_nonamer', '13_ref_nonamer', '13_nonamer_matched',
        '23_heptamer', '23_ref_heptamer', '23_heptamer_matched',
        '23_nonamer', '23_ref_nonamer', '23_nonamer_matched',
        '24_heptamer', '24_ref_heptamer', '24_heptamer_matched',
        '24_nonamer', '24_ref_nonamer', '24_nonamer_matched',
    ]

    for column in columns_to_merge:
        if column not in original_df.columns:
            original_df[column] = ''
    combined_df = pd.merge(new_df, original_df[columns_to_merge], on=['Reference', 'Old name-like', 'Start coord', 'End coord'], how='left')

    return combined_df


@log_error()
def apply_check_ref_rss(row, ref_rss_dict, config, options):
    """
    Applies the reference RSS check to each row in the DataFrame.
    Validates the segment against reference RSS sequences and updates the row with match status.
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    rss_variants = config['RSS_LAYOUT'].get(full, {}).keys()
    for rss_variant in list(rss_variants):
        if full in options:
            check_ref_rss(row, ref_rss_dict, rss_variant)

    return row


@log_error()
def create_dict(row, separated_segments, config):
    """
    Creates a dictionary entry for the RSS sequence based on the current row of the DataFrame.
    Adds the segment, region, and function to the separated_segments dictionary.
    """
    query, region, segment, function = row[["Old name-like", "Region", "Segment", "Function"]]
    combi = region + segment
    rss_variants = config['RSS_LAYOUT'].get(combi, {}).keys()

    for rss_variant in rss_variants:
        rss_sequence = fetch_sequence(row, segment, rss_variant, config)
        add_segment(query, f"{segment}_{rss_variant}", region, separated_segments, rss_sequence, function)

    return row


@log_error()
def create_RSS_files(df, RSS_directory, config):
    """
    Generates RSS files from the DataFrame and saves them to the specified directory.
    """
    if not RSS_directory.exists():
        separated_segments = {}
        df.apply(lambda row: create_dict(row, separated_segments, config), axis=1)
        write_fasta_file(separated_segments, RSS_directory)


@log_error()
def create_all_RSS_meme_files(cwd, df, config):
    """
    Generates all necessary RSS and MEME files for the analysis.
    """
    base = cwd / "RSS"
    df_100 = pd.read_excel(cwd / 'annotation' / 'annotation_report_known.xlsx')
    combined = pd.concat([df_100, df], axis=0)
    datasets = [df_100, df, combined]
    rss_filenames = ["reference_RSS", "new_RSS", "combined_RSS"]
    meme_directories = ["reference_meme", "new_meme", "complete_meme"]
    
    # First, create all the RSS files
    for dataset, rss_filename in zip(datasets, rss_filenames):
        create_RSS_files(dataset, base / rss_filename, config)
    
    # Now, collect all the RSS files and their corresponding output meme directories
    meme_tasks = []
    for rss_filename, meme_dir_name in zip(rss_filenames, meme_directories):
        RSS_directory = base / rss_filename
        rss_files = [rss_file for rss_file in RSS_directory.iterdir() if rss_file.suffix in ['.fasta', '.fa']]
        meme_directory = base / meme_dir_name
        make_dir(meme_directory)
        for rss_file in rss_files:
            meme_tasks.append({
                'rss_file': rss_file,
                'meme_directory': meme_directory,
                'config': config
            })
    
    # Process all MEME tasks with a progress bar
    create_meme_directory(meme_tasks)
    
    return base / meme_directories[0]


@log_error()
def update_df(df, ref_rss_dict, config, options):
    """
    Updates the DataFrame with base RSS components and checks against reference RSS sequences.
    Adds columns for heptamer and nonamer sequences and their match status.
    """
    df = df.apply(lambda row: add_base_rss_parts(row, config), axis=1)
    df = df.apply(lambda row: apply_check_ref_rss(row, ref_rss_dict, config, options), axis=1)
    return df


@log_error()
def create_rss_excel_file(cwd, final_df, no_split):
    """
    Process the DataFrame to remove non-best overlapping rows and export
    the remaining data into separate Excel files for 'Novel' and 'Known' statuses.
    """
    novel_df = final_df.query("Status == 'Novel'")
    known_df = final_df.query("Status == 'Known'")
    for df, filename in zip([novel_df, known_df], ["annotation_report_novel", "annotation_report_known"]):
        if not no_split:
            file_log.info("Creating individual sample excel files...")
            df.groupby("Sample").apply(lambda group: seperate_annotation(
                group, cwd / "annotation", f"{filename}_rss.xlsx"))
        wrtie_rss_excel_file(cwd, df, filename)


@log_error()
def wrtie_rss_excel_file(cwd, df, filename):
    """
    Creates an extended Excel file with the updated DataFrame.
    Merges the original annotation report with the new RSS data and saves it as an Excel file.
    """
    file_log.info(f"Generating {filename}_rss.xlsx!")
    df.to_excel(cwd / 'annotation' / f'{filename}_rss.xlsx', index=False)


@log_error()
def RSS_main(no_split=False, threads=8):
    """
    Main function of the RSS creation script.
    """
    complete_df = pd.DataFrame()
    cwd = Path.cwd()

    config = load_config(cwd / "config" / "config.yaml")
    options = set(config.get("RSS_LAYOUT", {}).keys())

    df1 = pd.read_excel(cwd / 'annotation' / 'annotation_report_novel.xlsx')
    df2 = pd.read_excel(cwd / 'annotation' / 'annotation_report_known.xlsx')
    ref_meme_directory = create_all_RSS_meme_files(cwd, df1, config)
    ref_rss_dict = make_reference_rss(ref_meme_directory, config)

    for df, filename in zip([df1, df2], ["annotation_report_novel", "annotation_report_known"]):
        df = update_df(df, ref_rss_dict, config, options)
        reference_df = combine_df(df, pd.read_excel(cwd / 'annotation' / f'{filename}.xlsx')).drop_duplicates()
        complete_df = pd.concat([complete_df, reference_df], ignore_index=True)
    final_df = remove_overlapping_segments(complete_df)
    create_rss_excel_file(cwd, final_df, no_split)


if __name__ == '__main__':
    RSS_main()
