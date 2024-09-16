import re
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq

from util import make_dir, load_config
from logger import custom_logger
from overlap import remove_overlapping_segments
pd.set_option('display.max_rows', None)

"""
Used CLI packages:
    1. yaml
    2. pandas
    3. openpyxl
    4. biopython
    5. meme
    6. imagemagick

Used Python packages:
    1. yaml
    2. pandas
    3. openpyxl
    4. biopython
"""
# Method for logging current states of the program.
logger = custom_logger(__name__)


def write_fasta_file(dictionary, folder):
    """
    Writes multiple FASTA files for each region and segment in the provided dictionary.
    Each file is named "{key}{region}.fasta" and stored in the specified folder.
    Headers in the FASTA files are made unique by appending an incrementing number.

    Args:
        dictionary (dict): Dictionary containing regions as keys and segments with their corresponding RSS sequences as values.
        folder (Path): Path to the directory where FASTA files will be saved.

    Raises:
        OSError: If there is an issue writing to the files, logs the error and raises an exception.
    """
    make_dir(folder)
    try:
        for key, value in dictionary.items():
            for region, segments in value.items():
                ffile = folder / f"{key}{region}.fasta"
                with open(ffile, 'w') as f:
                    for header, segment in segments.items():
                        for x, rss in enumerate(segment):
                            f.write(f">{header}-{x+1}\n{rss}\n")
    except OSError as e:
        logger.error(f"Failed to write to FASTA file {ffile}: {e}")
        raise


def calculate_position(position, variant_length, operation):
    """
    Calculates the adjusted start or end position for cutting the RSS sequence based on the given operation.
    It either returns an end coordinate ("end_plus") or a start coordinate ("start_minus").

    Args:
        position (int): The coordinate that either contains the start or end of an RSS sequence.
        variant_length (int): The length of the RSS sequence to be extracted.
        operation (str): Operation to perform, which can be "start_minus" or "end_plus".

    Returns:
        str: The adjusted coordinate as a string.

    Raises:
        ValueError: If an invalid operation is provided, logs the error and raises an exception.
    """
    if operation.endswith("plus"):
        return str(position + variant_length)
    elif operation.endswith("minus"):
        return str(position - variant_length)
    else:
        logger.error(f"Invalid operation provided: {operation}")
        raise ValueError(f"Invalid operation: {operation}")


def rss_type(start, end, combi, rss_variant, strand, config):
    """
    Determines the RSS coordinates based on the combination of segment and region, and whether the sequence is on the forward or reverse strand.
    It uses the start or end coordinate to calculate the appropriate RSS coordinates.

    Args:
        start (str): The start coordinate of the VDJ segment sequence.
        end (str): The end coordinate of the VDJ segment sequence.
        combi (str): Combination of the V, D, or J segment and the region.
        rss_variant (int): Type of RSS variant (e.g., 12 or 23).
        strand (str): Indicator of the sequence orientation ("+" for forward or "-" for reverse).

    Returns:
        list: A list containing two integers representing the start and end coordinates of the RSS variant.

    Raises:
        KeyError: If the combination or strand is not found in the configuration, logs the error and raises an exception.
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
        logger.error(f"Invalid layout or combination: {combi}, {strand}")
        raise KeyError(f"Invalid combination or strand: {combi}, {strand}")


def fetch_sequence(row, segment, rss_variant, config):
    """
    Extracts the RSS sequence from the specified row of the DataFrame based on the segment and RSS variant.
    The sequence is cut from the appropriate FASTA file based on the coordinates.
    If the sequence is on the reverse strand, it is reverse complemented.

    Args:
        row (pd.Series): Current row from the DataFrame containing relevant segment information.
        segment (str): Type of segment (e.g., V, D, or J).
        rss_variant (int): Type of RSS variant (e.g., 12 or 23).

    Returns:
        str: Extracted RSS sequence.

    Raises:
        FileNotFoundError: If the FASTA file is not found, logs the error and raises an exception.
        KeyError: If the segment or region combination is invalid, logs the error and raises an exception.
    """
    fasta, start, end, strand, region = row[[
        'Path', 'Start coord', 'End coord', 'Strand', 'Region']]
    try:
        with open(fasta, 'r') as fasta_file:
            record = SeqIO.read(fasta_file, 'fasta')
            rss_start, rss_stop = rss_type(
                start, end, f"{region}{segment}", rss_variant, strand, config)
            rss = record.seq[int(rss_start):int(rss_stop)]
            if strand == '-':
                rss = str(Seq(str(rss)).reverse_complement())
        return rss
    except FileNotFoundError as e:
        logger.error(f"FASTA file not found: {fasta}")
        raise
    except KeyError as e:
        logger.error(f"Invalid segment or region combination: {e}")
        raise


def add_to_dict(query, dictionary, rss):
    """
    Adds an RSS sequence to the specified dictionary under the given query.
    If the query does not exist, it creates a new entry.

    Args:
        query (str): Name of the segment.
        dictionary (dict): Dictionary to store the RSS sequences.
        rss (str): RSS sequence to add.

    Raises:
        TypeError: If the provided dictionary is not valid, logs the error and raises an exception.
    """
    try:
        dictionary.setdefault(query, []).append(rss)
    except TypeError as e:
        logger.error(f"Failed to add to dictionary: {e}")
        raise


def get_mers(rss, rss_variant, config):
    """
    Extracts the heptamer and nonamer sequences from the RSS based on the RSS variant.

    Args:
        rss (str): RSS sequence from which to extract the heptamer and nonamer.
        rss_variant (int): Type of RSS variant (e.g., 12 or 23).

    Returns:
        list: List containing the heptamer and nonamer sequences.

    Raises:
        KeyError: If the RSS variant is not found in the configuration, logs the error and raises an exception.
    """
    try:
        mers = config["RSS_MERS"].get(str(rss_variant), [])
        mer1, mer2 = mers
        return [rss[0:mer1], rss[-mer2:]]
    except KeyError as e:
        logger.error(f"RSS variant not found in configuration: {e}")
        raise


def add_to_row(row, mer1, mer2, rss_variant):
    """
    Adds the heptamer and nonamer sequences to the specified row of the DataFrame.
    Depending on the length of mer1, it determines which sequence represents the heptamer or nonamer.

    Args:
        row (pd.Series): Current row from the DataFrame.
        mer1 (list): List representing the heptamer or nonamer sequence.
        mer2 (list): List representing the heptamer or nonamer sequence.
        rss_variant (int): Type of RSS variant (e.g., 12 or 23).

    Returns:
        pd.Series: Updated row with the added heptamer and nonamer sequences.
    """
    heptamer, nonamer = (''.join(mer1), ''.join(mer2)) if len(
        mer1) == 7 else (''.join(mer2), ''.join(mer1))
    row[f"{rss_variant}_heptamer"], row[f"{
        rss_variant}_nonamer"] = heptamer, nonamer
    return row


def add_segment(query, segment, region, separated_segments, rss_sequence, function):
    """
    Adds the RSS sequence to the separated_segments dictionary based on the segment, region, and function.
    Only adds the RSS sequence if the function is not "P" (pseudogene).

    Args:
        query (str): Name of the segment.
        segment (str): Type of segment (e.g., V_12, D_23).
        region (str): Region associated with the segment.
        separated_segments (dict): Dictionary to store the RSS sequences.
        rss_sequence (str): RSS sequence to add.
        function (str): Functional status of the segment ("P" for pseudogene, "F/ORF" for functional).

    Raises:
        TypeError: If the provided dictionary is not valid, logs the error and raises an exception.
    """
    if function != "P":
        try:
            add_to_dict(query, separated_segments.setdefault(
                region, {}).setdefault(segment, {}), str(rss_sequence))
        except TypeError as e:
            logger.error(f"Failed to add segment to dictionary: {e}")
            raise


def add_base_rss_parts(row, config):
    """
    Adds base RSS components (heptamer and nonamer) to the row based on the segment type and RSS variant.
    Initializes the 12_heptamer, 12_nonamer, 23_heptamer, and 23_nonamer columns with empty strings.
    Populates these columns with the appropriate sequences based on the segment type and RSS variant.

    Args:
        row (pd.Series): Current row from the DataFrame.

    Returns:
        pd.Series: Updated row with the added RSS components.

    Raises:
        KeyError: If the segment or region combination is not found in the configuration, logs the error and raises an exception.
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    row["12_heptamer"], row["12_nonamer"], row["23_heptamer"], row["23_nonamer"] = "", "", "", ""
    rss_variants = config['RSS_LAYOUT'].get(full, {}).keys()

    try:
        for rss_variant in rss_variants:
            rss_sequence = fetch_sequence(row, segment, rss_variant, config)
            mer1, mer2 = get_mers(rss_sequence, rss_variant, config)
            add_to_row(row, mer1, mer2, rss_variant)
    except KeyError as e:
        logger.error(f"Invalid segment or region combination: {e}")
        raise

    return row


def run_meme(out, rss_file, rss_variant):
    """
    Executes the MEME suite to identify motifs in the RSS sequences.
    Chooses between a single-command or multi-command MEME run based on the number of sequences in the input file.

    Args:
        out (Path): Path to the output directory for MEME results.
        rss_file (Path): Path to the RSS input file.
        rss_variant (int): Length of the RSS variant (e.g., 12 or 23).

    Raises:
        CalledProcessError: If the MEME command fails, logs the error and raises an exception.
    """
    total_nucleotides = sum(len(record.seq)
                            for record in SeqIO.parse(rss_file, "fasta"))
    logger.debug(total_nucleotides)
    multi_command = f"meme {
        rss_file} -o {out} -dna -mod zoops -nmotifs 1 -minw {rss_variant} -maxsize {total_nucleotides}"
    single_command = f"meme {
        rss_file} -o {out} -dna -mod anr -nmotifs 1 -minw {rss_variant}"
    amount = subprocess.run(
        f"cat {rss_file} | egrep '^>' | wc -l", shell=True, capture_output=True)

    try:
        if int(amount.stdout) > 1:
            subprocess.run(multi_command, shell=True, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            subprocess.run(single_command, shell=True, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logger.error(f"MEME command failed: {e}")
        raise


def get_reference_mers(regex_string, rss_variant, config):
    """
    Extracts the heptamer and nonamer sequences from a MEME-generated regular expression string.

    Args:
        regex_string (str): Regular expression string containing the RSS sequence.
        rss_variant (int): Length of the RSS variant (e.g., 12 or 23).

    Returns:
        list: List containing the extracted heptamer and nonamer sequences.

    Raises:
        ValueError: If the regular expression string is invalid, logs the error and raises an exception.
    """
    try:
        rss = re.findall(r'\[[^\]]*\]|.', regex_string)
        mers = config["RSS_MERS"].get(str(rss_variant), [])
        mer1, mer2 = mers
        return rss[0:mer1], rss[-mer2:]
    except ValueError as e:
        logger.error(f"Invalid regular expression string: {e}")
        raise


def make_ref_dict(segment, ref_rss_dict, mer1, mer2):
    """
    Updates the reference RSS dictionary with the heptamer and nonamer sequences for a specific segment.
    Depending on the length of mer1, it determines whether it is the heptamer or nonamer.

    Args:
        segment (str): Type of segment (e.g., V, D, or J).
        ref_rss_dict (dict): Dictionary to store reference heptamers and nonamers.
        mer1 (list): List representing the heptamer or nonamer sequence.
        mer2 (list): List representing the heptamer or nonamer sequence.

    Returns:
        dict: Updated reference RSS dictionary.

    Raises:
        TypeError: If the provided segment or sequence is invalid, logs the error and raises an exception.
    """
    try:
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
    except TypeError as e:
        logger.error(f"Failed to update reference dictionary: {e}")
        raise


def create_meme_directory(meme_directory, RSS_directory, config):
    """
    Creates a directory for storing MEME results and runs the MEME suite on RSS sequences.
    Iterates over all RSS files in the specified directory, running MEME to generate motifs.

    Args:
        meme_directory (Path): Path to the directory where MEME results will be saved.
        RSS_directory (Path): Path to the directory containing RSS FASTA files.

    Raises:
        OSError: If the directory creation or file writing fails, logs the error and raises an exception.
    """
    make_dir(meme_directory)

    try:
        for rss_file in Path(RSS_directory).iterdir():
            stem = rss_file.stem
            rss_variant = stem.split("_")[-1]
            RSS_convert = config.get("RSS_LENGTH", {})
            out = meme_directory / stem
            meme = out / "meme.txt"
            if not meme.exists():
                run_meme(out, rss_file, RSS_convert[rss_variant])
    except OSError as e:
        logger.error(f"Failed to create meme directory or run MEME: {e}")
        raise


def make_reference_rss(ref_meme_directory, config):
    """
    Generates a reference dictionary of RSS heptamers and nonamers from MEME results.
    Parses the MEME output files to extract heptamer and nonamer sequences and stores them in a dictionary.

    Args:
        ref_meme_directory (Path): Path to the directory containing reference MEME results.

    Returns:
        dict: Dictionary containing reference heptamers and nonamers for all regions and segments.

    Raises:
        subprocess.CalledProcessError: If the command to extract MEME results fails, logs the error and raises an exception.
    """
    ref_rss_dict = {}

    try:
        for meme in ref_meme_directory.iterdir():
            meme_text = meme / "meme.txt"
            command = f'cat {meme_text} | egrep -A2 "regular expression"'
            result = subprocess.run(command, shell=True,
                                    capture_output=True, text=True)
            hits = result.stdout.replace(
                "-", "").replace("\t", "").strip().split("\n")
            hits = [hit for hit in hits if hit]
            if hits:
                split_stem = meme.stem.split("_")
                rss_variant = int(split_stem[1])
                val1, val2 = get_reference_mers(hits[1], rss_variant, config)
                make_ref_dict(meme.stem, ref_rss_dict, val1, val2)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to extract MEME results: {e}")
        raise

    return ref_rss_dict


def seperate_pattern(pattern):
    """
    Separates a regular expression pattern into individual components, preserving elements in brackets.

    Args:
        pattern (str): Regular expression pattern to separate.

    Returns:
        list: List of separated pattern components.
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

    Args:
        rss (str): RSS sequence to compare.
        pattern (str): Regular expression pattern to match against.

    Returns:
        int: Number of mismatches between the RSS sequence and the pattern.
    """
    pattern, rss = seperate_pattern(pattern), [*rss]
    return [bool(re.match(x, i)) for x, i in zip(pattern, rss)].count(False)


def check_ref_rss(row, ref_rss_dict, rss_variant):
    """
    Compares the heptamer and nonamer sequences of a segment with the reference sequences.
    Updates the DataFrame row with match status and reference sequences.

    Args:
        row (pd.Series): Current row from the DataFrame.
        ref_rss_dict (dict): Dictionary containing reference RSS sequences.
        rss_variant (int): Length of the RSS variant (e.g., 12 or 23).

    Returns:
        pd.Series: Updated row with match status and reference sequences.

    Raises:
        KeyError: If the segment or region combination is not found in the reference dictionary, logs the error and raises an exception.
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


def combine_df(original_df, new_df):
    """
    Merges the original DataFrame with the new DataFrame containing RSS heptamers and nonamers.
    Ensures that all relevant columns are included in the merged DataFrame.

    Args:
        original_df (pd.DataFrame): Original DataFrame with base data.
        new_df (pd.DataFrame): New DataFrame with added heptamers and nonamers.

    Returns:
        pd.DataFrame: Combined DataFrame with merged data.

    Raises:
        ValueError: If the merge operation fails, logs the error and raises an exception.
    """
    columns_to_merge = [
        'Reference', 'Old name-like', 'Start coord', 'End coord',
        '12_heptamer', '12_ref_heptamer', '12_heptamer_matched',
        '12_nonamer', '12_ref_nonamer', '12_nonamer_matched',
        '23_heptamer', '23_ref_heptamer', '23_heptamer_matched',
        '23_nonamer', '23_ref_nonamer', '23_nonamer_matched',
    ]

    try:
        for column in columns_to_merge:
            if column not in original_df.columns:
                original_df[column] = ''
        combined_df = pd.merge(new_df, original_df[columns_to_merge],
                               on=['Reference', 'Old name-like',
                                   'Start coord', 'End coord'],
                               how='left')
    except ValueError as e:
        logger.error(f"Failed to merge DataFrames: {e}")
        raise

    return combined_df


def apply_check_ref_rss(row, ref_rss_dict, config, options):
    """
    Applies the reference RSS check to each row in the DataFrame.
    Validates the segment against reference RSS sequences and updates the row with match status.

    Args:
        row (pd.Series): Current row of the DataFrame.
        ref_rss_dict (dict): Dictionary containing reference RSS sequences.

    Returns:
        pd.Series: Row with added columns indicating match status.

    Raises:
        KeyError: If the segment or region combination is not found in the configuration, logs the error and raises an exception.
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    rss_variants = config['RSS_LAYOUT'].get(full, {}).keys()

    try:
        for rss_variant in rss_variants:
            if full in options:
                return check_ref_rss(row, ref_rss_dict, rss_variant)
    except KeyError as e:
        logger.error(f"Invalid segment or region combination: {e}")
        raise

    return row


def create_dict(row, separated_segments, config):
    """
    Creates a dictionary entry for the RSS sequence based on the current row of the DataFrame.
    Adds the segment, region, and function to the separated_segments dictionary.

    Args:
        row (pd.Series): Current row of the DataFrame.
        separated_segments (dict): Dictionary to store separated segments and their RSS sequences.

    Returns:
        pd.Series: The original row, unchanged.

    Raises:
        KeyError: If the segment or region combination is not found in the configuration, logs the error and raises an exception.
    """
    query, region, segment, function = row[[
        "Old name-like", "Region", "Segment", "Function"]]
    combi = region + segment
    rss_variants = config['RSS_LAYOUT'].get(combi, {}).keys()

    try:
        for rss_variant in rss_variants:
            rss_sequence = fetch_sequence(row, segment, rss_variant, config)
            add_segment(query, f"{segment}_{rss_variant}", region,
                        separated_segments, rss_sequence, function)
    except KeyError as e:
        logger.error(f"Invalid segment or region combination: {e}")
        raise

    return row


def create_RSS_files(df, RSS_directory, config):
    """
    Generates RSS files from the DataFrame and saves them to the specified directory.
    Iterates over the DataFrame to create a dictionary of separated segments, then writes them to FASTA files.

    Args:
        df (pd.DataFrame): DataFrame containing the information needed to create RSS files.
        RSS_directory (Path): Path to the directory where RSS files will be saved.

    Raises:
        OSError: If the directory creation or file writing fails, logs the error and raises an exception.
    """
    if not RSS_directory.exists():
        separated_segments = {}

        try:
            df = df.apply(lambda row: create_dict(
                row, separated_segments, config), axis=1)
            write_fasta_file(separated_segments, RSS_directory)
        except OSError as e:
            logger.error(f"Failed to create RSS files: {e}")
            raise


def create_all_RSS_meme_files(cwd, df, config):
    """
    Generates all necessary RSS and MEME files for the analysis.
    Creates directories and files for reference, novel, and combined RSS sequences, and runs MEME to generate motifs.

    Args:
        cwd (Path): Path object representing the current working directory.
        df (pd.DataFrame): DataFrame containing information about novel sequences.

    Returns:
        Path: Path to the directory containing the reference MEME motifs.

    Raises:
        OSError: If the directory creation or file writing fails, logs the error and raises an exception.
    """
    base = cwd / "RSS"
    df_100 = pd.read_excel(cwd / 'annotation' / 'annotation_report_100%.xlsx')
    combined = pd.concat([df_100, df], axis=0)
    datasets = [df_100, df, combined]
    rss_filenames = ["reference_RSS", "new_RSS", "combined_RSS"]
    meme_directories = ["reference_meme", "new_meme", "complete_meme"]

    try:
        for dataset, rss_filename, meme_directory in zip(datasets, rss_filenames, meme_directories):
            create_RSS_files(dataset, base / rss_filename, config)
            create_meme_directory(base / meme_directory,
                                  base / rss_filename, config)
    except OSError as e:
        logger.error(f"Failed to create RSS or MEME files: {e}")
        raise

    return base / meme_directories[0]


def update_df(df, ref_rss_dict, config, options):
    """
    Updates the DataFrame with base RSS components and checks against reference RSS sequences.
    Adds columns for heptamer and nonamer sequences and their match status.

    Args:
        df (pd.DataFrame): DataFrame to be updated.
        ref_rss_dict (dict): Dictionary containing reference RSS sequences.

    Returns:
        pd.DataFrame: Updated DataFrame with additional columns.

    Raises:
        KeyError: If the segment or region combination is not found in the configuration, logs the error and raises an exception.
    """
    try:
        df = df.apply(lambda row: add_base_rss_parts(row, config), axis=1)
        df = df.apply(lambda row: apply_check_ref_rss(
            row, ref_rss_dict, config, options), axis=1)
    except KeyError as e:
        logger.error(f"Failed to update DataFrame: {e}")
        raise

    return df


def create_rss_excel_file(cwd, final_df):
    """
    Process the DataFrame to remove non-best overlapping rows and export
    the remaining data into separate Excel files for 'Novel' and 'Known' statuses.
    """
    novel_df = final_df.query("Status == 'Novel'")
    known_df = final_df.query("Status == 'Known'")
    for df, filename in zip([novel_df, known_df], ["annotation_report", "annotation_report_100%"]):
        wrtie_rss_excel_file(cwd, df, filename)


def wrtie_rss_excel_file(cwd, df, filename):
    """
    Creates an extended Excel file with the updated DataFrame.
    Merges the original annotation report with the new RSS data and saves it as an Excel file.

    Args:
        cwd (Path): Path object representing the current working directory.
        df (pd.DataFrame): Updated DataFrame with additional RSS data.
        filename (str): Base name of the file to be created.

    Raises:
        OSError: If the file creation fails, logs the error and raises an exception.
    """
    try:
        logger.info(f"Generating {filename}_plus.xlsx!")
        df.to_excel(cwd / 'annotation' /
                          f'{filename}_plus.xlsx', index=False)
    except OSError as e:
        logger.error(f"Failed to create Excel file: {e}")
        raise


def check_if_exists(filename):
    """
    Checks if a file exists at the specified path.
    Logs an error and exits the program if the file is not found.

    Args:
        filename (str or Path): Path to the file to check.

    Returns:
        str or Path: The filename if it exists.

    Raises:
        SystemExit: If the file does not exist, logs the error and exits the program.
    """
    if Path(filename).exists():
        return filename
    else:
        logger.error(
            f"The {filename} file does not exist, closing application!")
        sys.exit()


def RSS_main():
    """
    Main function of the RSS creation script.
    Loads the configuration and DataFrames, generates RSS and MEME files, and updates the annotation report.
    The process involves extracting RSS sequences, generating motifs, comparing them against reference sequences,
    and creating an extended annotation report with additional RSS data.

    Steps:
        1. Load configuration settings and annotation reports.
        2. Generate RSS files and run MEME to identify motifs.
        3. Compare the identified motifs against reference sequences.
        4. Update the DataFrame with RSS data and match status.
        5. Save the extended annotation report as an Excel file.

    Raises:
        Exception: If any step fails, logs the error and raises an exception.
    """
    try:
        complete_df = pd.DataFrame()
        cwd = Path.cwd()

        config = load_config(cwd)  # new
        options = set(config.get("RSS_LAYOUT", {}).keys())  # new

        df1 = pd.read_excel(check_if_exists(
            cwd / 'annotation' / 'annotation_report.xlsx'))
        df2 = pd.read_excel(check_if_exists(
            cwd / 'annotation' / 'annotation_report_100%.xlsx'))
        ref_meme_directory = create_all_RSS_meme_files(cwd, df1, config)
        ref_rss_dict = make_reference_rss(ref_meme_directory, config)
        for df, filename in zip([df1, df2], ["annotation_report", "annotation_report_100%"]):
            df = update_df(df, ref_rss_dict, config, options)
            reference_df = combine_df(df, pd.read_excel(
                cwd / 'annotation' / f'{filename}.xlsx')).drop_duplicates()
            complete_df = pd.concat(
                [complete_df, reference_df], ignore_index=True)
        final_df = remove_overlapping_segments(complete_df)
        create_rss_excel_file(cwd, final_df)
    except Exception as e:
        logger.error(f"Failed in RSS_main: {e}")
        raise


if __name__ == '__main__':
    RSS_main()
