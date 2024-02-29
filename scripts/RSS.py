import re
import yaml
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq

CONFIG = None
OPTIONS = None


def load_config(cwd):
    global CONFIG
    global OPTIONS
    with open(cwd / 'config' / 'config.yaml', 'r') as file:
        CONFIG = yaml.safe_load(file)
        OPTIONS = set(CONFIG.get("RSS_LAYOUT", {}).keys())


def create_directory(location):
    """
    Create an directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """
    Path(location).mkdir(parents=True, exist_ok=True)


def write_fasta_file(dictionary, folder):
    """
    Writes FASTA files for each region and segment in the dictionary, 
    creating files named "{key}{region}.fasta" in the "folder" directory.
    Iterates over the dictionary to write the different RSS and spacer 
    from the segments. Each header is unique because it has a 
    incrementing number for every duplicate entry.
    The first key in the dictionary represents a region type and the 
    second is the segment, and its associated value RSS.

    Args:
        dictionary (dict): dictionary with all the different 
        regions and segments as nested keys and as value a 
        list with their RSS and spacer. 
        folder (Path): Path of the base directory to save the fasta
        files to.
    """
    create_directory(folder)
    for key, value in dictionary.items():
        for region, segments in value.items():
            ffile = folder / f"{key}{region}.fasta"
            with open(ffile, 'w') as f:
                for header, segment in segments.items():
                    for x, rss in enumerate(segment):
                        f.write(f">{header}-{x+1}\n{rss}\n")


def calculate_position(position, length, operation):
    if operation.endswith("plus"):
        return str(position + length)
    elif operation.endswith("minus"):
        return str(position - length)
    return str(position)


def rss_type(start, end, combi, rss_variant, strand):
    """
    Calculates a dictionary for every possible segment variation based 
    on the start and end coordinates, the type of segment and RSS type. 
    Based on the given segment and RSS the right coordinates are chosen
    from the dictionary and is returned as list. 


    Args:
        start (str): The start coordinate of the VDJ segment sequence.
        end (str): The end coordinate of the VDJ segment sequence.
        segment (str): Type of segment either a V,D or J.
        rss_variant (int): Type of RSS variant either 12 or 23.
        strand (str): Indicator is the sequence.
        is a "normal" or reverse.

    Returns:
        list: A list containing two integers representing the 
        start and end coordinates of the chosen RSS variant.
    """
    lengths = CONFIG['RSS_LENGTH']
    variants = CONFIG['RSS_LAYOUT'].get(combi, {})
    variant_length = lengths.get(str(rss_variant), 0)
    layout = variants.get(str(rss_variant), {}).get(strand, "")
    if layout == "end_plus":
        return [str(end), calculate_position(
            end, variant_length, layout)]
    elif layout == "start_minus":
        return [calculate_position(
            start, variant_length, layout), str(start)]


def fetch_sequence(row, segment, rss_variant):
    """
    Extracts the needed "query, fasta, start, end, strand, region" 
    from the current row of the df. The needed fasta file is opened and
    and based on the start and end coordinates the new coordinates are 
    determined and cut out of the sequence. If the strand is "-", 
    the sequence is made reverse complement. In the end the query 
    (the segment's name) is returned alongside the rss sequence.

    Args:
        row (Series): Current row from the the df.
        segment (str): Type of segment either a V,D or J.
        rss_variant (int): Type of RSS variant either 12 or 23.

    Returns:
        query (str): The name of the segment.
        rss (str): The sequence of the RSS based of the coordinates
        of segment.
    """
    query, fasta, start, end, strand, region = row[[
        'Old name-like', 'Path', 'Start coord', 'End coord', 'Strand',
        'Region']]
    with open(fasta, 'r') as fasta_file:
        record = SeqIO.read(fasta_file, 'fasta')
        rss_start, rss_stop = rss_type(
            start, end, f"{region}{segment}", rss_variant, strand)
        rss = record.seq[int(rss_start):int(rss_stop)]
        if strand == '-':
            rss = str(Seq(str(rss)).reverse_complement())
    return query, rss


def add_to_dict(query, dictionary, rss):
    """
    Method for adding the RSS to the dictionary (separated_segments).
    It checks if a certain entry exists if this is not the case a 
    new list is created and the RSS is added to list. If the
    entry exists the RSS is directly added to the list.

    Args:
        query (str): The name of the segment.
        dictionary (dict): Dictionary containing the different
        regions, segments and queries. Every query contains a
        list with corresponding RSS (or the RSS is added in this function).
        rss (str): The sequence of the RSS.
    """
    dictionary.setdefault(query, []).append(rss)


def get_mers(segment, rss, rss_variant):
    """
    Calculates a dictionary with the different 
    heptamer and the nonamer possibilities as list with two values.
    Based on the "segment"and "rss_variant" the right heptamer
    and nonamer list is chosen and returned. If no valid entry
    is found, ["", ""] is returned.


    Args:
        segment (str): Type of segment either a V,D or J.
        rss (str): The sequence of the RSS.
        rss_variant (int): Type of RSS variant either 12 or 23.

    Returns:
        list: List containing the heptamer and nonamer
        of a given segment and RSS variant.
    """
    mers = CONFIG["RSS_MERS"].get(str(rss_variant), [])
    mer1, mer2 = mers
    return [rss[0:mer1], rss[-mer2:]]


def add_to_row(row, val1, val2, rss_variant):
    """
    Add two columns "{rss_variant}_heptamer", "{rss_variant}_nonamer"
    to the current row from the df. The values from "val1" and "val2"
    are individually joined and checked which of them are seven bases long,
    this would be the heptamer and retained value is the nonamer.

    Args:
        row (Series): Current row from the df.
        val1 (list): List containing either seven or nine elements, 
        which represent either the heptamer or the nonamer of the RSS.
        val2 (list): List containing either seven or nine elements, 
        which represent either the heptamer or the nonamer of the RSS.
        rss_variant (int): Type of RSS variant either 12 or 23.


    Returns:
        row (Series): Current row from the df with the extra 
        "{rss_variant}_heptamer", "{rss_variant}_nonamer" columns. 
        Which contains the reference heptamer and nonamer for this row.

    """
    heptamer, nonamer = (''.join(val1), ''.join(val2)) if len(
        val1) == 7 else (''.join(val2), ''.join(val1))
    row[f"{rss_variant}_heptamer"], row[f"{rss_variant}_nonamer"] = heptamer, nonamer
    return row


def add_segment(query, segment, region, separated_segments, rss_sequence, function):
    """
    Creates an entry in the "separated_segments" dictionary if 
    region and/or does not exists, except when the function is 'P'.
    Then call the "add_to_dict()" function with the "query" and 
    "rss_sequence" to add the RSS to the dictionary.

    Args:
        query (str): The name of the segment.
        segment (str): Type of segment either a V,D or J.
        region (str): The type of V, D or J. 
        For example beta and gamma.
        separated_segments (dict): Dictionary containing the different
        regions, segments and queries (if they already).
        rss_sequence (str): The sequence of the RSS.
        function (str): The functional status of the segment,
        either "P" or "F/ORF"
    """

    if function != "P":
        add_to_dict(query, separated_segments.setdefault(
            region, {}).setdefault(segment, {}), str(rss_sequence))


def add_base_rss_parts(row):
    """
    Gets the "region, segment, function" from the current row of the 
    df and determines the value for the following new columns
    ["12_heptamer", "12_nonamer", "23_heptamer", "23_nonamer"] of this segment.
    It first creates these columns with black ("") values, then checks 
    if the type of segment it is. If it is a D then all four columns are 
    filled in because the D contains the two types of RSS. If it is not a D,
    it is checked if the segment is either V or something else. 
    When it is a V, 23 is assigned to it otherwise a 12. 
    In the end all four columns are returned but depending 
    on the type certain columns are filled in.

    Args:
        row (Series): Current row in the df.
        separated_segments (dict): Dictionary containing the different
        regions, segments and queries (if they already).

    Returns:
        (Series): Row with the added columns ["12_heptamer", 
        "12_nonamer", "23_heptamer", "23_nonamer"]. Depending on 
        the type of segment certain columns are filled in. If not not 
        filled in, the value is "".
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    row["12_heptamer"], row["12_nonamer"], row["23_heptamer"], row["23_nonamer"] = "", "", "", ""
    rss_variants = CONFIG['RSS_LAYOUT'].get(full, {}).keys()
    for rss_variant in rss_variants:
        query, rss_sequence = fetch_sequence(row, segment, rss_variant)
        val1, val2 = get_mers(segment, rss_sequence, rss_variant)
        add_to_row(row, val1, val2, rss_variant)
    return row


def run_meme(out, rss_file: Path, rss_variant):
    """
    Generate the meme suite command based on the provided "rss_file" 
    and output (out) file. The constructed command is then 
    ran with subprocess.

    Args:
        out (Path): Path of the output file for the meme command result.
        rss_file (Path): Path of the RSS input file.
    """
    multi_command = f"meme {rss_file} -o {out} -dna -mod zoops -nmotifs 1 -minw {rss_variant}"
    single_command = f"meme {rss_file} -o {out} -dna -mod anr -nmotifs 1 -minw {rss_variant}"
    amount = subprocess.run(f"cat {rss_file} | egrep '^>' | wc -l",
                            shell=True, capture_output=True)
    if int(amount.stdout) > 1:
        subprocess.run(multi_command, shell=True, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        subprocess.run(single_command, shell=True, stdout=subprocess.PIPE)


def get_reference_mers(regex_string, segment, rss_variant):
    """
    Gets a string (regex_string) which contains the 
    full RSS sequence from the meme result. Then with help of a
    regular expressions, certain characters within brackets or 
    any single characters are matched and returned as list. 
    Then a dictionary is constructed with all the possibilities. 
    Based on the "segment" type and "rss_variant" the right list is chosen.
    Finally the two values are extracted from the list and returned as 
    "val1" and "val2".


    Args:
        regex_string (str): The meme regex string that contains
        the RSS including spacer.
        segment (str): Type of segment either a V,D or J.
        rss_variant (int): Type of RSS variant either 12 or 23.



    Returns:
        val1(list): List containing the first of the extracted elements out of the 
        dictionary, this can either be a value which is 7 long or 9 long.
        val1(list): List containing the second of the extracted elements out of the 
        dictionary, this can either be a value which is 7 long or 9 long.
    """
    rss = re.findall(r'\[[^\]]*\]|.', regex_string)
    mers = CONFIG["RSS_MERS"].get(str(rss_variant), [])
    mer1, mer2 = mers
    return rss[0:mer1], rss[-mer2:]


def make_ref_dict(segment, ref_rss, val1, val2):
    """
    Checks ifs a "segments" exists in the "ref_rss" dictionary if this 
    is not the case its creates a empty dictionary, it also checks 
    this for the "heptamer" and "nonamer". Then checks if the length of 
    "val1" is equal to seven or nine. Depending of this result either 
    "val1" or "val2" is appointed to "heptamer" or "nonamer".

    Args:
        segment (str): Type of segment either a V,D or J.
        ref_rss (dict): Dictionary which contains the segment, heptamer,
        nonamer. 
        val1(list): List containing the first of the extracted 
        elements out of the dictionary, this can either be a value 
        which is 7 long or 9 long.
        val1(list): List containing the second of the extracted 
        elements out of the dictionary, this can either be a value 
        which is 7 long or 9 long.

    Returns:
        _type_: _description_
    """
    if len(val1) == 7:
        ref_rss.setdefault(segment, {}).setdefault(
            "heptamer", ''.join(val1))
        ref_rss.setdefault(segment, {}).setdefault(
            "nonamer", ''.join(val2))
    else:
        ref_rss.setdefault(segment, {}).setdefault(
            "heptamer", ''.join(val2))
        ref_rss.setdefault(segment, {}).setdefault(
            "nonamer", ''.join(val1))
    return ref_rss


def create_meme_directory(cwd, meme_directory, RSS_directory):
    meme_directory = cwd / "RSS" / meme_directory
    RSS_directory = cwd / "RSS" / RSS_directory
    create_directory(meme_directory)
    for rss_file in Path(RSS_directory).iterdir():
        stem = rss_file.stem
        rss_variant = stem.split("_")[-1]
        RSS_convert = CONFIG.get("RSS_LENGTH", {})
        out = meme_directory / stem
        meme = out / "meme.txt"
        if not meme.exists():
            run_meme(out, rss_file, RSS_convert[rss_variant])
    return meme_directory


def make_referece_rss(cwd, RSS_directory):
    """
    Loops over the folder "cwd (current directory the user is in) / RSS".
    It determines the segments name based on the name (stem) of the fasta file.
    After this it is checked if the different meme.txt files exists. If they does
    not exists the the run_meme() command is executed. When the meme.txt is 
    found/or created the whole RSS regular expression is fetched, the 
    important heptamer and nonamer extracted and appended to the
    "ref_rss" dictionary. When all the files from the "cwd / RSS"
    are processed, the "ref_rss" is returned.


    Args:
        cwd (Path): The current directory path where the user currently 
        is in.

    Returns:
        ref_rss (dict): Directory that contains all the reference RSS 
        heptamers and nonamers for all the regions and segments.
    """
    ref_rss = {}
    meme_directory = "reference_meme"
    meme_directory: Path = create_meme_directory(
        cwd, meme_directory, RSS_directory)
    for meme in meme_directory.iterdir():
        meme_text = meme / "meme.txt"
        command = f'cat {meme_text} | egrep -A2 "regular expression"'
        result = subprocess.run(command, shell=True,
                                capture_output=True, text=True)
        hits = result.stdout.replace(
            "-", "").replace("\t", "").strip().split("\n")
        hits = [hit for hit in hits if hit]
        split_stem = meme.stem.split("_")
        segment, rss_variant = split_stem[0][-1], int(split_stem[1])
        val1, val2 = get_reference_mers(hits[1], segment, rss_variant)
        make_ref_dict(meme.stem, ref_rss, val1, val2)
    return ref_rss


def check_ref_rss(row, ref_rss, rss_variant):
    """
    Compares the heptamer and nonamers from the current row with the 
    generated reference heptamers and nonamers. When the reference and 
    current one match, a True is filled in otherwise a False. 
    These values are stored in the column "{rss_variant}_{i}_matched", 
    based on the rss variant and if the i is equal to a heptamer or
    nonamer, the right columns is filled in.


    Args:
        row (Series): Current row of the df. 
        ref_rss (dict): Directory that contains all the reference RSS 
        heptamers and nonamers for all the regions and segments.
        rss_variant (int): Type of RSS variant either 12 or 23.

    Returns:
        _type_: _description_
    """
    region, segment = row[["Region", "Segment"]]
    ref = ref_rss[f"{region}{segment}_{rss_variant}"]
    for i in ["heptamer", "nonamer"]:
        ref_seq = ref[i]
        query_seq = row[f"{rss_variant}_{i}"]
        matches = re.findall(ref_seq, query_seq)
        row[f"{rss_variant}_ref_{i}"] = ref_seq
        if matches:
            row[f"{rss_variant}_{i}_matched"] = True
        else:
            row[f"{rss_variant}_{i}_matched"] = False
    return row


def combine_df(original_df, new_df):
    """
    Takes the original df created from the "annotation_report.xlsx"
    file and append the new df with the newly found heptamers and nonamers 
    for every potential segment and merge it on the original df.

    Args:
        original_df (DataFrame): The original df containing the original
        new found segments without the heptamers and nonamers of
        the new segments.
        new_df (DataFrame): New df frame containing the newly found 
        heptamers and nonamers of the newly found segments.

    Returns:
        combined_df (DataFrame): Newly created df containing the old df 
        information and the new df which contains the heptamers and 
        nonamers of the newly found segments.
    """
    columns_to_merge = [
        'Reference', 'Old name-like',
        '12_heptamer', '12_ref_heptamer', '12_heptamer_matched',
        '12_nonamer', '12_ref_nonamer', '12_nonamer_matched',
        '23_heptamer', '23_ref_heptamer', '23_heptamer_matched',
        '23_nonamer', '23_ref_nonamer', '23_nonamer_matched',
    ]
    combined_df = pd.merge(new_df, original_df[columns_to_merge],
                           on=['Reference', 'Old name-like'],
                           how='left')
    return combined_df


def apply_check_ref_rss(row, ref_rss):
    """
    Gets the segment from the row and subtracts the last 
    character, which is the segment itself. Then verifies
    the "segment" and "rss_variant" runs the "check_ref_rss()"

    Args:
        row (Series): Current row of the df.
        ref_rss (dict): Directory that contains all the reference RSS 
        heptamers and nonamers for all the regions and segments.
        rss_variant (int): Type of RSS variant either 12 or 23.

    Returns:
        row (Series): Row with extra added columns.
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    rss_variants = CONFIG['RSS_LAYOUT'].get(full, {}).keys()
    for rss_variant in rss_variants:
        if full in OPTIONS:
            return check_ref_rss(row, ref_rss, rss_variant)
    return row


def create_dict(row, separated_segments):
    region, segment, function = row[["Region", "Segment", "Function"]]
    full = region + segment
    rss_variants = CONFIG['RSS_LAYOUT'].get(full, {}).keys()
    for rss_variant in rss_variants:
        query, rss_sequence = fetch_sequence(row, segment, rss_variant)
        add_segment(query, f"{segment}_{rss_variant}", region,
                    separated_segments, rss_sequence, function)
    return row


def create_RSS_files(cwd, df, directory):
    separated_segments = {}
    df = df.apply(lambda row: create_dict(
        row, separated_segments), axis=1)
    rss_path = cwd / "RSS" / directory
    if not rss_path.exists():
        write_fasta_file(separated_segments, rss_path)


def create_all_RSS_meme_files(cwd, df):
    ref_RSS_file, new_RSS_file, combined_RSS = "reference_RSS", "new_RSS", "combined_RSS"
    df_100 = pd.read_excel(cwd / 'annotation' / 'annotation_report_100%.xlsx')
    combined = pd.concat([df_100, df], axis=0)
    create_RSS_files(cwd, df_100, ref_RSS_file)
    create_RSS_files(cwd, df, new_RSS_file)
    create_RSS_files(cwd, combined, combined_RSS)
    create_meme_directory(cwd, "new_meme", new_RSS_file)
    create_meme_directory(cwd, "complete_meme", combined_RSS)
    return ref_RSS_file


def RSS_main():
    """
    Main function of the RSS creation script. It first sets 
    a cwd (current directory the user is in) object 
    based current location. Based on this cwd some 
    input and output directories and files are set.
    Then loads the a DataFrame (df) based on RSS_report.xlsx and 
    initiates the first dictionary (separated_segment) which contains 
    the different region, segments, query (name of the potential new
    segment) and the respective RSS sequence. After this, new columns 
    are added to the df, such as "Segment", "Region" and the
    new RSS heptamers and nonamers columns. 
    Additionally is checked if the "RSS" exists. 
    If this is not the case it is created with the 
    write_fasta_file() function. Now the reference RSS heptamers and 
    nonamers dictionary is created to validate against. 
    In the end the validation between the potential new segments 
    heptamers and nonamers and the reference ones is done. This new df 
    is merged with the old df from "annotation_report.xlsx" to form and 
    new extended df and is written to "annotation_report_plus.xlsx".
    """
    cwd = Path.cwd()
    load_config(cwd)
    df = pd.read_excel(cwd / 'annotation' / 'RSS_report.xlsx')
    ref_RSS_file = create_all_RSS_meme_files(cwd, df)
    df = df.apply(lambda row: add_base_rss_parts(row), axis=1)
    ref_rss = make_referece_rss(cwd, ref_RSS_file)
    df = df.apply(
        lambda row: apply_check_ref_rss(
            row, ref_rss), axis=1)
    filename = "annotation_report"
    final_df = combine_df(df, pd.read_excel(
        cwd / 'annotation' / f'{filename}.xlsx')).drop_duplicates()
    final_df.to_excel(cwd / 'annotation' /
                      f'{filename}_plus.xlsx', index=False)


if __name__ == '__main__':
    RSS_main()
