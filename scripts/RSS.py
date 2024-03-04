import re
import yaml
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq

# Global config settings set from config/config.yaml
CONFIG = None
OPTIONS = None


def load_config(cwd):
    """
    Reads the config file from config/config.yaml. 
    It stores the content as a dictionary in CONFIG. It also parses 
    all the potential region and segment combination in OPTIONS.

    Args:
        cwd (Path): Path object with the current working directory.
    """
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


def calculate_position(position, variant_length, operation):
    """
    Determines the start or end position of where do cut the RSS.
    It checks the operation and based on the value.
    It either return a end coordinate ("end_plus") or start coordinate
    ("start_minus). The coordinate is returned as string. 
    Args:
        position (int): Coordinate that either contains the 
        end or start of a RSS.
        variant_length (int): The length of needed RSS.
        operation (str): Indicator on which operation the perform.
        Can either be start_minus or end_plus.

    Returns:
        str: Needed coordinate as string
    """
    if operation.endswith("plus"):
        return str(position + variant_length)
    elif operation.endswith("minus"):
        return str(position - variant_length)
    return str(position)


def rss_type(start, end, combi, rss_variant, strand):
    """
    Calculates a list with coordinates on where the RSS is positioned. 
    It used either the start or end coordinate based on the value in 
    the config file. The rss length and layout are retreived from 
    the config.yaml file. Next the type of segment/region (combi) 
    combination is determined and the which RSS type it has. 
    Finally the way of cutting for this combi is fetched and 
    the coordinates are calculated and returned.



    Args:
        start (str): The start coordinate of the VDJ segment sequence.
        end (str): The end coordinate of the VDJ segment sequence.
        combi (str): Combination of V,D or J and region.
        rss_variant (int): Type of RSS variant either.
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
        rss_variant (int): Type of RSS variant.

    Returns:
        query (str): The name of the segment.
        rss (str): The sequence of the RSS based of the coordinates
        of segment.
    """
    fasta, start, end, strand, region = row[[
        'Path', 'Start coord', 'End coord', 'Strand', 'Region']]
    with open(fasta, 'r') as fasta_file:
        record = SeqIO.read(fasta_file, 'fasta')
        rss_start, rss_stop = rss_type(
            start, end, f"{region}{segment}", rss_variant, strand)
        rss = record.seq[int(rss_start):int(rss_stop)]
        if strand == '-':
            rss = str(Seq(str(rss)).reverse_complement())
    return rss


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


def get_mers(rss, rss_variant):
    """
    Calculates the heptamer and the nonamer based on the 
    RSS_MERS in the config file based on the rss_variant. 
    The  heptamer and nonamer is returned as a list.

    Args:
        rss (str): The sequence of the RSS.
        rss_variant (int): Type of RSS variant.

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
        rss_variant (int): Type of RSS variant.


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
        region (str): The destination where the VDJs are coming from.
        separated_segments (dict): Dictionary containing the different
        regions, segments and rss sequence.
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
    12_heptamer, 12_nonamer, 23_heptamer, 23_nonamer of this segment.
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
        regions, segments and rss sequence.

    Returns:
        (Series): Row with the added columns 12_heptamer, 
        12_nonamer, 23_heptamer, 23_nonamer. Depending on 
        the type of segment certain columns are filled in. If not not 
        filled in, the value is "".
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    row["12_heptamer"], row["12_nonamer"], row["23_heptamer"], row["23_nonamer"] = "", "", "", ""
    rss_variants = CONFIG['RSS_LAYOUT'].get(full, {}).keys()
    for rss_variant in rss_variants:
        rss_sequence = fetch_sequence(row, segment, rss_variant)
        val1, val2 = get_mers(rss_sequence, rss_variant)
        add_to_row(row, val1, val2, rss_variant)
    return row


def run_meme(out, rss_file, rss_variant):
    """
    Generate the meme suite commands based on the provided rss file,
    output file and rss_file. The constructed command is then 
    ran with subprocess. With a egrep and wc -l is calculated how many 
    headers are present in the file. If it is equal to 1 single command 
    is used otherwise the multi command.

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


def get_reference_mers(regex_string, rss_variant):
    """
    Gets a string (regex_string) which contains the 
    full RSS sequence from the meme result. Then with help of a
    regular expressions, certain characters within brackets or 
    any single characters are matched and returned as list. 
    With the rss mers and rss variant the heptamer and nonamer 
    are extracted from the regex list as separate lists
    mer1 and mer2. Finally both lists are returned.

    Args:
        regex_string (str): The meme regex string that contains
        the RSS including spacer.
        rss_variant (int): Type of RSS variant.

    Returns:
        list: List containing the first of the extracted elements out of the 
        regex list, this can either be a value which is 7 long or 9 long.
        list: List containing the second of the extracted elements out of the 
        regex list, this can either be a value which is 7 long or 9 long.
    """
    rss = re.findall(r'\[[^\]]*\]|.', regex_string)
    mers = CONFIG["RSS_MERS"].get(str(rss_variant), [])
    mer1, mer2 = mers
    return rss[0:mer1], rss[-mer2:]


def make_ref_dict(segment, ref_rss_dict, mer1, mer2):
    """
    Checks ifs a segments exists in the ref rss dict dictionary if this 
    is not the case its creates a empty dictionary, it also checks 
    this for the "heptamer" and "nonamer". Then checks if the length of 
    "mer1" is equal to seven or nine. Depending of this result either 
    "val1" or "val2" is appointed to "heptamer" or "nonamer".

    Args:
        segment (str): Type of segment either a V,D or J.
        ref_rss_dict (dict): Dictionary which contains the heptamers and
        nonamers for a the different segments. 
        val1(list): List containing the first of the extracted 
        elements out of the dictionary, this can either be a value 
        which is 7 long or 9 long.
        val1(list): List containing the second of the extracted 
        elements out of the dictionary, this can either be a value 
        which is 7 long or 9 long.

    Returns:
        ref_rss_dict (dict): Dictionary which contains the heptamers and
        nonamers for a the different segments. Plus this current 
        segments heptamer and nonamer.
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


def create_meme_directory(meme_directory, RSS_directory):
    """
    Creates a directory for the meme directory. Next it loops over the 
    fasta files (ffile) in rss directory and creates the needed 
    output directory path, based on the stem of the current ffile. 
    From the name of the ffile the rss type is also determined and 
    converted to the right length with rss length. 
    After this it is checked if the meme.txt file exists. If this is 
    not the case meme suite is run with the output directory, 
    the file itself and the rss length.

    Args:
        meme_directory (Path): Path of the meme directory.
        RSS_directory (Path): Path to the directory where the RSS fasta 
        ffiles are stored.
    """
    create_directory(meme_directory)
    for rss_file in Path(RSS_directory).iterdir():
        stem = rss_file.stem
        rss_variant = stem.split("_")[-1]
        RSS_convert = CONFIG.get("RSS_LENGTH", {})
        out = meme_directory / stem
        meme = out / "meme.txt"
        if not meme.exists():
            run_meme(out, rss_file, RSS_convert[rss_variant])


def make_reference_rss(ref_meme_directory):
    """
    Loops over the ref meme directory.
    It first determines the rss variant based on the name (stem) 
    of the ref meme subdirectory. If the file meme.txt is found in the 
    subdirectory the whole RSS regular expression is fetched, the 
    important heptamer and nonamer extracted and appended to the
    ref rss dictionary. When all the files from the ref meme sub directory
    are processed, the ref_rss_dict is returned.


    Args:
        ref_meme_directory (Path): The  directory of the reference meme 
        directory, where the needed reference meme sub directories 
        are located.

    Returns:
        ref_rss_dict (dict): Directory that contains all the reference RSS 
        heptamers and nonamers for all the regions and segments.
    """
    ref_rss_dict = {}
    for meme in ref_meme_directory.iterdir():
        meme_text = meme / "meme.txt"
        command = f'cat {meme_text} | egrep -A2 "regular expression"'
        result = subprocess.run(command, shell=True,
                                capture_output=True, text=True)
        hits = result.stdout.replace(
            "-", "").replace("\t", "").strip().split("\n")
        hits = [hit for hit in hits if hit]
        split_stem = meme.stem.split("_")
        rss_variant = int(split_stem[1])
        val1, val2 = get_reference_mers(hits[1], rss_variant)
        make_ref_dict(meme.stem, ref_rss_dict, val1, val2)
    return ref_rss_dict


def check_ref_rss(row, ref_rss_dict, rss_variant):
    """
    Compares the heptamer and nonamers from the current row with the 
    generated reference heptamers and nonamers. When the reference and 
    current one match, a True is filled in otherwise a False. 
    These values are stored in the column "{rss_variant}_{i}_matched", 
    based on the rss variant and if the i is equal to a heptamer or
    nonamer, the right columns is filled in.


    Args:
        row (Series): Current row of the df. 
        ref_rss_dict (dict): Directory that contains all the reference RSS 
        heptamers and nonamers for all the regions and segments.
        rss_variant (int): Type of RSS variant.

    Returns:
        Series: Current row with the filled in matched column.
    """
    region, segment = row[["Region", "Segment"]]
    ref = ref_rss_dict[f"{region}{segment}_{rss_variant}"]
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


def apply_check_ref_rss(row, ref_rss_dict):
    """
    Gets the segment from the row and subtracts the last 
    character, which is the segment itself. Then verifies
    the "segment" and "rss_variant" runs the "check_ref_rss()"

    Args:
        row (Series): Current row of the df.
        ref_rss_dict (dict): Directory that contains all the reference RSS 
        heptamers and nonamers for all the regions and segments.
        rss_variant (int): Type of RSS variant.

    Returns:
        row (Series): Row with extra added columns.
    """
    region, segment = row["Region"], row["Segment"]
    full = region + segment
    rss_variants = CONFIG['RSS_LAYOUT'].get(full, {}).keys()
    for rss_variant in rss_variants:
        if full in OPTIONS:
            return check_ref_rss(row, ref_rss_dict, rss_variant)
    return row


def create_dict(row, separated_segments):
    """
    Fetched the region, segment and function from the current row of the 
    df. Then it creates a combi of the region and segment to get 
    the complete variant. Of this variant the rss variant(s) are determined.
    Next loop over the rss variant(s) to get the rss sequence and name (query).
    Finally the sequence and name, segment and region are added to
    separated segments dict. 

    Args:
        row (Series): The current row of the df.
        separated_segments (dict): Dictionary containing the different
        regions, segments and rss sequence.

    Returns:
        row (Series): Current row of the df.
    """
    query, region, segment, function = row[[
        "Old name-like", "Region", "Segment", "Function"]]
    combi = region + segment
    rss_variants = CONFIG['RSS_LAYOUT'].get(combi, {}).keys()
    for rss_variant in rss_variants:
        rss_sequence = fetch_sequence(row, segment, rss_variant)
        add_segment(query, f"{segment}_{rss_variant}", region,
                    separated_segments, rss_sequence, function)
    return row


def create_RSS_files(df, RSS_directory):
    """
    Generates the separated segments dict if RSS directory is not 
    created. It also parses the separated segments dict to generate 
    the RSS fasta file.

    Args:
        df (DataFrame): Dataframe containing all the information needed 
        to create the RSS files.
        RSS_directory (Path): Path to the RSS directory the files 
        need to be stored.
    """
    if not RSS_directory.exists():
        separated_segments = {}
        df = df.apply(lambda row: create_dict(
            row, separated_segments), axis=1)
        write_fasta_file(separated_segments, RSS_directory)


def create_all_RSS_meme_files(cwd, df):
    """
    Set a base path, where the RSS and meme directories need to be 
    saved (cwd / RSS). Then create or fetch all the needed dfs, such as 
    the df containing information about novel sequences, non-novel sequences 
    and a combination of novel and non-novel sequences. 
    All the dfs are stored in a list.
    For the different dfs the suffixes are determined and stored in two list.
    Next it loops over all three lists with zip to create
    the RSS and meme files.


    Args:
        cwd (Path): Path object with the current working directory.
        df (DataFrame): Dataframe containing all the information of the
        novel sequences. 

    Returns:
        Path: The full path to the directory containing the novel meme 
        sequencing motifs.
    """
    base = cwd / "RSS"
    df_100 = pd.read_excel(cwd / 'annotation' / 'annotation_report_100%.xlsx')
    combined = pd.concat([df_100, df], axis=0)
    datasets = [df_100, df, combined]
    rss_filenames = ["reference_RSS", "new_RSS", "combined_RSS"]
    meme_directories = ["reference_meme", "new_meme", "complete_meme"]
    for dataset, rss_filename, meme_directory in zip(datasets, rss_filenames, meme_directories):
        create_RSS_files(dataset, base / rss_filename)
        create_meme_directory(base / meme_directory, base / rss_filename)
    return base / meme_directories[0]


def RSS_main():
    """
    Main function of the RSS creation script. It first sets 
    a cwd (current directory the user is in) object 
    based current location. Based on this cwd some 
    input and output directory and file paths are set.
    Then loads the a DataFrame (df) based on RSS_report.xlsx.
    First the needed whole RSSs are extracted and saved. Based on these 
    RSS, the meme motifs are generated. There are three types of 
    RSSs and meme motif being generated. The first is based 
    non-novel sequences, the second on novel-sequences and the last is 
    based on both novel and non-novel. The non-novel sequence meme motifs 
    are used to create the reference meme dictionary. After this, new columns 
    are added to the df, such as segment and region and the
    new RSS heptamers and nonamers columns. Now the reference RSS heptamers and 
    nonamers meme dictionary is created to use for validation.
    In the end the validation between the potential new segments 
    heptamers and nonamers and the reference ones is done. This new df 
    is merged with the old df from "annotation_report.xlsx" to form and 
    new extended df and is written to "annotation_report_plus.xlsx".
    """
    cwd = Path.cwd()
    load_config(cwd)
    df = pd.read_excel(cwd / 'annotation' / 'RSS_report.xlsx')
    ref_meme_directory = create_all_RSS_meme_files(cwd, df)
    df = df.apply(lambda row: add_base_rss_parts(row), axis=1)
    ref_rss_dict = make_reference_rss(ref_meme_directory)
    df = df.apply(
        lambda row: apply_check_ref_rss(
            row, ref_rss_dict), axis=1)
    filename = "annotation_report"
    final_df = combine_df(df, pd.read_excel(
        cwd / 'annotation' / f'{filename}.xlsx')).drop_duplicates()
    final_df.to_excel(cwd / 'annotation' /
                      f'{filename}_plus.xlsx', index=False)


if __name__ == '__main__':
    RSS_main()
