import re
import logging
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path


# Method for logging current states of the program.
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class MappingFiles:
    """
    Creating a class for all the needed files for a certain mapping type 
    based on a certain prefix, which is either
    "bowtie" or "bowtie2" or "minimap2".

    Parameters:
        -

    Returns:
        - 
    """

    def __init__(self, prefix, index, beddir):
        self.bowtie_db = index / f"{prefix}_index"
        self.sam = beddir / f"{prefix}.sam"
        self.bam = beddir / f"{prefix}_sorted.bam"
        self.bed = beddir / f"mapped_{prefix}.bed"


def get_sequence(line, fasta):
    """
    Takes a line from a bedfile, which is an list and a fasta file.
    It takes the first three elements from the line.
    Based on the coordinates from the bed file line list and the
    fasta file the sequence is cut out with SeqIO.


    Args:
        line (list): list of a bedfile line containing reference,
        start and stop.
        fasta (str): Path of the fasta file location.

    Returns:
        str: Sequence of a based on the coordinates from the bed file line.
    """
    reference, start, stop = line[0:3]
    sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    specific_sequence = sequences[reference].seq
    return str(specific_sequence[int(start): int(stop)])


def get_region_and_segment(name):
    """
    Takes in a name of a potential segment and find the region and
    segment name in that. The name must contain a part which start 
    with eiter a "TR" or "LOC". If a part starts with TR the region 
    and segment is returned otherwise return "LOC" and "-"

    Args:
        name (str): Potential name of a segment.

    Returns:
        str, str: Return either a region name, segment name as string 
        or LOC and - as str.
    """
    prefix = [i for i in name.split("_") if i.startswith(("TR", "LOC"))][0]
    prefix = re.sub(r"[0-9]", "", prefix)
    if prefix.startswith("TR"):
        return prefix[0:3], prefix[3]
    else:
        return "LOC", "-"


def parse_bed(file_path, accuracy, fasta, mapping_type):
    """

    Reads the content of a bed file and loops over the content line by 
    line to get all the neccessery information such as the sequence, 
    the used score/accuracy from 1 to 100, bed_file path it self, 
    the mapping type this can either be "bowtie" or "bowtie2" or "minimap2", region, 
    segment and the path of the fasta file.
    It returns a list with all the different found entries of the bedfile

    Args:
        file_path (str): The path of the bedfile 
        accuracy (int): Score between 1 and 100
        fasta (str): The path of a fasta file
        mapping_type (str): The type of mapping that is being used, 
        either "bowtie", "bowtie2" or "minimap2".

    Returns:
        entries (list[list]): A nested list which contains all the 
        information from this bedfile.
    """
    entries = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            region, segment = get_region_and_segment(line[3])
            line.extend([get_sequence(line, fasta),
                        accuracy, str(file_path), mapping_type, region, segment, fasta])
            entries.append(line)
    return entries


def run_command(command):
    """
    Tries to run a certain command with subprocess, if this fails it 
    throws a logging error saying failed with exit code.

    Args:
        command (str): A string containg the command.
    """
    try:
        subprocess.run(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Command '{command}' failed with exit code {e.returncode}")


def make_bowtie2_command(acc, bowtie_db, rfasta, sam_file):
    """
    Configures the bowtie2 command.
    It adjusts some parameters based on a given acc.
    - N: is either a 1 if "acc" is below 50, otherwise 0.
    - L: Calculates an integer seed length based on "acc".
    - score_min: Is a minimum score threshold for alignments,
    based on the "acc".
    In the end the command is returned with the right parameters, 
    bowtie_db and files.


    Args:
        acc (int): Score between 1 and 100
        bowtie_db (str): Bowtie database path
        rfasta (Path): Reference fasta file path.
        sam_file (str): Sam file path

    Returns:
        str: A constructed bowtie2 command.
    """
    N = 1 if acc < 50 else 0
    L = int(15 + (acc / 100) * 5)
    score_min_base = -0.1 + (acc / 100) * 0.08
    score_min = f"L,0,{score_min_base:.2f}"
    command = f"bowtie2 -N {N} -L {L} --score-min {score_min} -f -x {bowtie_db} -U {rfasta} -S {sam_file}"
    return command


def make_bowtie_command(acc, bowtie_db, rfasta, sam_file):
    """
    Configures the bowtie command.
    It adjust some parameters based on the given "acc".
    - v: These are the amount of mismatches, 3 if acc <= 33, 2 in acc <= 66 and 1 < 100.
    In the end the bowtie command is constructed and returned with the right parameters, 
    bowtie_db and files


    Args:
        acc (int): Score between 1 and 100
        bowtie_db (str): Bowtie database path
        rfasta (str): Reference fasta file path.
        sam_file (str): Sam file path


    Returns:
        str: A constructed bowtie command.
    """
    mismatches = 3 if acc <= 33 else (2 if acc <= 66 else 1 if acc < 100 else 0)
    command = f"bowtie -v {mismatches} -m 1 -f -x {bowtie_db} {rfasta} -S {sam_file}"
    return command


def make_minimap2_command(acc, ffile, rfasta, sam_file):
    """
    Constructs the minimap2 command. Based on the accuracy (acc) / 
    score, fasta file (ffile), reference fasta file (rfasta) and sam file.

    Args:
        acc (int): Score between 1 and 100
        ffile (Path): Region fasta file path.
        rfasta (Path): Reference fasta file path.
        sam_file (Path): Sam file path

    Returns:
        command (str): Constructed minimap2 command.
    """
    command = f"minimap2 -a -m {acc} -t 10 {ffile} {rfasta} > {sam_file}"
    return command


def all_commands(files: MappingFiles, fasta_file, rfasta, acc, mapping_type):
    """
    Creates all the neccesessary commands needed for a certain 
    accuracy/score combined with the fasta_file 
    (which is a fasta file of a region of interest). It also check what 
    the mapping type is to construct the right command and to remove 
    the build command because minimap2 does require it.  

    Args:
        files (MappingFiles): class containg all the input and output files.
        fasta_file (str): Fasta file path of a region of interest.
        rfasta (str): Reference fasta file path.
        acc (int): Score between 1 and 100
        mapping_type (str): The type of mapping that is being used, 
        either "bowtie", "bowtie2", or "minimap2".

    Returns:
        command (list): List of all mapping commands suited to run
        the right mapping variant.
    """
    if mapping_type == "bowtie2":
        bowtie_command = make_bowtie2_command(
            acc, files.bowtie_db, rfasta, files.sam)
    elif mapping_type == "minimap2":
        bowtie_command = make_minimap2_command(
            acc, fasta_file, rfasta, files.sam)
    else:
        bowtie_command = make_bowtie_command(
            acc, files.bowtie_db, rfasta, files.sam)
    commands = [
        f"{mapping_type}-build {fasta_file} {files.bowtie_db}",
        bowtie_command,
        f"samtools sort -o {files.bam} {files.sam}",
        f"samtools index {files.bam}",
        f"bedtools bamtobed -i {files.bam} > {files.bed}",
    ]
    if mapping_type == "minimap2":
        del commands[0]
    return commands


def make_df(all_entries):
    """
    Make a pandas dataframe (df) based on the given list of entries. 
    Duplicates based on a subset of ["name", "start", "stop"]
    are removed. Also the index of the df is reset.

    Args:
        all_entries (list[list]): A nested list which contains all the 
        information from this bedfile.

    Returns:
        _type_: _description_
    """
    headers = [
        "reference",
        "start",
        "stop",
        "name",
        "score",
        "strand",
        "sequence",
        "accuracy",
        "file",
        "tool",
        "region",
        "segment",
        "fasta-file"

    ]
    df = pd.DataFrame(all_entries, columns=headers)
    unique_combinations = df.drop_duplicates(subset=["name", "start", "stop"])
    return unique_combinations.reset_index(drop=True)


def run(cwd, indir, outdir, rfasta, beddir, acc, mapping_type):
    """
    Loop over the "indir" directory which contains all
    the region of interest fasta files. Create the index directories,
    runs alignment commands if BED file is absent, and parses BED file entries. 
    It also logs the number of parsed entries or warns the user if a 
    certain BED file is not present. In the end yield a list with all 
    parsed entries or an empty list per FASTA file.

    Args:
        cwd (Path): A path of the current folder location.
        indir (Path): A path to input directory.
        outdir (Path): A path to the base bowtie output directory.
        rfasta (Path): Reference fasta file path.
        beddir (Path): A path to the bed input directory.
        acc (int): Score between 1 and 100.
        mapping_type (str): The type of bowtie that is being used, 
        either bowtie or bowtie2.

    Yields:
        entries (list[list]): A nested list which contains all the 
        information from this bedfile.
    """

    for fasta in indir.glob("*.fasta"):
        prefix = fasta.stem
        index = outdir / prefix
        if mapping_type != "minimap2":
            create_directory(index)
        files = MappingFiles(prefix, index, beddir)
        if not files.bed.exists():
            for command in all_commands(files, fasta, rfasta, acc, mapping_type):
                run_command(command)
        if files.bed.exists():
            entries = parse_bed(files.bed, acc, fasta, mapping_type)
            logging.info(f"Parsed {len(entries)} entries from {files.bed}")
            yield entries
        else:
            logging.warning(f"Required file missing: {files.bed}")
            yield list()


def create_directory(location):
    """
    Create an directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """
    Path(location).mkdir(parents=True, exist_ok=True)


def mapping_main(mapping_type):
    """
    Main function to run the mapping script. It takes a mapping type in 
    as an argument to determine which type script to run. It created a 
    cwd path object and based on this it sets some input and output
    directories/files. It runs the run function
    for a accuracy/score (acc) of 1 to 100 to create
    all the different entries. Each run is saved in a directory with
    the following format: "{acc}%acc".
    When all entries are created and fetched a dataframe (df) is made
    and returned.

    Args:
    mapping_type (str): The type of mapping that is being used, 
    either "bowtie", "bowtie2", or "minimap2".


    Returns:
        df (DataFrame): A df containg all the entries.
    """
    cwd = Path.cwd()
    outdir = cwd / f"{mapping_type}_db"
    indir = cwd / "contig"
    rfasta = cwd / "library" / "retained.fasta"
    start, stop = 100, 0
    all_entries = []
    for acc in range(start, stop - 1, -1):
        beddir = cwd / mapping_type / f"{acc}%acc"
        create_directory(beddir)
        for current_entry in run(cwd, indir, outdir, rfasta, beddir, acc, mapping_type):
            all_entries.extend(current_entry)
    df = make_df(all_entries)
    return df


if __name__ == "__main__":
    mapping_type = "minimap2"
    mapping_main(mapping_type)
