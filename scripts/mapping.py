import re
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from logger import custom_logger

"""
Used python packages:
    1. yaml
    2. pandas
    3. openpyxl
    4. biopython

Used CLI packages:
    1. minimap2
    2. bowtie2
    3. bowtie
    4. samtools
    5. bedtools
"""

# Method for logging the current states of the program.
logger = custom_logger(__name__)


class MappingFiles:
    """
    A class representing the file paths required for the mapping process.

    This class initializes the paths for various files needed in the mapping
    process, including the Bowtie/Bowtie2 database, SAM, BAM, and BED files.

    Args:
        prefix (str): A prefix used to name the files.
        index (Path): Path to the directory where index files are stored.
        beddir (Path): Path to the directory where BED files are stored.
    """

    def __init__(self, prefix, index, beddir):
        self.bowtie_db = index / f"{prefix}_index"
        self.sam = beddir / f"{prefix}.sam"
        self.bam = beddir / f"{prefix}_sorted.bam"
        self.bed = beddir / f"mapped_{prefix}.bed"


def get_sequence(line, fasta):
    """
    Extracts a sequence from a FASTA file based on coordinates in a BED file.

    This function takes a line from a BED file, extracts the sequence
    specified by the coordinates (start and stop) and reference, and
    returns the sequence.

    Args:
        line (list): A list representing a line from a BED file, 
        containing reference, start, and stop coordinates.
        fasta (str): Path to the FASTA file.

    Returns:
        str: The sequence extracted from the FASTA file based on 
        the coordinates in the BED file line.
    """
    reference, start, stop = line[0:3]
    sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    specific_sequence = sequences[reference].seq
    return str(specific_sequence[int(start): int(stop)])


def get_region_and_segment(name, cell_type):
    """
    Determines the region and segment from a potential segment name.

    This function identifies the region and segment in the provided name
    by checking for a part that starts with the specified cell type.

    Args:
        name (str): The name of a potential segment.
        cell_type (str): The type of cell (e.g., TR, IG) to look for in the name.

    Returns:
        str, str: The region and segment names, or "other" and "-" if not found.
    """
    prefix = [i for i in name.split("_") if i.startswith(cell_type)]
    if prefix:
        prefix = re.sub(r"[0-9]", "", prefix[0])
        if prefix.startswith(cell_type):
            return prefix[0:3], prefix[3]
    else:
        return "other", "-"


def parse_bed(file_path, accuracy, fasta, mapping_type, cell_type):
    """
    Parses a BED file and extracts relevant information for each entry.

    This function reads a BED file line by line, extracting information
    such as sequence, accuracy, file path, mapping type, region, and segment.
    The function returns a list of all the entries found in the BED file.

    Args:
        file_path (str): Path to the BED file.
        accuracy (int): Accuracy score between 1 and 100.
        fasta (str): Path to the reference FASTA file.
        mapping_type (str): The type of mapping tool used (`bowtie`, `bowtie2`, or `minimap2`).
        cell_type (str): The type of cell (e.g., TR, IG).

    Returns:
        list[list]: A nested list containing information for each entry in the BED file.
    """
    entries = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            region, segment = get_region_and_segment(line[3], cell_type)
            line.extend([get_sequence(line, fasta),
                         accuracy, str(file_path), mapping_type, region, segment, fasta])
            entries.append(line)
    return entries


def run_command(command):
    """
    Executes a shell command using subprocess.

    This function tries to run a specified shell command using subprocess.
    If the command fails, it logs an error with the exit code.

    Args:
        command (str): The shell command to execute.
    """
    try:
        subprocess.run(command, shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Command '{command}' failed with exit code {e.returncode}")


def make_bowtie2_command(acc, bowtie_db, rfasta, sam_file, threads):
    """
    Configures the Bowtie2 command for mapping sequences.

    This function adjusts several Bowtie2 parameters based on the provided accuracy score (`acc`):
    - `N`: The number of mismatches allowed in the seed alignment. It is set to 1 if `acc` is below 50, otherwise 0.
    - `L`: The seed length, calculated as an integer based on the accuracy. It is determined by the formula `L = 15 + (acc / 100) * 5`.
    - `score_min`: The minimum alignment score threshold, calculated based on the accuracy. It is set using the formula `score_min = L,0,(-0.1 + (acc / 100) * 0.08)`.

    The constructed command includes these parameters along with the Bowtie2 database, reference FASTA file, and the output SAM file paths.

    Args:
        acc (int): Accuracy score between 1 and 100.
        bowtie_db (str): Path to the Bowtie2 database index.
        rfasta (Path): Path to the reference FASTA file.
        sam_file (str): Path to the output SAM file.
        threads (int): Number of threads to use for the mapping process.

    Returns:
        str: A fully configured Bowtie2 command string.
    """
    N = 1 if acc < 50 else 0
    L = int(15 + (acc / 100) * 5)
    score_min_base = -0.1 + (acc / 100) * 0.08
    score_min = f"L,0,{score_min_base:.2f}"
    command = f"bowtie2 -p {threads} -N {N} -L {L} --score-min {score_min} -f -x {bowtie_db} -U {rfasta} -S {sam_file}"
    return command


def make_bowtie_command(acc, bowtie_db, rfasta, sam_file, threads):
    """
    Configures the Bowtie command for mapping sequences.

    This function adjusts the following Bowtie parameters based on the provided accuracy score (`acc`):
    - `v`: The number of allowed mismatches. It is set based on the accuracy score as follows:
        - 3 mismatches if `acc` <= 33
        - 2 mismatches if `acc` <= 66
        - 1 mismatch if `acc` < 100
        - 0 mismatches if `acc` = 100

    The constructed command includes these parameters along with the Bowtie database, reference FASTA file, and the output SAM file paths.

    Args:
        acc (int): Accuracy score between 1 and 100.
        bowtie_db (str): Path to the Bowtie database index.
        rfasta (str): Path to the reference FASTA file.
        sam_file (str): Path to the output SAM file.
        threads (int): Number of threads to use for the mapping process.

    Returns:
        str: A fully configured Bowtie command string.
    """
    mismatches = 3 if acc <= 33 else (2 if acc <= 66 else 1 if acc < 100 else 0)
    command = f"bowtie -p {threads} -v {mismatches} -m 1 -f -x {bowtie_db} {rfasta} -S {sam_file}"
    return command


def make_minimap2_command(acc, ffile, rfasta, sam_file, threads):
    """
    Configures the Minimap2 command for mapping sequences.

    The following parameters are used in the construction of the command:
    - `-a`: Output in SAM format.
    - `-m {acc}`: Sets the minimum alignment score required to output an alignment.
    - `-t {threads}`: Number of threads to use.

    The command maps the sequences in the query FASTA file (`ffile`) against the reference FASTA file (`rfasta`) and outputs the results to the specified SAM file.

    Args:
        acc (int): Accuracy score between 1 and 100.
        ffile (Path): Path to the query FASTA file.
        rfasta (Path): Path to the reference FASTA file.
        sam_file (Path): Path to the output SAM file.
        threads (int): Number of threads to use for the mapping process.

    Returns:
        str: A fully configured Minimap2 command string.
    """
    command = f"minimap2 -a -m {acc} -t {threads} {ffile} {rfasta} > {sam_file}"
    return command


def all_commands(files: MappingFiles, fasta_file, rfasta, acc, mapping_type, threads):
    """
    Generates all the necessary commands to run a sequence mapping process using the specified tool.

    Based on the mapping type (`bowtie`, `bowtie2`, or `minimap2`), the appropriate commands are generated:
    - If using `bowtie` or `bowtie2`, an index is first created using `bowtie-build` or `bowtie2-build`.
    - The sequences are then aligned, and the results are sorted and converted to BAM format using `samtools`.
    - Finally, the BAM file is converted to BED format using `bedtools`.

    Args:
        files (MappingFiles): An instance of the `MappingFiles` class containing paths to input and output files.
        fasta_file (str): Path to the query FASTA file containing sequences to be mapped.
        rfasta (str): Path to the reference FASTA file.
        acc (int): Accuracy score between 1 and 100.
        mapping_type (str): The mapping tool to be used (`bowtie`, `bowtie2`, or `minimap2`).
        threads (int): Number of threads to use for the mapping process.

    Returns:
        list: A list of shell command strings to execute the mapping process.
    """
    if mapping_type == "bowtie2":
        bowtie_command = make_bowtie2_command(
            acc, files.bowtie_db, rfasta, files.sam, threads)
    elif mapping_type == "minimap2":
        bowtie_command = make_minimap2_command(
            acc, fasta_file, rfasta, files.sam, threads)
    else:
        bowtie_command = make_bowtie_command(
            acc, files.bowtie_db, rfasta, files.sam, threads)

    commands = [
        f"{mapping_type}-build {fasta_file} {files.bowtie_db}",
        bowtie_command,
        f"samtools sort -o {files.bam} {files.sam}",
        f"samtools index {files.bam}",
        f"bedtools bamtobed -i {files.bam} > {files.bed}",
    ]

    if mapping_type == "minimap2":
        del commands[0]  # Remove the index building command for minimap2

    return commands


def make_df(all_entries):
    """
    Creates a pandas DataFrame from a list of entries.

    This function takes a list of entries (each representing information
    from a BED file line) and constructs a pandas DataFrame.

    Args:
        all_entries (list[list]): A nested list where each sublist contains 
        information for an entry.

    Returns:
        pd.DataFrame: A pandas DataFrame containing all entries.
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
    return df.reset_index(drop=True)


def run(indir, outdir, rfasta, beddir, acc, mapping_type, cell_type, threads):
    """
    Runs the mapping process for each FASTA file in the input directory.

    This function loops over the files in the specified input directory,
    creates the necessary index directories, runs alignment commands if
    the corresponding BED file is absent, and parses the resulting BED
    file entries.

    Args:
        indir (Path): Path to the input directory containing FASTA files.
        outdir (Path): Path to the output directory where indices are stored.
        rfasta (Path): Path to the reference FASTA file.
        beddir (Path): Path to the directory where BED files are stored.
        acc (int): Accuracy score between 1 and 100.
        mapping_type (str): The mapping tool to be used (`bowtie`, `bowtie2`, or `minimap2`).
        cell_type (str): The type of cell (e.g., TR, IG).
        threads (int): Number of threads to use for the mapping process.

    Yields:
        list[list]: A nested list containing parsed entries for each FASTA file.
    """
    for fasta in indir.glob("*.fasta"):
        prefix = fasta.stem
        index = outdir / prefix
        if mapping_type != "minimap2":
            create_directory(index)
        files = MappingFiles(prefix, index, beddir)
        if not files.bed.exists():
            for command in all_commands(files, fasta, rfasta, acc, mapping_type, threads):
                run_command(command)
        if files.bed.exists():
            entries = parse_bed(files.bed, acc, fasta, mapping_type, cell_type)
            logger.info(f"Parsed {len(entries)} entries from {files.bed}")
            yield entries
        else:
            logger.warning(f"Required file missing: {files.bed}")
            yield list()


def create_directory(location):
    """
    Creates a directory if it does not already exist.

    This function checks if a specified directory exists, and if not,
    it creates the directory (including any necessary parent directories).

    Args:
        location (str): Path to the directory to create.
    """
    Path(location).mkdir(parents=True, exist_ok=True)


def mapping_main(mapping_type, cell_type, input_dir, library, threads, start=100, stop=70):
    """
    Main function to run the mapping process for a specified tool.

    This function sets up the paths for input and output directories, and
    runs the mapping process for a range of accuracy scores. The results
    are accumulated into a pandas DataFrame, which is then returned.

    Args:
        mapping_type (str): The mapping tool to be used (`bowtie`, `bowtie2`, or `minimap2`).
        cell_type (str): The type of cell (e.g., TR, IG).
        input_dir (str): The input directory containing the region of interest FASTA files.
        library (str): Path to the reference library FASTA file.
        threads (int): Number of threads to use for the mapping process.
        start (int, optional): Starting accuracy score (default is 100).
        stop (int, optional): Stopping accuracy score (default is 70).

    Returns:
        pd.DataFrame: A DataFrame containing all entries from the mapping process.
    """
    cwd = Path.cwd()
    outdir = cwd / "mapping" / f"{mapping_type}_db"
    indir = cwd / input_dir
    rfasta = cwd / library
    all_entries = []
    for acc in range(start, stop - 1, -1):
        beddir = cwd / "mapping" / mapping_type / f"{acc}%acc"
        create_directory(beddir)
        for current_entry in run(indir, outdir, rfasta, beddir, acc, mapping_type, cell_type, threads):
            all_entries.extend(current_entry)
    df = make_df(all_entries)
    return df


if __name__ == "__main__":
    mapping_type = "minimap2"
    mapping_main(mapping_type)
