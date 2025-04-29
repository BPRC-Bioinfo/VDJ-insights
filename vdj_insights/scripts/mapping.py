"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import re
import subprocess
import pandas as pd
from pathlib import Path

from Bio import SeqIO

from .util import make_dir, calculate_available_resources
from .logger import console_logger, file_logger


console_log = console_logger(__name__)
file_log = file_logger(__name__)


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


def parse_bed(file_path, fasta, mapping_type, cell_type):
    """
    Parses a BED file and extracts relevant information for each entry.

    This function reads a BED file line by line, extracting information
    such as sequence, accuracy, file path, mapping type, region, and segment.
    The function returns a list of all the entries found in the BED file.

    Args:
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
            line.extend([get_sequence(line, fasta), mapping_type, fasta])
            entries.append(line)
    return entries


def make_bowtie2_command(bowtie_db, rfasta, sam_file, threads):
    """
    Configures the Bowtie2 command for mapping sequences.

    This function adjusts several Bowtie2 parameters based on the provided accuracy score (`acc`):
    - `N`: The number of mismatches allowed in the seed alignment. It is set to 1 if `acc` is below 50, otherwise 0.
    - `L`: The seed length, calculated as an integer based on the accuracy. It is determined by the formula `L = 15 + (acc / 100) * 5`.
    - `score_min`: The minimum alignment score threshold, calculated based on the accuracy. It is set using the formula `score_min = L,0,(-0.1 + (acc / 100) * 0.08)`.

    The constructed command includes these parameters along with the Bowtie2 database, reference FASTA file, and the output SAM file paths.

    Args:
        bowtie_db (str): Path to the Bowtie2 database index.
        rfasta (Path): Path to the reference FASTA file.
        sam_file (str): Path to the output SAM file.
        threads (int): Number of threads to use for the mapping process.

    Returns:
        str: A fully configured Bowtie2 command string.
    """
    command = f"bowtie2 --end-to-end --very-sensitive -p {threads} --score-min L,0,-0.5 -f -x {bowtie_db} -U {rfasta} -S {sam_file}"
    #command = f"bowtie2 --end-to-end --very-sensitive -p {threads}  -f -x {bowtie_db} -U {rfasta} -S {sam_file}"
    return command


def make_bowtie_command(bowtie_db, rfasta, sam_file, threads):
    """
    Configures the Bowtie command for mapping sequences.
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
    #command = f"bowtie -a -n 2 -p {threads} -M 10 --strata -f -x {bowtie_db} {rfasta} -S {sam_file}"
    command = f"bowtie -k 5 -a --strata -p {threads} -f -x {bowtie_db} {rfasta} -S {sam_file}"
    #command = f"bowtie --best -n 2 -p {threads} -M 10 --strata -f -x {bowtie_db} {rfasta} -S {sam_file}"

    #command = f"bowtie --strata k 3 -p {threads} -f -x {bowtie_db} {rfasta} -S {sam_file}"
    return command


def make_minimap2_command(ffile, rfasta, sam_file, threads):
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
    command = f"minimap2 -a -m 70 -t {threads} {ffile} {rfasta} > {sam_file}"
    return command


def all_commands(files: MappingFiles, fasta_file, rfasta, mapping_type, threads):
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
        bowtie_command = make_bowtie2_command(files.bowtie_db, rfasta, files.sam, threads)
    elif mapping_type == "minimap2":
        bowtie_command = make_minimap2_command(fasta_file, rfasta, files.sam, threads)
    else:
        bowtie_command = make_bowtie_command(files.bowtie_db, rfasta, files.sam, threads)

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
        "tool",
        "fasta-file"
    ]
    df = pd.DataFrame(all_entries, columns=headers)
    return df.reset_index(drop=True)


def run_single_task(fasta, outdir, rfasta, mapping_type, cell_type, threads_per_process):
    """
    Runs the mapping process for each FASTA file in the input directory.

    This function loops over the files in the specified input directory,
    creates the necessary index directories, runs alignment commands if
    the corresponding BED file is absent, and parses the resulting BED
    file entries.

    Args:
        outdir (Path): Path to the output directory where indices are stored.
        rfasta (Path): Path to the reference FASTA file.
        beddir (Path): Path to the directory where BED files are stored.
        acc (int): Accuracy score between 1 and 100.
        mapping_type (str): The mapping tool to be used (`bowtie`, `bowtie2`, or `minimap2`).
        cell_type (str): The type of cell (e.g., TR, IG).
        threads (int): Number of threads to use for the mapping process.

    Returns:
        list[list]: A nested list containing parsed entries for each FASTA file.
    """
    try:
        beddir = outdir.parent / mapping_type
        make_dir(beddir)
        prefix = fasta.stem
        index = outdir / prefix
        if mapping_type != "minimap2":
            make_dir(index)
        files = MappingFiles(prefix, index, beddir)
        if not files.bed.exists():
            for command in all_commands(files, fasta, rfasta, mapping_type, threads_per_process):
                subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if files.bed.exists():
            entries = parse_bed(files.bed, fasta, mapping_type, cell_type)
            file_log.info(f"Parsed {len(entries)} entries from {files.bed}")
            return entries
        else:
            file_log.warning(f"Required file missing: {files.bed}")
            return []
    except Exception as e:
        file_log.error(f"Error in run_single_task for {fasta}: {e}")
        return []


def mapping_main(mapping_type, cell_type, input_dir, library, threads):
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
    outdir = cwd / "tmp/mapping" / f"{mapping_type}_db"
    indir = cwd / input_dir
    rfasta = library
    all_entries = []
    #tasks = list(indir.glob("*.fasta"))
    assembly_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in (indir).glob(ext)]

    total_tasks = len(assembly_files)

    max_jobs = calculate_available_resources(max_cores=threads, threads=2, memory_per_process=2)
    with ThreadPoolExecutor(max_workers=max_jobs) as executor:
        futures = [executor.submit(run_single_task, fasta, outdir, rfasta, mapping_type, cell_type, 2) for fasta in assembly_files]
        with tqdm(total=total_tasks, desc=f'Mapping library with {mapping_type}:', unit="file") as pbar:
            for future in as_completed(futures):
                try:
                    entries = future.result()
                    all_entries.extend(entries)
                except Exception as e:
                    file_log.error(f"Task resulted in an exception: {e}")
                pbar.update(1)
    df = make_df(all_entries)
    return df


if __name__ == "__main__":
    pass
