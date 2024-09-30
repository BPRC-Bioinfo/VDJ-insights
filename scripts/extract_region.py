from pathlib import Path
import re
import subprocess
from typing import Tuple, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from util import make_dir
from property import log_error

from logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)



@log_error()
def make_record_dict(fasta: str | Path) -> dict:
    """
    Creates a dictionary of SeqIO records from a given FASTA file. 
    The dictionary allows for efficient access to sequences by their IDs.

    Args:
        fasta (str or Path): Path to the FASTA file.

    Returns:
        dict: A dictionary where the keys are sequence IDs and the values are SeqRecord objects.

    Raises:
        Exception: If the FASTA file cannot be read or parsed, logs the error and raises an exception.
    """
    with open(fasta, 'r') as fasta_file:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return record_dict


def write_seq(record_dict: dict, name: str, start: int, stop: int, output: str | Path):
    """
    Writes a sequence from a given record dictionary to an output FASTA file, 
    using specified start and stop coordinates.

    Args:
        record_dict (dict): Dictionary containing SeqRecord objects indexed by sequence ID.
        name (str): ID of the sequence to write.
        start (int): Start coordinate of the region to extract.
        stop (int): End coordinate of the region to extract.
        output (Path): Path to the output FASTA file.

    Raises:
        Exception: If the sequence cannot be written, logs the error and raises an exception.
    """
    with open(output, 'w') as out_file:
        record = record_dict[name]
        out_file.write(f">{output.stem}\n{record.seq[start:stop]}")
    console_log.info(f"Sequence written to {output}")


@log_error()
def get_best_coords(sam_list: list[str]):
    """
    Determines the best coordinates from a list of SAM file entries based on the lowest bitwise flag value.

    Args:
        sam_list (list): List of strings representing SAM file entries, including bitwise flags and coordinates.

    Returns:
        tuple: A tuple containing:
            - list: The best coordinates from the SAM list.
            - str: The contig name associated with the best coordinates.

    Raises:
        Exception: If the best coordinates cannot be determined, logs the error and raises an exception.
    """
    bitwise_flag = float("inf")
    best = []
    for i in range(0, len(sam_list), 4):
        sublist = sam_list[i:i+4]
        sam_bitwise = int(sublist[0])
        if sam_bitwise < bitwise_flag:
            bitwise_flag, contig_name, best = sam_bitwise, sublist[1], sublist
    return best[2:], contig_name



def get_positions_and_name(sam: str | Path, first: str, second: str, record_dict: Dict[str, SeqRecord]):
    """
    Extracts the positions and contig name from a SAM file for given flanking genes.
    Handles reverse strand mapping and assigns missing genes to the telomere if necessary.

    Args:
        sam (Path): Path to the SAM file.
        first (str): First flanking gene.
        second (str): Second flanking gene.
        record_dict (dict): Dictionary containing SeqRecord objects indexed by sequence ID.

    Returns:
        tuple: A tuple containing:
            - list: A list of coordinates for the flanking genes or telomeres.
            - list: A list of contig names associated with these coordinates.

    Warns:
        Logs a warning if positions and names cannot be extracted.
    """
    try:
        coords, name, best_coords, contig_name = [], [], [], ""
        awk = r"awk '{if($1 !~ /^@/ && $6 !~ /\*/){print $2, $3, $4, $4 + length($10) - 1}}'"
        genes = [first, second]
        commands = [
            f"egrep '{first}' {sam} | {awk}" if first != "-" else None,
            f"egrep '{second}' {sam} | {awk}" if second != "-" else None
        ]

        flag = None

        for i, gene in enumerate(genes):
            if gene == "-":
                if i == 0:
                    if flag is not None and flag & 16:
                        record = record_dict[contig_name]
                        coords.append(len(record.seq))
                    else:
                        coords.append(0)
                else:
                    if flag is not None and flag & 16:
                        coords.append(0)
                    else:
                        record = record_dict[contig_name]
                        coords.append(len(record.seq))
            else:
                if commands[i]:
                    line = subprocess.run(
                        commands[i], shell=True, capture_output=True, text=True)
                    if line.returncode != 0:
                        console_log.error(f"Failed to run command: {commands[i]}")
                        continue
                    sam_list = line.stdout.strip().split()
                    if sam_list:
                        flag = int(sam_list[0])
                        best_coords, contig_name = get_best_coords(sam_list)
                        coords.extend([int(coord) for coord in best_coords])
                        name.append(contig_name)

        if len(coords) == 2:
            region_length = abs(coords[1] - coords[0])
            if region_length < 100:
                console_log.warning(f"Region is too short: {region_length} base pairs. No region extracted.")
                return [], []

        if not coords or not name:
            console_log.warning("No coordinates or contig names could be found. Region could not be extracted.")
            return [], []

        return coords, name

    except Exception as e:
        console_log.warning(f"Failed to get positions and name from SAM file: {e}")
        return [], []


def extract(cwd: str | Path, assembly_fasta: str | Path, directory : str | Path, first: str, second: str, sample: str, haplotype: str):
    """
    Extracts a sequence from an assembly FASTA file based on flanking genes,
    and writes it to an output FASTA file.

    Args:
        cwd (Path): The current working directory.
        assembly_fasta (Path): Path to the assembly FASTA file.
        directory (Path):
        first (str): First flanking gene.
        second (str): Second flanking gene.
        sample (str): Sample identifier.
        haplotype (str): Haplotype identifier.
    Warns:
        Logs a warning if no region can be extracted instead of raising an exception.
    """
    outfile = directory / f"{sample}_{first}_{second}_{haplotype}.fasta"
    if not outfile.is_file():
        sam = cwd / "mapped_genes" / assembly_fasta.with_suffix(".sam").name

        record_dict = make_record_dict(assembly_fasta)
        coords, name = get_positions_and_name(sam, first, second, record_dict)

        if coords:
            if len(set(name)) == 1:
                console_log.info(f"Extracting region: {first}, {second}, {name[0]}, {min(coords)}, {max(coords)}")
                write_seq(record_dict, name[0], min(coords), max(coords), outfile)
            else:
                console_log.warning("Broken region detected, unable to create a valid region.")
        else:
            console_log.warning(f"No coordinates found for {first} and {second}. Region could not be extracted.")



def clean_filename(filename: str):
    """
    Cleans the filename by removing common prefixes or suffixes that are not part of the key identifiers, 
    but preserves legitimate accession codes like GCA or GCF.

    Args:
        filename (str): The original filename to clean.

    Returns:
        tuple: A tuple containing:
            - str: The cleaned filename without irrelevant prefixes or suffixes.
            - str: The accession code if found.
    """
    unwanted_terms = ["unmasked", "filtered", "trimmed", "masked"]
    accession_code_pattern = re.compile(r'(GCA|GCF|DRR|ERR)_?\d{6,9}(\.\d+)?')
    accession_match = accession_code_pattern.search(filename)
    if accession_match:
        accession_code = accession_match.group(0)
        filename = filename.replace(accession_code, "")
    else:
        accession_code = ""

    for term in unwanted_terms:
        filename = filename.replace(term, "")
    filename = re.sub(r'\.', '_', filename)
    filename = re.sub(r'_+', '_', filename).strip('_')

    return filename, accession_code


def parse_name(filename: str | Path) -> Tuple[str, str, str]:
    """
    Parses the filename to extract the chromosome, sample (accession code or custom ID), and haplotype information,
    handling any order of parts and missing values.

    Args:
        filename (Path): The file name to parse.

    Returns:
        tuple: A tuple containing:
            - str: Chromosome identifier (empty if not found).
            - str: Sample identifier (accession code, if found, or custom identifier).
            - str: Haplotype identifier (default to "hap1" if not found).
    """
    cleaned_stem, accession_code = clean_filename(filename.stem)
    name_part = cleaned_stem.split("_")
    chrom, sample, haplotype = "", accession_code, "hap1"
    chrom_pattern = re.compile(r'chr\d+|chr[XY]')
    haplotype_pattern = re.compile(r'hap\d+')

    for part in name_part:
        if chrom_pattern.match(part):
            chrom = part
        elif haplotype_pattern.match(part):
            haplotype = part
        elif not sample:
            sample = part

    return chrom, sample, haplotype


def region_main(flanking_genes: list[str], assembly_dir=""):
    """
    Main function that processes SAM files to create region-specific assembly files
    based on flanking genes specified in the configuration.

    Args:
        flanking_genes (list): List of flanking genes used to define regions.
        assembly_dir (str, optional): Directory containing the assembly FASTA files. Defaults to "".

    Raises:
        Exception: If the main region processing fails, logs the error and raises an exception.
    """
    cwd = Path.cwd()
    directory = cwd / "region"
    make_dir(directory)

    extensions = ["*.fna", "*.fasta", "*.fa"]

    for first, second in zip(*[iter(flanking_genes)]*2):
        assembly_files = [file for ext in extensions for file in Path(assembly_dir).glob(ext)]
        for assembly in assembly_files:
            chrom, sample, haplotype = parse_name(assembly)
            extract(cwd, assembly, directory, first, second, sample, haplotype)
    if any(directory.iterdir()):
        console_log.info("Region extraction completed successfully")
    else:
        console_log.error("No regions where extracted")
        raise
