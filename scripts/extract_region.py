from pathlib import Path
import yaml
import subprocess
from Bio import SeqIO
from logger import custom_logger

"""
Used Python packages:
    1. yaml
    2. biopython
"""
# Method for logging the current states of the program.
logger = custom_logger(__name__)


def make_dir(dir):
    """
    Creates a directory and any necessary parent directories if they do not already exist.

    Args:
        dir (str or Path): Path of the directory to create.

    Raises:
        Exception: If the directory cannot be created, logs the error and raises an exception.
    """
    try:
        Path(dir).mkdir(parents=True, exist_ok=True)
        logger.debug(f"Directory created or already exists: {dir}")
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")
        raise


def make_record_dict(fasta):
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
    try:
        with open(fasta, 'r') as fasta_file:
            record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
            return record_dict
    except Exception as e:
        logger.error(f"Failed to create record dictionary for {fasta}: {e}")
        raise


def load_config(cwd):
    """
    Loads a configuration file (config.yaml) located in the 'config' directory of the current working directory.

    Args:
        cwd (Path): The current working directory.

    Returns:
        dict: Dictionary containing the configuration settings, including chromosomes of interest and their respective flanking genes.

    Raises:
        Exception: If the configuration file cannot be loaded, logs the error and raises an exception.
    """
    try:
        with open(cwd / "config" / "config.yaml") as f:
            config = yaml.safe_load(f)
            logger.info("Config file loaded successfully")
            return config
    except Exception as e:
        logger.error(f"Failed to load config file: {e}")
        raise


def write_seq(record_dict, name, start, stop, out):
    """
    Writes a sequence from a given record dictionary to an output FASTA file, 
    using specified start and stop coordinates.

    Args:
        record_dict (dict): Dictionary containing SeqRecord objects indexed by sequence ID.
        name (str): ID of the sequence to write.
        start (int): Start coordinate of the region to extract.
        stop (int): End coordinate of the region to extract.
        out (Path): Path to the output FASTA file.

    Raises:
        Exception: If the sequence cannot be written, logs the error and raises an exception.
    """
    try:
        with open(out, 'w') as out_file:
            record = record_dict[name]
            out_file.write(f">{out.stem}\n{record.seq[start:stop]}")
            logger.info(f"Sequence written to {out}")
    except Exception as e:
        logger.error(f"Failed to write sequence to {out}: {e}")
        raise


def get_best_coords(sam_list):
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
    try:
        bitwise_flag = float("inf")
        best = []
        for i in range(0, len(sam_list), 4):
            sublist = sam_list[i:i+4]
            sam_bitwise = int(sublist[0])
            if sam_bitwise < bitwise_flag:
                bitwise_flag, contig_name, best = sam_bitwise, sublist[1], sublist
        return best[2:], contig_name
    except Exception as e:
        logger.error(f"Failed to get best coordinates: {e}")
        raise


def get_positions_and_name(sam, first, second, record_dict):
    """
    Extracts the positions and contig name from a SAM file for given flanking genes.

    Args:
        sam (Path): Path to the SAM file.
        first (str): First flanking gene.
        second (str): Second flanking gene.
        record_dict (dict): Dictionary containing SeqRecord objects indexed by sequence ID.

    Returns:
        tuple: A tuple containing:
            - list: A list of coordinates for the flanking genes.
            - list: A list of contig names associated with these coordinates.

    Raises:
        Exception: If positions and names cannot be extracted, logs the error and raises an exception.
    """
    try:
        coords, name, best_coords, contig_name = [], [], [], ""
        for c, gene in enumerate([first, second]):
            awk = r"awk '{if($1 !~ /^@/ && $6 !~ /\*/){print $2, $3, $4, $4 + length($10) - 1}}'"
            command = f"cat {sam} |  egrep '{gene}' | {awk}"
            if c == 0 and not gene:
                coords.append(0)
            elif c == 1 and not gene:
                record = record_dict[contig_name]
                coords.append(len(record.seq))
            else:
                line = subprocess.run(command, shell=True,
                                      capture_output=True, text=True)
                sam_list = line.stdout.strip().split()
                if sam_list:
                    best_coords, contig_name = get_best_coords(sam_list)
                coords.extend([int(i) for i in best_coords]
                              ), name.append(contig_name)
        return coords, name
    except Exception as e:
        logger.error(f"Failed to get positions and name from SAM file: {e}")
        raise


def extract(cwd, assembly_fasta, directory, first, second, sample, haplotype, config):
    """
    Extracts a sequence from an assembly FASTA file based on flanking genes, 
    and writes it to an output FASTA file.

    Args:
        cwd (Path): The current working directory.
        assembly_fasta (Path): Path to the assembly FASTA file.
        first (str): First flanking gene.
        second (str): Second flanking gene.
        sample (str): Sample identifier.
        haplotype (str): Haplotype identifier.
        config (dict): Dictionary containing regions of interest and their associated flanking genes.

    Raises:
        Exception: If the sequence cannot be extracted, logs the error and raises an exception.
    """
    outfile = directory / f"{sample}_{first}_{second}_{haplotype}.fasta"
    if not outfile.is_file():
        try:
            sam = cwd / "mapped_genes" / assembly_fasta.with_suffix(".sam").name
            record_dict = make_record_dict(assembly_fasta)
            coords, name = get_positions_and_name(
                sam, first, second, record_dict)
            if coords:
                if len(set(name)) == 1:
                    logger.info(
                        f"Extracting region: {first}, {second}, {name[0]}, {min(coords)}, {max(coords)}")
                    write_seq(record_dict, name[0], min(
                        coords), max(coords), outfile)
                else:
                    logger.warning("Broken region, can't create region!")
        except Exception as e:
            logger.error(f"Failed to extract region: {e}")
            raise


def create_name(filename: Path):
    """
    Parses the filename to extract the chromosome, sample, and haplotype information.

    Args:
        filename (Path): The file name to parse.

    Returns:
        tuple: A tuple containing:
            - str: Chromosome identifier.
            - str: Sample identifier.
            - str: Haplotype identifier.
    """
    name_part = filename.stem.split("_")
    chrom, sample, haplotype = "", "", ""

    for part in name_part:
        if part.startswith("hap"):
            haplotype = part
        elif part.startswith("chr"):
            chrom = part
        else:
            sample = part

    return chrom, sample, haplotype


def region_main(flanking_genes, assembly_dir=""):
    """
    Main function that processes SAM files to create region-specific FASTA files 
    based on flanking genes specified in the configuration.

    Args:
        flanking_genes (list): List of flanking genes used to define regions.
        assembly_dir (str, optional): Directory containing the assembly FASTA files. Defaults to "".

    Raises:
        Exception: If the main region processing fails, logs the error and raises an exception.
    """
    try:
        cwd = Path.cwd()
        directory = cwd / "region"
        make_dir(directory)
        config = load_config(cwd)
        for first, second in zip(*[iter(flanking_genes)]*2):
            extensions = ["*.fna", "*.fasta", "*.fa"]
            fasta_files = [file for ext in extensions for file in Path(
                assembly_dir).glob(ext)]
            for assembly in fasta_files:
                chrom, sample, haplotype = create_name(assembly)
                extract(cwd, assembly, directory, first,
                        second, sample, haplotype, config)
        logger.info("Region extraction completed successfully")
    except Exception as e:
        logger.error(f"Failed in region_main: {e}")
        raise


if __name__ == "__main__":
    pass
