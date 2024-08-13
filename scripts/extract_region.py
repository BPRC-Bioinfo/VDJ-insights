from pathlib import Path
import yaml
import subprocess
from Bio import SeqIO
from logger import custom_logger

"""
Used python packages:
    1. yaml
    2. biopython
"""
# Method for logging the current states of the program.
logger = custom_logger(__name__)


def make_dir(dir):
    """
    Create a directory if not existing.

    Args:
        location (str): Path of the directory to create.
    """
    try:
        Path(dir).mkdir(parents=True, exist_ok=True)
        logger.debug(f"Directory created or already exists: {dir}")
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")


def make_record_dict(fasta):
    try:
        with open(fasta, 'r') as fasta_file:
            record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
            return record_dict
    except Exception as e:
        logger.error(f"Failed to create record dictionary for {fasta}: {e}")
        raise


def load_config(cwd):
    """
    Load the config yaml file "config.yaml" and return the dictionary from that file.

    Args:
        cwd (Path): Path object from the user's current cwd (current directory the user is in).

    Returns:
        dict: Dictionary containing the chromosomes of interest and their respective start and end flanking genes.
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
    Write the sequence to the output file.

    Args:
        record_dict (dict): Dictionary containing the SeqIO information for a certain fasta file.
        name (str): Name of the region the current region is on.
        start (int): Start coordinate of the region of interest.
        stop (int): End coordinate of the region of interest.
        out (Path): Path of the output directory.
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
    Get the best coordinates from the sam list.

    Args:
        sam_list (list): List containing the bitwise flag(s), coordinates, and contig name(s), for a flanking gene.

    Returns:
        best[2:] (Slice): Slice of the best hit in the sam list containing only the coordinates.
        contig_name (str): The name of the contig of the best hit in the sam list.
    """
    try:
        bitwise_flag = float("inf")
        best = list()
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
    Get the positions and name from the sam file.

    Args:
        sam (Path): Path to the sam file.
        first (str): First flanking gene.
        second (str): Second flanking gene.
        record_dict (dict): Dictionary containing the SeqIO information for a certain fasta file.

    Returns:
        coords (list): List containing the coordinates of the flanking gene(s), begin or end of the region.
        name (list): List containing the name of the regions.
    """
    try:
        coords, name, best_coords, contig_name = list(), list(), list(), ""
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
    Extract the sequence and write it to the output file.

    Args:
        cwd (Path): Path object from the user's current cwd (current directory the user is in).
        assembly_fasta (Path): Path to the assembly fasta file.
        first (str): First flanking gene.
        second (str): Second flanking gene.
        sample (str): Sample identifier.
        haplotype (str): Haplotype identifier.
        config (dict): Dictionary containing the chromosomes with their region of interest and their start and stop flanking genes.
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
    Main function of this script. It takes a sam file of the region and creates the different region fasta files based on the flanking genes, located in "config/config.yaml".

    Args:
        flanking_genes (list): List of flanking genes.
        assembly_dir (str): Directory containing the assembly fasta files.
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
