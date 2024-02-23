from pathlib import Path
import yaml
import subprocess
from Bio import SeqIO


def make_dir(dir):
    """
    Create an directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """

    Path(dir).mkdir(parents=True, exist_ok=True)


def make_record_dict(fasta):
    with open(fasta, 'r') as fasta_file:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        return record_dict


def load_config(cwd):
    """
    Load the config yaml file "config.yaml" and return the dictionary
    from that file.

    Args:
        cwd (Path): Path object from the users current cwd (current directory
        the user is in).

    Returns:
        dict: Dictionary containing the chromosomes of interest and 
        their respective start and end flanking genes. 
    """
    with open(cwd / "config" / "config.yaml") as f:
        return yaml.safe_load(f)["FLANKING"]


def write_seq(record_dict, name, start, stop, out):
    """
    Get the record_dict from the used fasta file (record_dict). The name of contig, start and 
    stop coordinates (name, start and stop) and output ("out") file. 
    It parses the name of the contig to record_dict and fetches the needed contig sequence out of the record_dict.
    The sequence and header are parsed and written to the fasta file 
    "contig/{accession}_{region}_hap{hap}.fasta".

    Args:
        record_dict (dict): Dictionary containing the SeqIO information
        for a certain fasta file.
        name (str): Name of the contig the current region is on.
        start (int): Start coordinate of the region of interest. 
        stop (int): End coordinate of the region of interest. 
        out (Path): Path of the output directory.
    """
    with open(out, 'w') as out:
        record = record_dict[name]
        out.write(f">{record.id}\n{record.seq[start:stop]}")


def get_positions_and_name(genes, sam, record_dict):
    """
    Initiate a list (coords) for the coordinates and list for the 
    contig name (contig_list). As well a global contig_name is set to 
    store the contig name. Then loop over genes dict to get the position 
    (start or end) and the flanking gene name. First there is checked if there 
    is no flaking gene name. If this is the case for a start position 
    a 0 appended to the coord list, if this is for the end position 
    the end of the contig is determined based on the contig name of
    the start contig and appended to coord list. If the flanking gene 
    name is present the start, end and name are retrieved from sam
    with help of awk and appended to the coord and name list.
    The coord and name list are returned.


    Args:
        genes (dict): Dictionary containing the start and end 
        flanking genes for a certain region.
        sam (Path): Path to the sam file.
        record_dict (dict): Dictionary containing the SeqIO information
        for a certain fasta file.

    Returns:
        coords (list): List containing the coordinates of the 
        flanking gene(s), begin or end of the contig.
        name (list): List containing the name of the contigs
    """
    coords, name = list(), list()
    contig_name = ""
    for position, flank in genes.items():
        if position == "start" and flank == "":
            coords.append(0)
        elif position == "end" and flank == "":
            record = record_dict[contig_name]
            coords.append(len(record.seq))
        else:
            awk = "awk '{if($1 !~ /^@/ && $6 !~ /\*/){print $3, $4, $4 + length($10) - 1}}'"
            command = f"cat {sam} | egrep '{flank}' | {awk}"
            line = subprocess.run(command, shell=True,
                                  capture_output=True, text=True)
            splitted = line.stdout.strip().split()
            contig_name = splitted[0]
            coords.extend([int(i) for i in splitted[1:]]
                          ), name.append(contig_name)
    return coords, name


def process(cwd, chrom, hap, sam, config):
    """
    Loops over the config dictionary with the parsed chromosome. 
    It first parses the given sam file to make_record_dict to construct a list 
    containing the name of the contig ("name"), the start and stop coordinates 
    based on the start and stop coordinates. 
    Then is checked if the contig names are on the same contig and the 
    region of interest in on the same contig. If this is not the case a 
    message is given, that is was not able to construct a region. 
    Otherwise it parses the name of the contig 
    and min and max values from the coords list.

    Args:
        cwd (Path): Path object from the users current cwd (current directory
        the user is in).
        chrom (str): The chromosome abbreviation, that is used for this file.
        hap (str): The haplo type number of this same file.
        sam (Path): Path of the input sam file.
        config (Dict): Dictionary containing the chromosomes with their 
        region of interest and their start and stop flanking genes.
    """
    ffile = cwd / "converted" / "gfatofasta" / \
        f"{chrom}_EAW_hap{hap}.fasta"
    record_dict = make_record_dict(ffile)
    for region, value in config[chrom].items():
        coords, name = get_positions_and_name(value, sam, record_dict)
        if len(set(name)) == 1:
            directory = cwd / "contig"
            make_dir(directory)
            outfile = directory / f"EAW_{region}_hap{hap}.fasta"
            write_seq(record_dict, name[0], min(coords), max(coords), outfile)
        else:
            print("Broken contig, can't create region!")


def main():
    """
    Main function of this script. I takes a sam file of the region and 
    creates the different region fasta files based on the flanking genes, 
    located in "config/config.yaml". Then parses the
    chromosome abbreviation and haplotype based from the name of the file.
    Finally the sam file is processed and the region files are made.
    """
    cwd = Path.cwd()
    sam = Path(cwd / snakemake.input[0])
    inter = sam.stem.split("_")
    chrom, hap = inter[0], inter[2][-1]
    config = load_config(cwd)
    process(cwd, chrom, hap, sam, config)


if __name__ == "__main__":
    main()
