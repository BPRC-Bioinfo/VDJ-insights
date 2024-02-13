from pathlib import Path
import json
import subprocess
from Bio import SeqIO


def load_config(cwd):
    """
    Load the config json file "flanking.json" and return the dictionary
    from that file.

    Args:
        cwd (Path): Path object from the users current cwd (current directory
        the user is in).

    Returns:
        dict: Dictionary containing the chromosomes of interest and 
        their respective start and end flanking genes. 
    """
    with open(cwd / "config" / "flanking.json") as f:
        return json.load(f)


def write_seq(name, fasta, start, stop, out):
    """
    Get the header of the contig the region is placed ("location"), the start and 
    stop coordinates and the in ("fasta") and output ("out") files. 
    It reads the fasta file and fetches the needed contig out of the fasta file.
    The sequence and header are parsed and writen to the fasta file 
    "contig/{accession}_{region}_{hap}.fasta".

    Args:
        name (str): Name of the contig the current region is on.
        fasta (Path): Path of the input sam file.
        start (int): Start coordinate of the region of interest. 
        stop (int): End coordinate of the region of interest. 
        out (Path): Path of the output directory.
    """
    with open(fasta, 'r') as fasta_file, open(out, 'w') as out:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        record = record_dict[name]
        out.write(f">{record.id}\n{record.seq[start:stop]}")


def process(cwd, chrom, hap, fasta, config):
    """
    Loops over the config dictionary with the parsed chromosome. 
    It first parses the given sam file to construct a list 
    contaning the name of the contig ("name"), the start and stop coordinates 
    based on the start and stop coordinates. The name is stored in "name" 
    list and coordinates in "coords".
    Then is checked if the contig names are on the same conting and the 
    region of interest in on the same contig. If this is not the case a 
    message is given, that is was not able to construct a region. 
    Otherwise it parses the name of the contig 
    and min and max values from the coords list.

    Args:
        cwd (Path): Path object from the users current cwd (current directory
        the user is in).
        chrom (str): The chromosome abbreviation, that is used for this file.
        hap (str): The haplo type number of this same file.
        fasta (Path): Path of the input sam file.
        config (Dict): Dictionary contaning the chromosomes with their 
        region of interest and their start and stop flanking genes.
    """
    for region, value in config[chrom].items():
        coords, name = list(), list()
        for _, flank in value.items():
            awk = "awk '{if($1 !~ /^@/ && $6 !~ /\*/){print $3, $4, $4 + length($10) - 1}}'"
            command = f"cat {fasta} | egrep '{flank}' | {awk}"
            line = subprocess.run(command, shell=True,
                                  capture_output=True, text=True)
            splitted = line.stdout.strip().split()
            coords.extend([int(i) for i in splitted[1:]]
                          ), name.append(splitted[0])
        if len(set(name)) == 1:
            ffile = cwd / "converted" / "gfatofasta" / \
                f"{chrom}_{hap}.p.fasta"
            directory = cwd / "contig2"
            outfile = directory / f"EAW_{region}_{hap}.fasta"
            write_seq(name[0], ffile, min(coords), max(coords), outfile)
        else:
            print("Broken contig, can't create region!")


def main():
    """
    Main function of this script. I takes a sam file of the region and 
    creates the different region fasta files based on the flanking genes, 
    located in "config/flanking.json". Then parses the
    chromosome abbreviation and haplotype based from the name of the file.
    Finally the sam file is processed and the region files are made.
    """
    print(snakemake.wildcards)
    cwd = Path.cwd()
    fasta = Path(snakemake.input[0])
    inter = fasta.stem.split("_")
    chrom, hap = inter[0], inter[2].split(".")[0]
    config = load_config(cwd)
    process(cwd, chrom, hap, fasta, config)


if __name__ == "__main__":
    main()
