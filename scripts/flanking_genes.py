import pandas as pd
import subprocess
import os

def make_new_QC(info_files):
    flanking_info = []
    for file in os.listdir(info_files):
        if "hap1" in file or "hap2" in file:
            accession = file.split('_')[1]
            file_string = f"{current_path}/{info_files}/{file}"
            with open(file_string, 'r') as f:
                for line in f:
                    flanking_info.append(line.strip())
    return flanking_info, accession


def get_matching_lines(flanking_genes, pattern):
    """
    Retrieve specific fields from lines in a file that contain a given pattern.

    :param file_path: Path to the file.
    :param pattern: The pattern to search for in each line.
    :return: A list of strings with selected fields from matching lines.
    """
    return [f"{id}\t{contig}" for line in flanking_genes if pattern in line for id, t1, contig, t2 in [line.strip().split("\t")]]

def process_locations(file_path):
    """
    Process the locations file to create a combinations dictionary.

    :param file_path: Path to the locations file.
    :return: A dictionary of combinations.
    """
    combinations = {}
    with open(file_path, "r") as f:
        next(f)
        chrid = None
        for line in f:
            line = line.strip()
            if line:
                if not line.startswith("Chr"):
                    chrid = line
                    combinations[chrid] = []
                else:
                    flanking = line.split()[-1]
                    combinations[chrid].append(flanking)
    return combinations

def make_df(combinations, flanking_genes):
    tot = {}
    for key, genes in combinations.items():
        for gene in genes:
            matching_lines = get_matching_lines(flanking_genes, gene)
            for i, line in enumerate(matching_lines[:2]):
                tot.setdefault(key, {}).setdefault(f"h{i+1}", []).append(line)
    rows = []
    for group, subdict in tot.items():
        for h_key, values in subdict.items():
            for value in values:
                id, contig = value.split("\t")
                gene, geneid, chromosome = id.split(":")
                rows.append([group, h_key[-1], gene, geneid.split("=")[-1], chromosome.split("=")[-1], contig])
    return pd.DataFrame(rows, columns=["Group", "Haplotype", "Gene", "GeneID", "Chromosome", "Contig"])


def unique(df):
    unique_combinations = df[['Group', 'Haplotype', 'Chromosome', 'Contig']].drop_duplicates()
    unique_combinations.reset_index(drop=True, inplace=True)
    unique_combinations = unique_combinations.groupby(["Group", "Haplotype", "Chromosome"])
    unique_combinations = unique_combinations["Contig"].apply(list).reset_index()
    return unique_combinations


def get_regions_files(unicom, accession):
    for ttype, haplo, chrom, contigs in unicom.values.tolist():
        if len(contigs) == 1:
            file_string = f"{current_path}/{bed_files}/chr{chrom}_{accession}_hap{haplo}.p_flank_aligned.bed"
            with open(file_string, "r") as f:
                for region in get_region(f):
                    yield region


def get_region(file):
    region = {}
    for line in file:
        contig, start, stop, gene = line.split()[0:4]
        if contig not in region:
            region[contig] = [start, stop]
        else:
            region[contig].extend([start, stop])

    for key, value in region.items():
        sorted_values = sorted(value, key=int)
        yield [key, sorted_values[1], sorted_values[2]]


def process_data(location_file, output_dir, info_files, bed_files):
    """
    Main processing function for genomic data.

    :param location_file: Path to the location file.
    :param output_dir: Directory to output the results.
    :param info_files: Directory containing QC files.
    :param bed_files: Directory containing BED files.
    :return: A dictionary of contig information.
    """
    os.makedirs(output_dir, exist_ok=True)
    flanking_genes, accession = make_new_QC(info_files)
    combinations = process_locations(location_file)
    df = make_df(combinations, flanking_genes)
    unicom = unique(df)
    
    contigs = {}
    for region in get_regions_files(unicom, accession):
        contigs[region[0]] = [region[1], region[2]]
        with open(os.path.join(output_dir, f"{region[0]}.bed"), "w") as bed:
            bed.write("\t".join(region))
    
    for ttype, haplo, chrom, contig in unicom.values.tolist():
        try:
            ttype = ttype.lower().replace(" & ", "-")
            ffile = f"converted/gfatofasta/chr{chrom}_{accession}_hap{haplo}.p.fasta"
            bedfile = os.path.join(output_dir, f"{contig[0]}.bed")
            output_region = os.path.join(output_dir, f"{accession}_{ttype}_hap{haplo}.fasta")
            command = f"bedtools getfasta -fi {ffile} -bed {bedfile} -fo {output_region}"
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {command}. Removing {output_region}!")
            os.remove(output_region)
            
    
    return contigs

if __name__ == "__main__":
    current_path = os.getcwd()
    location_file = snakemake.input.location_file
    output_dir = snakemake.output.output_dir
    info_files = snakemake.input.qc_folder
    bed_files = snakemake.input.flank_folder

    # Execute the main genomic data processing function
    contigs = process_data(location_file, output_dir, info_files, bed_files)