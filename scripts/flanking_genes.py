import json
import pandas as pd
import subprocess
import os

contig_file = snakemake.input.contig_file
location_file = snakemake.input.location_file
output_dir = snakemake.output.output_dir

def get_matching_lines(file_path, pattern):
    """
    Retrieve specific fields from lines in a file that contain a given pattern.

    :param file_path: Path to the file.
    :param pattern: The pattern to search for in each line.
    :return: A list of strings with selected fields from matching lines.
    """
    with open(file_path, 'r') as file:
        return [f"{id}\t{contig}" for line in file if pattern in line for id, t1, contig, t2 in [line.strip().split("\t")]]

def process_locations(file_path):
    """
    Process the locations file to create a combinations dictionary.

    :param file_path: Path to the locations file.
    :return: A dictionary of combinations.
    """
    combinations = {}
    with open(file_path, "r") as f:
        next(f)  # Skip the first line
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

def make_df(combinations):
    tot = {}
    for key, genes in combinations.items():
        for gene in genes:
            matching_lines = get_matching_lines(contig_file, gene)
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


def write_contig_files(unicom):    
    for ttype, haplo, chrom, contigs in unicom.values.tolist():
        ttype = ttype.lower().replace(" & ", "-")
        for contig in contigs:
            input_file = f"converted/gfatofasta/chr{chrom}_EAW_hap{haplo}.p.fasta"
            output_file = f"contig/chr{chrom}_EAW_{ttype}_h{haplo}.fasta"
            command = f"awk '/{contig}/{{flag=1;print;next}}/^>/{{flag=0}}flag'  {input_file} >> {output_file}" 
            subprocess.call(command, shell=True)
               

if __name__ == "__main__":
    os.makedirs(output_dir, exist_ok=True)
    combinations = process_locations(location_file)
    df = make_df(combinations)
    unicom = unique(df)
    write_contig_files(unicom)
