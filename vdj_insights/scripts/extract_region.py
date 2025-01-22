"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import re
import subprocess
from typing import Union
from typing import Tuple
import json
from tqdm import tqdm

from Bio.Seq import Seq

from .util import make_dir, calculate_available_resources
from .property import log_error

from .logger import console_logger, file_logger


console_log = console_logger(__name__)
file_log = file_logger(__name__)


def write_seq(seq: str, output: Union[str, Path]):
    """
    Writes a sequence to an output FASTA file.
    """
    with open(output, 'w') as out_file:
        out_file.write(f">{output.stem}\n{seq}")
    file_log.info(f"Sequence written to {output}")


@log_error()
def get_best_coords(sam_list: list[str]):
    """
    Determines the best coordinates from a list of SAM file entries based on the lowest bitwise flag value.
    """
    bitwise_flag = float("inf")
    best = []
    for i in range(0, len(sam_list), 4):
        sublist = sam_list[i:i+4]
        sam_bitwise = int(sublist[0])
        if sam_bitwise < bitwise_flag:
            bitwise_flag, contig_name, best = sam_bitwise, sublist[1], sublist
    return best[2:], contig_name


def get_length_contig(sam_file: Union[str, Path], contig: str) -> str:
    cmd = f'grep "^@SQ" {sam_file} | awk \'$2 ~ /{contig}/ {{print $3}}\' | cut -d":" -f2'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    length_contig = result.stdout.strip("\n")
    return int(length_contig)


def get_positions_and_name(sam: Union[str, Path], first: str, second: str):
    """
    Extracts the positions and contig name from a SAM file for given flanking genes.
    Handles reverse strand mapping and assigns missing genes to the telomere if necessary.
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
                        coords.append(get_length_contig(sam, contig_name))
                    else:
                        coords.append(0)
                else:
                    if flag is not None and flag & 16:
                        coords.append(0)
                    else:
                        coords.append(get_length_contig(sam, contig_name))
                name.append(contig_name)
            else:
                if commands[i]:
                    line = subprocess.run(commands[i], shell=True, capture_output=True, text=True)
                    if line.returncode != 0:
                        file_log.error(f"Failed to run command: {commands[i]}")
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
                file_log.warning(f"Region is too short: {region_length} base pairs. No region extracted.")
                return [], []

        if not coords or not name:
            file_log.warning("No coordinates or contig names could be found. Region could not be extracted.")
            return [], []

        return coords, name

    except Exception as e:
        file_log.warning(f"Failed to get positions and name from SAM file: {e}")
        return [], []


def extract(cwd: Union[str, Path], assembly_fasta: Union[str, Path], directory : Union[str, Path], first: str, second: str, sample: str, haplotype: str):
    """
    Extracts a sequence from an assembly FASTA file based on flanking genes,
    and writes it to an output FASTA file.
    """
    outfile = directory / f"{sample}_{first}_{second}_{haplotype}.fasta"
    log_data = {}
    sam = cwd / "mapped_genes" / assembly_fasta.with_suffix(".sam").name

    coords, name = get_positions_and_name(sam, first, second)
    if coords:
        if len(set(name)) == 1:
            contig_name = name[0]
            if not outfile.is_file():
                file_log.info(f"Extracting region: {first}, {second}, {contig_name}, {min(coords)}, {max(coords)}")

                cmd = f"samtools faidx {assembly_fasta} {contig_name}:{min(coords)+1}-{max(coords)}"
                result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                concatenated_sequence = "".join(result.stdout.strip().splitlines()[1:])
                write_seq(str(Seq(concatenated_sequence)), outfile)

            log_data = {
                "flanking_regions": f"{first}-{second}",
                "Contig": contig_name
            }
        else:
            file_log.warning(
                f"Broken region detected, unable to create a valid region.: {assembly_fasta.name}, {first}, {second}, {name[0]}, {name[1]}, {min(coords)}, {max(coords)}")
            log_data = {
                "flanking_regions": f"{first}-{second}",
                "5_Contig": name[0],
                "3_Contig": name[1]
            }
    else:
        file_log.warning(f"No coordinates found for {first} and {second}. Region could not be extracted. {assembly_fasta.name}")
    return log_data if log_data else None


def clean_filename(filename: str):
    """
    Cleans the filename by removing common prefixes or suffixes that are not part of the key identifiers,
    but preserves legitimate accession codes like GCA or GCF.
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


def parse_name(filename: Union[str, Path]) -> Tuple[str, str, str]:
    """
    Parses the filename to extract the chromosome, sample (accession code or custom ID), and haplotype information,
    handling any order of parts and missing values.
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


def region_main(flanking_genes: list[str], assembly_dir: Union[str, Path], threads: int):
    """
    Main function that processes SAM files to create region-specific assembly files
    based on flanking genes specified in the configuration.
    """
    cwd = Path.cwd()
    directory = cwd / "region"
    make_dir(directory)

    assembly_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in Path(assembly_dir).glob(ext)]

    output_json = {}
    tasks = []

    for first, second in zip(*[iter(flanking_genes)] * 2):
        for assembly in assembly_files:
            chrom, sample, haplotype = parse_name(assembly)
            tasks.append((cwd, assembly, directory, first, second, sample, haplotype))

    max_jobs = calculate_available_resources(max_cores=threads, threads=4, memory_per_process=12)
    total_tasks = len(tasks)

    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {executor.submit(extract, *task): task for task in tasks}
        with tqdm(total=total_tasks, desc='Extracting regions', unit='task') as pbar:
            for future in as_completed(futures):
                log_data = future.result()
                if log_data:
                    assembly_name = futures[future][1].name
                    flanking_key = log_data["flanking_regions"]
                    if assembly_name not in output_json:
                        output_json[assembly_name] = {}
                    if flanking_key not in output_json[assembly_name]:
                        output_json[assembly_name][flanking_key] = {}
                    if "Contig" in log_data:
                        output_json[assembly_name][flanking_key]["5_Contig"] = log_data["Contig"]
                        output_json[assembly_name][flanking_key]["3_Contig"] = log_data["Contig"]
                        output_json[assembly_name][flanking_key]["assembly_type"] = "Complete"
                    else:
                        output_json[assembly_name][flanking_key]["5_Contig"] = log_data["5_Contig"]
                        output_json[assembly_name][flanking_key]["3_Contig"] = log_data["3_Contig"]
                        output_json[assembly_name][flanking_key]["assembly_type"] = "Fragmented"
                pbar.update(1)

    log_file = cwd / "broken_regions.json"
    with open(log_file, 'w') as f:
        json.dump(output_json, f, indent=4)

    if any(directory.iterdir()):
        file_log.info("Region extraction completed successfully")
    else:
        file_log.error("No regions were extracted")
        raise Exception("No regions extracted.")


if __name__ == '__main__':
    region_main(["TMEM121", "-", "RPIA", "LSP1P4", "GNAZ", "TOP3B"], "/mnt/nanopore/Jaimy_intern/test_one", 8)
