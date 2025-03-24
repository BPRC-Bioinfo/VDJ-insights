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
import pandas as pd

from .util import make_dir, calculate_available_resources
from .property import log_error

from .logger import console_logger, file_logger


console_log = console_logger(__name__)
file_log = file_logger(__name__)


def get_length_contig(sam_file: Union[str, Path], contig: str) -> int:
    cmd = f'grep "^@SQ" {sam_file} | awk \'$2 == "SN:{contig}" {{print $3}}\' | cut -d":" -f2'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    length_contig = result.stdout.strip()
    return int(length_contig)


def filter_best_mapq(df: pd.DataFrame) -> pd.DataFrame:
    df = df.loc[df.groupby("RNAME")["MAPQ"].idxmax()]
    return df


def get_positions_and_name(sam: Union[str, Path], first: str, second: str) -> tuple[list[tuple[str, int, int]], str, str]:

    sam_file = pd.read_csv(sam, sep="\t", header=None, comment='@', dtype=str, usecols=[0, 1, 2, 3, 4], names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ"])
    sam_file["MAPQ"] = sam_file["MAPQ"].astype(int)

    df = sam_file[sam_file["QNAME"].str.contains(first, na=False) | sam_file["QNAME"].str.contains(second, na=False)]

    extraction_regions = []
    if not df.empty:
        for rname, group in df.groupby("RNAME"):
            first_subset = group[group["QNAME"].str.contains(first, na=False)]
            first_subset = filter_best_mapq(first_subset)

            second_subset = group[group["QNAME"].str.contains(second, na=False)]
            second_subset = filter_best_mapq(second_subset)

            if not first_subset.empty and not second_subset.empty:
                first_start = first_subset["POS"].min()
                first_end = first_subset["POS"].max()
                second_start = second_subset["POS"].min()
                second_end = second_subset["POS"].max()
                start = min(first_start, second_start)
                end = max(first_end, second_end)
                extraction_regions.append((rname, start, end, first, second))
            if not first_subset.empty and second_subset.empty:
                start = first_subset["POS"].min()
                end = get_length_contig(sam, rname)
                extraction_regions.append((rname, start, end, first,  "-"))
            if first_subset.empty and not second_subset.empty:
                end = second_subset["POS"].min()
                extraction_regions.append((rname, 0, end, "-", second))
    return extraction_regions


def extract(cwd: Union[str, Path], assembly_fasta: Union[str, Path], directory : Union[str, Path], first: str, second: str, sample: str, immuno_region: str):
    """
    Extracts a sequence from an assembly FASTA file based on flanking genes,
    and writes it to an output FASTA file.
    """
    if second == "":
        second = "-"
    log_data = {}
    sam = cwd / "mapped_genes" / assembly_fasta.with_suffix(".sam").name

    extraction_regions = get_positions_and_name(sam, first, second)
    print(extraction_regions)
    for contig, start, end, flanking_gene_one, flanking_gene_second in extraction_regions:
        print(contig, start, end, flanking_gene_one, flanking_gene_second)
        output_file = directory / f"{sample}_{flanking_gene_one}_{flanking_gene_second}_{immuno_region}.fasta"
        if not output_file.is_file():
            file_log.info(f"Extracting region: {flanking_gene_one}, {flanking_gene_second}, {contig}, {min(start)}, {max(end)}")
            cmd = f"samtools faidx {assembly_fasta} {contig}:{start}-{end}"
            result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            concatenated_sequence = "".join(result.stdout.strip().splitlines()[1:])
            with open(output_file, 'w') as file:
                file.write(f">{output_file.stem}\n{str(Seq(concatenated_sequence))}")

            log_data = {
                "Region": immuno_region,
                "Contig": contig,
                "5_flanking_gene": flanking_gene_one,
                "3_flanking_gene": flanking_gene_second,
                "5_coords": int(start),
                "3_coords": int(end),
            }

    return log_data if log_data else None

@log_error()
def region_main(flanking_genes: dict[list[str]], assembly_dir: Union[str, Path], threads: int):
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

    for region, extract_flanking_genes in flanking_genes.items():
        for assembly in assembly_files:
            tasks.append((cwd, assembly, directory, extract_flanking_genes[0], extract_flanking_genes[1], assembly.stem, region))
    max_jobs = calculate_available_resources(max_cores=threads, threads=4, memory_per_process=12)
    total_tasks = len(tasks)
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {executor.submit(extract, *task): task for task in tasks}
        with tqdm(total=total_tasks, desc='Extracting regions', unit='task') as pbar:
            for future in as_completed(futures):
                log_data = future.result()
                print(log_data)
                """

                if log_data:
                    assembly_name = futures[future][1].name
                    immuno_region = log_data["Region"]
                    if assembly_name not in output_json:
                        output_json[assembly_name] = {}
                    if immuno_region not in output_json[assembly_name]:
                        output_json[assembly_name][immuno_region] = {}
                    output_json[assembly_name][immuno_region]["5'-contig"] = log_data.get("5'-contig", None)
                    output_json[assembly_name][immuno_region]["3'-contig"] = log_data.get("3'-contig", None)
                    output_json[assembly_name][immuno_region]["5'-flanking_gene"] = log_data.get("5_flanking_gene", None)
                    output_json[assembly_name][immuno_region]["3'-flanking_gene"] = log_data.get("3_flanking_gene",None)
                    output_json[assembly_name][immuno_region]["5'-Coords"] = log_data.get("5_coords", None)
                    output_json[assembly_name][immuno_region]["3'-Coords"] = log_data.get("3_coords", None)
                    output_json[assembly_name][immuno_region]["Extraction status"] = "Complete" if log_data.get("5_flanking_gene", None) == log_data.get("3_flanking_gene",None) else "Fragmented"
                pbar.update(1)
                """
    log_file = cwd / "broken_regions.json"
    with open(log_file, 'w') as f:
        json.dump(output_json, f, indent=4)

    if any(directory.iterdir()):
        file_log.info("Region extraction completed successfully")
    else:
        file_log.error("No regions were extracted")
        raise Exception("No regions extracted.")
