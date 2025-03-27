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
    if sam_file.empty:
        return [], first, second
    sam_file["MAPQ"] = sam_file["MAPQ"].astype(int)

    filtered_first = sam_file[sam_file["QNAME"].str.contains(first, na=False)] if first != "-" else pd.DataFrame(columns=["QNAME", "FLAG", "RNAME", "POS", "MAPQ"])
    filtered_second = sam_file[sam_file["QNAME"].str.contains(second, na=False)] if second != "-" else pd.DataFrame(columns=["QNAME", "FLAG", "RNAME", "POS", "MAPQ"])

    if filtered_first.empty and first != "-":
        file_log.warning(f"Flanking gene {first} not found in SAM file.")
        first = None
    if filtered_second.empty and second != "-":
        file_log.warning(f"Flanking gene {second} not found in SAM file.")
        second = None

    if not filtered_first.empty:
        filtered_first["POS"] = filtered_first["POS"].astype(int)
    if not filtered_second.empty:
        filtered_second["POS"] = filtered_second["POS"].astype(int)

    extraction_regions = []
    if not filtered_first.empty:
        for rname in filtered_first["RNAME"].unique():
            first_subset = filtered_first[filtered_first["RNAME"] == rname]
            first_subset = filter_best_mapq(first_subset)

            first_start = first_subset["POS"].min()
            first_end = first_subset["POS"].max()

            if filtered_second is not None and rname in filtered_second["RNAME"].values:
                second_subset = filtered_second[filtered_second["RNAME"] == rname]
                second_subset = filter_best_mapq(second_subset)

                second_start = second_subset["POS"].min()
                second_end = second_subset["POS"].max()

                start = min(first_start, second_start)
                end = max(first_end, second_end)
            elif second == "-" or filtered_second is not None:
                end = get_length_contig(sam, rname)
                start = first_start
                flag = int(first_subset["FLAG"].values[0])
                if flag == 16:
                    extraction_regions.append((rname, 0, start))
                else:
                    second = "-"
                    extraction_regions.append((rname, start, end))

                return extraction_regions, first, second
            else:
                continue
            extraction_regions.append((rname, start, end))


    elif not filtered_second.empty:
        for rname in filtered_second["RNAME"].unique():
            second_subset = filtered_second[filtered_second["RNAME"] == rname]
            second_subset = filter_best_mapq(second_subset)

            second_end = second_subset["POS"].max()
            extraction_regions.append((rname, 0, second_end))
    return extraction_regions, first, second


def extract(cwd: Union[str, Path], assembly_fasta: Union[str, Path], directory : Union[str, Path], first: str, second: str, sample: str, immuno_region: str):
    """
    Extracts a sequence from an assembly FASTA file based on flanking genes,
    and writes it to an output FASTA file.
    """
    if second == "":
        second = "-"
    log_data = {}
    sam = cwd / "tmp/mapped_genes" / assembly_fasta.with_suffix(".sam").name

    coords, first, second = get_positions_and_name(sam, first, second)

    if len(coords) > 1:
        coords = [coords[0]]
    if coords:
        contig_name = coords[0][0]
        start = coords[0][1]
        end = coords[0][2]

        output_path = directory / f"{sample}_{first}_{second}_{contig_name}_{immuno_region}.fasta"

        if len(coords) == 1:
            if not output_path.is_file():
                file_log.info(f"Extracting region: {first}, {second}, {contig_name}, {min(coords)}, {max(coords)}")

                cmd = f"samtools faidx {assembly_fasta} {contig_name}:{start}-{end}"
                result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                concatenated_sequence = "".join(result.stdout.strip().splitlines()[1:])
                with open(output_path, 'w') as file:
                    file.write(f">{output_path.stem}\n{str(Seq(concatenated_sequence))}")

            log_data = {
                "Region": immuno_region,
                "Contig": contig_name,
                "5'-contig": contig_name,
                "3'-contig": contig_name,
                "5_flanking_gene": first,
                "3_flanking_gene": second,
                "5_coords": int(start),
                "3_coords": int(end),
            }
        else:
            log_data = {
                "Region": immuno_region,
                "5'-contig": contig_name,
                "3'-contig": contig_name,
                "5_flanking_gene": first,
                "3_flanking_gene": second,
            }
    else:
        file_log.warning(f"No coordinates found for {first} and {second}. Region could not be extracted. {assembly_fasta.name}")
    return log_data if log_data else None

@log_error()
def region_main(flanking_genes: dict[list[str]], assembly_dir: Union[str, Path], threads: int):
    """
    Main function that processes SAM files to create region-specific assembly files
    based on flanking genes specified in the configuration.
    """
    cwd = Path.cwd()
    directory = cwd / "tmp/region"
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
                    output_json[assembly_name][immuno_region]["Extraction status"] = "Complete" if log_data.get("Contig") else "Fragmented"
                pbar.update(1)
    log_file = cwd / "annotation/broken_regions.json"
    with open(log_file, 'w') as f:
        json.dump(output_json, f, indent=4)

    if any(directory.iterdir()):
        file_log.info("Region extraction completed successfully")
    else:
        file_log.error("No regions were extracted")
        raise Exception("No regions extracted.")
