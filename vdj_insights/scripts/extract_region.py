"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import subprocess
from typing import Union
import json
from tqdm import tqdm

from Bio.Seq import Seq
import pandas as pd

from .util import make_dir, calculate_available_resources
from .property import log_error

from .logger import console_logger, file_logger
import warnings
import re

warnings.simplefilter(action='ignore', category=FutureWarning)

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def get_length_contig(sam_file: Union[str, Path], contig: str) -> int:
    cmd = f'grep "^@SQ" {sam_file} | awk \'$2 == "SN:{contig}" {{print $3}}\' | cut -d":" -f2'
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output = result.stdout.strip()
    if output.isdigit():
        return int(output)


def filter_best_ms(df: pd.DataFrame) -> pd.DataFrame:
    return df.loc[df.groupby("RNAME")["MS"].idxmax()]


def get_positions_and_name(sam: Union[str, Path], first: str, second: str) -> tuple[list[tuple[str, int, int]], str, str]:
    try:
        sam_file = pd.read_csv(sam, sep="\t", header=None, comment='@', dtype=str, usecols=[0, 1, 2, 3, 4, 12], names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "MS"])
        if second != "-":
            filterd_sam = sam_file[(sam_file["QNAME"].str.contains(first, na=False) | (sam_file["QNAME"].str.contains(second, na=False))) & (sam_file["FLAG"].astype(int) != 4)]
        else:
            filterd_sam = sam_file[(sam_file["QNAME"].str.contains(first, na=False)) & (sam_file["FLAG"].astype(int) != 4)]
        extraction_regions = []
        if not filterd_sam.empty:
            filterd_sam["MS"] = filterd_sam["MS"].str.extract(r"ms:i:(\d+)").fillna(0).astype(int)
            filterd_sam["MAPQ"] = filterd_sam["MAPQ"].astype(int)
            filterd_sam["POS"] = filterd_sam["POS"].astype(int)
            if second != "-":
                for rname, group in filterd_sam.groupby("RNAME"):
                    first_subset = group[group["QNAME"].str.contains(first, na=False)]
                    first_subset = filter_best_ms(first_subset)

                    second_subset = group[group["QNAME"].str.contains(second, na=False)]
                    second_subset = filter_best_ms(second_subset)

                    if not first_subset.empty and not second_subset.empty:
                        first_start = int(first_subset["POS"].min())
                        first_end = int(first_subset["POS"].max())
                        second_start = int(second_subset["POS"].min())
                        second_end = int(second_subset["POS"].max())
                        start = min(first_start, second_start)
                        end = max(first_end, second_end)
                        extraction_regions.append((rname, start, end, first, second))
                        return extraction_regions

                first_subset = filterd_sam[filterd_sam["QNAME"].str.contains(first, na=False)]
                first_subset = filter_best_ms(first_subset)
                if not first_subset.empty:
                    first_contig = first_subset["RNAME"].astype(str).iloc[0]
                    start = first_subset["POS"].iloc[0]
                    end = get_length_contig(sam, first_contig)
                    extraction_regions.append((first_contig, start, end, first, "-"))


                second_subset = filterd_sam[filterd_sam["QNAME"].str.contains(second, na=False)]
                second_subset = filter_best_ms(second_subset)
                if not second_subset.empty:
                    second_contig = second_subset["RNAME"].astype(str).iloc[0]
                    end = int(second_subset["POS"].astype(int).iloc[0])
                    extraction_regions.append((second_contig, 1, end, "-", second))

            else:
                first_subset = filterd_sam[filterd_sam["QNAME"].str.contains(first, na=False)]
                first_subset = filter_best_ms(first_subset)
                firts_contig = first_subset["RNAME"].astype(str).iloc[0]
                start = int(first_subset["POS"].astype(int).iloc[0])
                end = get_length_contig(sam, firts_contig)
                reverse_strand = (first_subset["FLAG"].astype(int).iloc[0] & 16) != 0
                if reverse_strand:
                    extraction_regions.append((firts_contig, 1, start, first, "-"))
                else:
                    extraction_regions.append((firts_contig, start, end, first, "-"))
        return extraction_regions
    except:
        console_log.error(f"Error processing SAM file: {sam.stem}")
        return []


def extract(cwd: Union[str, Path], assembly_fasta: Union[str, Path], directory : Union[str, Path], first: str, second: str, sample: str, immuno_region: str, verbose: bool):
    """
    Extracts a sequence from an assembly FASTA file based on flanking genes,
    and writes it to an output FASTA file.
    """
    if second == "":
        second = "-"

    all_log_data = []

    sam = cwd / "tmp" / "mapped_genes" / assembly_fasta.with_suffix(".sam").name
    extraction_regions = get_positions_and_name(sam, first, second)
    for contig, start, end, flanking_gene_one, flanking_gene_second in extraction_regions:
        Extraction_status = "extracted" if first == flanking_gene_one and second == flanking_gene_second else "fragmented"

        output_file = directory / f"{sample}__{flanking_gene_one}__{flanking_gene_second}__{contig}__{Extraction_status}__{immuno_region}.fasta"
        if not output_file.is_file():
            cmd = f"samtools faidx {assembly_fasta} {contig}:{start}-{end}"
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            concatenated_sequence = "".join(result.stdout.strip().splitlines()[1:])

            with open(output_file, 'w') as file:
                file.write(f">{output_file.stem}\n{str(Seq(concatenated_sequence))}")

            n_contigs = len(list(re.finditer(r'N{100,}', str(Seq(concatenated_sequence)).upper()))) + 1
            log_data = {
                "Region": immuno_region,
                "Contig": contig,
                "5'-flanking gene": flanking_gene_one,
                "3'-flanking gene": flanking_gene_second,
                "5 coords": int(start),
                "3 coords": int(end),
                "Extraction status": Extraction_status.capitalize(),
                "Assembly": str(assembly_fasta),
                "Assembly file": str(assembly_fasta.name),
                "Output": str(output_file),
                "Output file": str(output_file.name),
                "Contig counts": n_contigs,
                "Amount basepairs": len(concatenated_sequence),
            }
            all_log_data.append(log_data)

    return all_log_data



@log_error()
def region_main(flanking_genes: dict[list[str]], assembly_dir: Union[str, Path], threads: int, verbose: bool):
    """
    Main function that processes SAM files to create region-specific assembly files
    based on flanking genes specified in the configuration.
    """
    cwd = Path.cwd()
    directory = cwd / "tmp/region"
    make_dir(directory)

    assembly_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in Path(assembly_dir).glob(ext)]

    log_file = cwd / "broken_regions.json"
    if log_file.is_file():
        with open(log_file) as f:
            output_json = json.load(f)
    else:
        output_json = {}

    tasks = []

    for region, extract_flanking_genes in flanking_genes.items():
        for assembly in assembly_files:
            tasks.append((cwd, assembly, directory, extract_flanking_genes[0], extract_flanking_genes[1], assembly.stem, region, verbose))
    max_jobs = calculate_available_resources(max_cores=threads, threads=4, memory_per_process=12)
    total_tasks = len(tasks)
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {executor.submit(extract, *task): task for task in tasks}
        with tqdm(total=total_tasks, desc='Extracting regions', unit='task') as pbar:
            for future in as_completed(futures):
                log_data = future.result()
                if len(log_data) > 0:
                    assembly_name = futures[future][1].name
                    for entry in log_data:
                        immuno_region = entry["Region"]
                        if assembly_name not in output_json:
                            output_json[assembly_name] = {}
                        if immuno_region not in output_json[assembly_name]:
                            output_json[assembly_name][immuno_region] = []
                        if entry not in output_json[assembly_name][immuno_region]:
                            output_json[assembly_name][immuno_region].append(entry)
                pbar.update(1)

                with open(log_file, 'w') as f:
                    json.dump(output_json, f, indent=4)

    with open(log_file, 'w') as f:
        json.dump(output_json, f, indent=4)

    if any(directory.iterdir()):
        file_log.info("Region extraction completed successfully")
    else:
        file_log.error("No regions were extracted")
        raise Exception("No regions extracted.")
