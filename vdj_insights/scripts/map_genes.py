"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

import time
from pathlib import Path
import shutil
import subprocess
from typing import Union
from tqdm import tqdm
from time import sleep
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed

from .util import make_dir, unzip_file, calculate_available_resources
from .property import log_error

from .logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def download_flanking_genes(gene: str, path: Path, species="Homo sapiens") -> None:
    """
    Downloads flanking gene sequences for the specified gene and species using the NCBI datasets CLI.
    Extracts the sequences to a specified directory and processes them into a FASTA file.

    Args:
        gene (str): The gene symbol for which to download flanking sequences.
        path (Path): Directory where the downloaded files will be saved.
        species (str, optional): Species name for the gene. Defaults to "Homo sapiens".

    Raises:
        subprocess.CalledProcessError: If the command to download the genes fails.
        Exception: If an unexpected error occurs during file processing, logs the error and raises an exception.
    """

    path = Path(path)
    output_zip = path / f"{gene}.zip"
    output_fna = path / f"{gene}.fna"
    if not output_fna.is_file():
        command = f'datasets download gene symbol {gene} --taxon "{species}" --include gene --filename {output_zip}'

        max_retries = 10
        retry_delay = 5
        for attempt in range(1, max_retries + 1):
            try:
                subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                break
            except subprocess.CalledProcessError as e:
                print(e.stderr)
                if attempt < max_retries:
                    time.sleep(retry_delay)

        unzip_file(output_zip, path)

        ncbi_fna_path = path / "ncbi_dataset" / "data" / "gene.fna"
        ncbi_fna_path.rename(output_fna)
        output_zip.unlink()

        readme_path = path / "README.md"
        if readme_path.exists():
            readme_path.unlink()

        shutil.rmtree(path / "ncbi_dataset")
        sleep(2)
        file_log.info(f"Downloaded and processed flanking genes for {gene}")


@log_error()
def combine_genes(path: Union[str, Path], flanking_output: Union[str, Path]) -> Path:
    """
    Combines all gene sequences in a specified directory into a single FASTA file.

    Args:
        path (Path): Directory containing the individual gene FASTA files to combine.
        flanking_output (Path): Path to the output FASTA file where all gene sequences will be combined.

    Returns:
        Path: The path to the combined FASTA file.

    Raises:
        Exception: If the genes cannot be combined, logs the error and raises an exception.
    """
    path = Path(path)
    flanking_output = Path(flanking_output)

    if not flanking_output.is_file():
        with open(flanking_output, 'w') as outfile:
            gene_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in path.glob(ext)]
            for gene_file in gene_files:
                if gene_file.is_file() and gene_file != flanking_output:
                    with open(gene_file, 'r') as fasta:
                        for record in SeqIO.parse(fasta, 'fasta'):
                            outfile.write(f">{record.description.replace(' ', '_')}\n{record.seq}\n")
                    file_log.info(f'Added {gene_file} to {flanking_output}')
        file_log.info(f'All files in {path} have been combined into {flanking_output}')
    return flanking_output


def map_flanking_genes(output_dir: Path, flanking_genes: Path, assembly_file: Path, threads : int = 8) -> None:
    """
    Maps flanking genes to an assembly file using Minimap2, producing a SAM file.

    Args:
        output_dir (Path): Current working directory.
        flanking_genes (Path): Path to the FASTA file containing flanking genes.
        assembly_file (Path): Path to the assembly FASTA file to which flanking genes will be mapped.
        threads (int, optional): Number of threads to use for Minimap2. Defaults to 8.

    Raises:
        subprocess.CalledProcessError: If the Minimap2 command fails.
        Exception: If an unexpected error occurs during the mapping process, logs the error and raises an exception.
    """
    sam_file = output_dir / assembly_file.with_suffix(".sam").name
    if not sam_file.is_file():
        command = (f'minimap2 -ax asm5 -t {threads} {assembly_file} {flanking_genes} > {sam_file}')
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        file_log.info(f"Mapped flanking genes from {flanking_genes} to {assembly_file}")


def map_main(flanking_genes: dict[list[str]], assembly_dir: Union[str, Path], species: str, threads: int) -> None:
    """
    Main function that coordinates the downloading, combining, and mapping of flanking genes to assemblies.

    Args:
        flanking_genes (list): List of gene symbols representing the flanking genes to be processed.
        assembly_dir (str or Path): Directory containing assembly FASTA files.
        species (str): Species name for which the flanking genes are to be downloaded and processed.

    Raises:
        Exception: If the main mapping process fails, logs the error and raises an exception.
    """
    cwd = Path.cwd()
    flanking_genes_dir = cwd / "flanking_genes"
    map_flanking_genes_dir = cwd / "mapped_genes"
    make_dir(flanking_genes_dir)
    make_dir(map_flanking_genes_dir)

    console_log.info(f"Downloading flanking genes for {species}")

    flanking_genes_list = [gene for genes in flanking_genes.values() for gene in genes]
    for gene in flanking_genes_list:
        if gene != "-" and len(gene) > 1:
            download_flanking_genes(gene, flanking_genes_dir, species)

    gene_output = combine_genes(flanking_genes_dir, flanking_genes_dir / "all_genes.fna")
    assembly_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in (cwd / assembly_dir).glob(ext)]
    console_log.info(f"Number of assembly files: {len(assembly_files)}")

    tasks = [
        (map_flanking_genes_dir, gene_output, Path(assembly_file), 4)
        for assembly_file in assembly_files
    ]
    total_tasks = len(tasks)

    max_jobs = calculate_available_resources(max_cores=threads, threads=4, memory_per_process=12)
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {executor.submit(map_flanking_genes, *arg): arg for arg in tasks}
        with tqdm(total=total_tasks, desc="Mapping flanking genes", unit='task') as pbar:
            for future in as_completed(futures):
                future.result()
                pbar.update(1)

    file_log.info("Extract main process completed successfully")
