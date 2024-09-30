from pathlib import Path
import shutil
import subprocess
from time import sleep
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed

from util import make_dir, unzip_file, calculate_available_resources
from property import log_error

from logger import console_logger, file_logger

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
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        unzip_file(output_zip, path)

        ncbi_fna_path = path / "ncbi_dataset" / "data" / "gene.fna"
        ncbi_fna_path.rename(output_fna)
        output_zip.unlink()

        readme_path = path / "README.md"
        if readme_path.exists():
            readme_path.unlink()

        shutil.rmtree(path / "ncbi_dataset")
        sleep(2)
        console_log.info(f"Downloaded and processed flanking genes for {gene}")


@log_error()
def combine_genes(path: str | Path, flanking_output: str | Path) -> Path:
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
            extensions = ["*.fna", "*.fasta", "*.fa"]
            gene_files = [file for ext in extensions for file in path.glob(ext)]

            for gene_file in gene_files:
                if gene_file.is_file() and gene_file != flanking_output:
                    with open(gene_file, 'r') as fasta:
                        for record in SeqIO.parse(fasta, 'fasta'):
                            outfile.write(f">{record.description.replace(' ', '_')}\n{record.seq}\n")
                    console_log.info(f'Added {gene_file} to {flanking_output}')
        console_log.info(f'All files in {path} have been combined into {flanking_output}')
    return flanking_output


def map_flanking_genes(output_dir: Path, flanking_genes: Path, assembly_file: Path, threads=8) -> None:
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
    awk_command = "awk '$1 ~ /^@/ || $5 == 60 {print $1, $2, $3, $4, $5}' | awk '!seen[$3,$4]++'"
    if not sam_file.is_file():
        command = (
            f'minimap2 -ax asm5 --secondary=no -t {threads} '
            f'{assembly_file} {flanking_genes} | {awk_command} > {sam_file}'
        )
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        console_log.info(f"Mapped flanking genes from {flanking_genes} to {assembly_file}")


def map_main(flanking_genes: list[str], assembly_dir: str | Path, species: str, memory_per_process=8, buffer_percentage=20) -> None:
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
    make_dir(flanking_genes_dir)
    for gene in flanking_genes:
        if gene != "-":
            download_flanking_genes(gene, flanking_genes_dir, species)

    gene_output = combine_genes(flanking_genes_dir, flanking_genes_dir / "all_genes.fna")

    assembly_dir = cwd / assembly_dir
    map_flanking_genes_dir = cwd / "mapped_genes"
    make_dir(map_flanking_genes_dir)

    extensions = ["*.fna", "*.fasta", "*.fa"]
    assembly_files = [file for ext in extensions for file in assembly_dir.glob(ext)]

    max_jobs = calculate_available_resources(max_cores=24, threads=4, memory_per_process=12)

    args = [
        (map_flanking_genes_dir, gene_output, Path(assembly_file), 4)
        for assembly_file in assembly_files
    ]
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {executor.submit(map_flanking_genes, *arg): arg for arg in args}
        for future in as_completed(futures):
            future.result()

    console_log.info("Extract main process completed successfully")


if __name__ == "__main__":
    map_main(["SALL2", "DAD1", "MGAM2", "EPHB6", "EPDR1", "VPS41"], path, "Homo sapiens")
