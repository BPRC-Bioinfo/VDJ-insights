from pathlib import Path
import shutil
import subprocess
from time import sleep
import zipfile
from Bio import SeqIO
from logger import custom_logger
"""
Used python packages:
    1. biopython
    
Used CLI packages:
    1. ncbi datasets
    2. minimap2
"""
# Method for logging the current states of the program.
logger = custom_logger(__name__)


def make_dir(dir):
    """
    Create a directory if not existing.

    Args:
        location (str): Path of the directory to create.
    """
    try:
        Path(dir).mkdir(parents=True, exist_ok=True)
        logger.debug(f"Directory created or already exists: {dir}")
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")


def unzip_file(file_path, dir):
    extract_to_path = Path(dir)
    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to_path)
        logger.info(f'Extracted {file_path} to {extract_to_path}')
    except Exception as e:
        logger.error(f"Failed to extract {file_path} to {extract_to_path}: {e}")


def download_flanking_genes(gene, dir: Path, species="homo sapiens"):
    """
    Downloads and extracts flanking gene sequences for the specified gene and species.

    :param gene: Gene symbol to download.
    :param dir: Directory where the files will be saved.
    :param species: Species name, default is "homo sapiens".
    """
    dir = Path(dir)
    output_zip = dir / f"{gene}.zip"
    output_fna = dir / f"{gene}.fna"
    if not output_fna.is_file():
        try:
            command = f'datasets download gene symbol {gene} --taxon "{species.capitalize()}" --include gene --filename {output_zip}'
            result = subprocess.run(
                command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            unzip_file(output_zip, dir)
            ncbi_fna_path = dir / "ncbi_dataset" / "data" / "gene.fna"
            ncbi_fna_path.rename(output_fna)
            output_zip.unlink()
            readme_path = dir / "README.md"
            if readme_path.exists():
                readme_path.unlink()
            shutil.rmtree(dir / "ncbi_dataset")
            sleep(2)
            logger.info(f"Downloaded and processed flanking genes for {gene}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error occurred while running command: {e}")
            logger.error(e.stdout.decode())
            logger.error(e.stderr.decode())
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")


def combine_genes(dir, flanking_output):
    dir = Path(dir)
    flanking_output = Path(flanking_output)

    if not dir.is_dir():
        make_dir(dir)

    if not flanking_output.is_file():
        try:
            with open(flanking_output, 'w') as outfile:
                # Using a single list comprehension to gather files with desired extensions
                extensions = ["*.fna", "*.fasta", "*.fa"]
                gene_files = [
                    file for ext in extensions for file in dir.glob(ext)]

                for gene_file in gene_files:
                    if gene_file.is_file() and gene_file != flanking_output:
                        with open(gene_file, 'r') as fasta:
                            for record in SeqIO.parse(fasta, 'fasta'):
                                outfile.write(
                                    f">{record.description.replace(' ', '_')}\n{record.seq}\n")
                        logger.info(f'Added {gene_file} to {flanking_output}')
            logger.info(
                f'All files in {dir} have been combined into {flanking_output}')
        except Exception as e:
            logger.error(f"Failed to combine genes: {e}")
    return flanking_output


def map_flanking_genes(cwd, flanking_genes: Path, assembly_file: Path, threads=8):
    output_dir = cwd / f"mapped_genes"
    if not output_dir.is_dir():
        make_dir(output_dir)
    sam_file = output_dir / assembly_file.with_suffix(".sam").name
    awk_command = "awk '$1 ~ /^@/ || $5 == 60 {print $1, $2, $3, $4, $5}' | awk '!seen[$3,$4]++'"
    if not sam_file.is_file():
        try:
            command = (
                f'minimap2 -ax asm5 --secondary=no -t {threads} '
                f'{assembly_file} {flanking_genes} | {awk_command} > {sam_file}'
            )
            subprocess.run(
                command, shell=True, check=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            logger.info(
                f"Mapped flanking genes from {flanking_genes} to {assembly_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error occurred while running command: {e}")
            logger.error(e.stdout.decode())
            logger.error(e.stderr.decode())
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")


def extract_main(flanking_genes, assembly_dir, species):
    try:
        cwd = Path.cwd()
        flanking_genes_dir = cwd / "flanking_genes"
        make_dir(flanking_genes_dir)
        for gene in flanking_genes:
            download_flanking_genes(gene, flanking_genes_dir, species)
        gene_output = combine_genes(
            flanking_genes_dir, flanking_genes_dir / "all_genes.fna")
        assembly_dir = cwd / assembly_dir
        extensions = ["*.fna", "*.fasta", "*.fa"]
        assembly_files = [
            file for ext in extensions for file in assembly_dir.glob(ext)]
        for assembly_file in assembly_files:
            map_flanking_genes(cwd, gene_output, Path(assembly_file))
        logger.info("Extract main process completed successfully")
    except Exception as e:
        logger.error(f"Failed in extract_main: {e}")


if __name__ == "__main__":
    extract_main()
