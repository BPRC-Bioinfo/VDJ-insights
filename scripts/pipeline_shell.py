import argparse
import re
import shutil
import sys
import zipfile
from time import sleep
from Bio import SeqIO
from pathlib import Path
import subprocess
from logger import custom_logger
import yaml

# Method for logging the current states of the program.
logger = custom_logger(__name__)

CONFIG = {}


def load_config(config_file):
    if config_file.exists():
        with open(config_file, 'r') as file:
            return yaml.safe_load(file)
    else:
        logger.error("No configuration file provided, closing application!")
        sys.exit()


def validate_read_files(file_path):
    data_path = Path(file_path)
    validate_files(file_path)
    if data_path.suffixes != ['.fastq', '.gz']:
        logger.error(
            f"Invalid file type for {file_path}. Must be a '.fastq.gz' file.")
        raise argparse.ArgumentTypeError(
            f"The file {file_path} must be a '.fastq.gz' file.")
    logger.info(f"Validated read file: {file_path}")
    return data_path


def validate_files(file_path):
    data_path = Path(file_path)
    if not data_path.is_file():
        logger.error(f"File does not exist or is a directory: {file_path}")
        raise argparse.ArgumentTypeError(
            f"The file {file_path} does not exist or is a directory.")
    logger.info(f"Validated file: {file_path}")
    return data_path


def validate_reference(value):
    fasta_extensions = {".fasta", ".fna"}
    if Path(value).suffix in fasta_extensions:
        validate_files(value)
        logger.info(f"Validated reference genome file: {value}")
        return value

    accession_patterns = [
        r"^[A-Z]{1,2}\d{5,6}$",
        r"^[A-Z]{2}_\d{6,}\.\d+$",
        r"^[A-Z]{3}_\d{9}\.\d+$",
        r"^[A-Z]{2}\d{6,}\.\d+$",
    ]
    if any(re.match(pattern, value) for pattern in accession_patterns):
        logger.info(f"Validated reference genome accession: {value}")
        return value

    logger.error(
        f"Invalid reference genome input: '{value}'. Must be a .fasta or .fna file, or a valid accession code.")
    raise argparse.ArgumentTypeError(
        f"Invalid reference genome input: '{value}'. Must be a .fasta or .fna file, or a valid accession code.")


def validate_flanking_genes(value):
    flanking_genes = [gene.strip().upper() if gene.strip() !=
                      '-' else '' for gene in value.split(',')]
    if len(flanking_genes) % 2 == 1:
        logger.error(
            f"The specified flanking genes: {flanking_genes} should be even numbers.")
        raise argparse.ArgumentTypeError(
            f"The specified flanking genes: {flanking_genes} should be even numbers (e.g., 2, 4, 6, 8) rather than odd (e.g., 1, 3, 5).")
    logger.info(f"Validated flanking genes: {flanking_genes}")
    return flanking_genes


def validate_chromosome(value):
    chromosomes = [chromosome.strip() for chromosome in value.split(',')]
    valid_chromosomes = set(map(str, range(1, 23))) | {"X", "Y"}

    if all(chromosome in valid_chromosomes for chromosome in chromosomes):
        logger.info(f"Validated chromosomes: {chromosomes}")
        return chromosomes
    else:
        logger.error(
            f"Invalid chromosome list: '{value}'. All values must be integers between 1-22, or 'X', 'Y'.")
        raise argparse.ArgumentTypeError(
            f"Invalid chromosome list: '{value}'. All values must be integers between 1-22, or 'X', 'Y'.")


def argparser_setup():
    parser = argparse.ArgumentParser(description="Process sequencing data.")

    reads_group = parser.add_argument_group(
        'reads', 'Arguments related to read files')
    reads_group.add_argument('-ont', '--nanopore', required=True, type=validate_read_files,
                             help='Path to the Oxford Nanopore Technologies reads file.')
    reads_group.add_argument('-pb', '--pacbio', required=True, type=validate_read_files,
                             help='Path to the Pacific Biosciences reads file.')

    reference_group = parser.add_argument_group(
        'reference', 'Arguments related to the reference genome')
    reference_group.add_argument('-ref', '--reference', type=validate_reference, required=False,
                                 help='Path to the reference genome file if one is already present (optional).')

    analysis_group = parser.add_argument_group(
        'analysis', 'Arguments related to analysis settings')
    analysis_group.add_argument('-s', '--species', type=str.capitalize,
                                required=True, help='Species name, e.g., Homo sapiens.')
    analysis_group.add_argument('-f', '--flanking-genes', type=validate_flanking_genes,
                                help='Comma-separated list of flanking genes, e.g., MGAM2,EPHB6. Add them as pairs.')
    analysis_group.add_argument('-r', '--receptor-type', required=True, type=str.upper, choices=[
                                'TR', 'IG'], help='Type of receptor to analyze: TR (T-cell receptor) or IG (Immunoglobulin).')
    analysis_group.add_argument('-c', '--chromosomes', type=validate_chromosome,
                                required='--default' not in sys.argv, help='List of chromosomes where TR or IG is located.')
    analysis_group.add_argument('-t', '--threads', type=int,
                                required=False, default=8, help='Amount of threads to run the analysis.')
    analysis_group.add_argument('--default', action='store_true',
                                help='Use default settings. Cannot be used with -f/--flanking-genes or -c/--chromosomes.')

    args = parser.parse_args()

    if args.default and (args.flanking_genes or args.chromosomes):
        parser.error(
            "--default cannot be used with -f/--flanking-genes or -c/--chromosomes")

    return args


def make_dir(dir):
    Path(dir).mkdir(parents=True, exist_ok=True)
    logger.info(f"Created directory: {dir}")


def key_sort(s):
    """
    Sort strings with numbers in a natural way (e.g., 1, 2, 10 instead of 1, 10, 2).
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split('(\d+)', s)]


def split_chromosomes(cwd, fasta_path):
    global CONFIG
    output_dir = cwd / 'reference' / 'chromosomes'
    genome_dir = cwd / 'reference' / 'genome'
    new_genome_file = genome_dir / 'new_reference.fasta'
    if not new_genome_file.is_file():
        [make_dir(i) for i in [output_dir, genome_dir]]
        write_chromosomes(fasta_path, output_dir, new_genome_file)
    else:
        for chromosome in sorted(output_dir.glob("reference_chr*.fasta"), key=lambda x: key_sort(x.stem)):
            chromosome_number = chromosome.stem.replace('reference_chr', '')
            CONFIG.setdefault('ALL_CHROMOSOMES', []).append(chromosome_number)


def write_chromosomes(fasta_path, output_dir, new_genome_file):
    files = {}
    files["all"] = new_genome_file.open('a')
    for record in SeqIO.parse(fasta_path, 'fasta'):
        process_record(record, output_dir, files)
    close_files(files)
    logger.info(f"Files have been written to: {output_dir}")


def process_record(record, output_dir, files):
    global CONFIG
    splitted = record.description.split(' ')
    try:
        chromosome, chromosome_number = extract_chromosome_info(splitted)
        if chromosome not in files:
            CONFIG.setdefault('ALL_CHROMOSOMES', []).append(
                chromosome_number.rstrip(","))
            file_path = output_dir / f"{chromosome}.fasta"
            file_path_txt = output_dir / f"{chromosome}.txt"
            files[chromosome] = file_path.open('a')
            files[f"{chromosome}_txt"] = file_path_txt.open('a')
        SeqIO.write(record, files[chromosome], 'fasta')
        files[f"{chromosome}_txt"].write(f"{record.id}\n")
        SeqIO.write(record, files['all'], 'fasta')
    except (StopIteration, IndexError) as e:
        logger.warning(
            f"Error processing record: {record.description}. Error: {e}")


def extract_chromosome_info(splitted):
    chromosome, chromosome_number = next(
        (f"reference_chr{splitted[i+1].rstrip(',')}", splitted[i+1])
        for i, x in enumerate(splitted)
        if x == "chromosome"
    )
    return chromosome, chromosome_number


def close_files(files):
    for file in files.values():
        file.close()


def rename_files(cwd: Path, file: Path, read_type: str):
    moved_dir = cwd / 'downloads'
    make_dir(moved_dir)
    new_file = get_new_file_path(file, read_type, moved_dir)
    logger.info(f"Renamed file {file} to {new_file}")
    return file.stem.split('_')[0].upper(), new_file


def get_new_file_path(file, read_type, moved_dir):
    conversion = {"ONT": "nanopore", "NANOPORE": "nanopore",
                  "PACBIO": "pacbio", "PB": "pacbio"}
    stripped_extensions = file.with_suffix('').with_suffix('')
    sample = stripped_extensions.stem.split('_')[0].upper()
    new_sample_name = sample if sample not in conversion else "SAMPLE"
    return moved_dir / f"{new_sample_name}_{read_type}.fastq.gz"


def filter_and_move_files(cwd: Path, nanopore_file: Path, pacbio_file: Path):
    ont_sample, ont_file = rename_files(cwd, nanopore_file, 'nanopore')
    pb_sample, pb_file = rename_files(cwd, pacbio_file, 'pacbio')
    move_files(nanopore_file, pacbio_file, ont_file,
               pb_file, ont_sample, pb_sample)


def move_files(nanopore_file, pacbio_file, ont_file, pb_file, ont_sample, pb_sample):
    if ont_sample != pb_sample:
        ont_file = ont_file.with_name(f"SAMPLE_ont.fastq.gz")
        pb_file = pb_file.with_name(f"SAMPLE_pacbio.fastq.gz")
    shutil.move(nanopore_file, ont_file)
    shutil.move(pacbio_file, pb_file)
    logger.info(
        f"Moved {nanopore_file} to {ont_file} and {pacbio_file} to {pb_file}")


def unzip_file(file_path, dir):
    extract_to_path = Path(dir)
    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to_path)
        logger.info(f'Extracted {file_path} to {extract_to_path}')
    except Exception as e:
        logger.error(f"Failed to extract {file_path} to {extract_to_path}: {e}")


def download_flanking_genes(gene, dir: Path, species):
    dir = Path(dir)
    output_zip = dir / f"{gene}.zip"
    output_fna = dir / f"{gene}.fna"
    if not output_fna.is_file():
        run_download_command(gene, species, output_zip, output_fna, dir)


def run_download_command(gene, species, output_zip, output_fna, dir):
    try:
        command = f'datasets download gene symbol {gene} --taxon "{species.capitalize()}" --include gene --filename {output_zip}'
        result = subprocess.run(
            command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process_downloaded_files(output_zip, dir, output_fna)
        logger.info(f"Downloaded and processed flanking genes for {gene}")
    except subprocess.CalledProcessError as e:
        log_subprocess_error(e)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")


def process_downloaded_files(output_zip, dir, output_fna):
    unzip_file(output_zip, dir)
    ncbi_fna_path = dir / "ncbi_dataset" / "data" / "gene.fna"
    ncbi_fna_path.rename(output_fna)
    cleanup_downloaded_files(output_zip, dir)


def cleanup_downloaded_files(output_zip, dir):
    output_zip.unlink()
    readme_path = dir / "README.md"
    if readme_path.exists():
        readme_path.unlink()
    shutil.rmtree(dir / "ncbi_dataset")
    sleep(2)


def log_subprocess_error(e):
    logger.error(f"Error occurred while running command: {e}")
    logger.error(e.stdout.decode())
    logger.error(e.stderr.decode())


def download_reference_genome(genome_code, reference_dir: Path):
    reference_dir = Path(reference_dir)
    make_dir(reference_dir)
    output_zip = reference_dir / f"{genome_code}.zip"
    output_fna = reference_dir / f"reference.fna"
    if not output_fna.is_file():
        run_reference_download_command(
            genome_code, output_zip, reference_dir, output_fna)
    return output_fna


def run_reference_download_command(genome_code, output_zip, reference_dir, output_fna):
    try:
        command = f'datasets download genome accession {genome_code} --include genome --filename {output_zip}'
        logger.info(f"Running command: {command}")
        result = subprocess.run(
            command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process_reference_files(output_zip, reference_dir,
                                output_fna, genome_code)
    except subprocess.CalledProcessError as e:
        log_subprocess_error(e)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")


def process_reference_files(output_zip, reference_dir, output_fna, genome_code):
    unzip_file(output_zip, reference_dir)
    ncbi_fna_path = (reference_dir / "ncbi_dataset" /
                     "data" / genome_code).glob('*.fna')
    ncbi_fna_path = list(ncbi_fna_path)[0]
    ncbi_fna_path.rename(output_fna)
    cleanup_downloaded_files(output_zip, reference_dir)
    logger.info(f"Downloaded the reference genome: {genome_code}")


def deep_merge(d1, d2):
    for key, value in d2.items():
        if key in d1:
            if isinstance(d1[key], dict) and isinstance(value, dict):
                deep_merge(d1[key], value)
            else:
                d1[key] = value
        else:
            d1[key] = value


def loop_flanking_genes(cwd, args):
    receptor_chromosomes = set()
    flanking_genes_dir = prepare_flanking_genes_directory(cwd)
    species_dict = get_species_dict(cwd, args)
    flanking_genes = species_dict.get(
        args.receptor_type, {}).get('FLANKING_GENES', [])
    args.flanking_genes = flanking_genes
    for gene in args.flanking_genes:
        process_flanking_gene(gene, flanking_genes_dir,
                              args.species, receptor_chromosomes)
    args.chromosomes = list(receptor_chromosomes)


def prepare_flanking_genes_directory(cwd):
    default_setting_file = cwd / '_config' / 'species.yaml'
    flanking_genes_dir = cwd / "flanking_genes"
    make_dir(flanking_genes_dir)
    return flanking_genes_dir


def get_species_dict(cwd, args):
    default_setting_file = cwd / '_config' / 'species.yaml'
    default_dict = load_config(default_setting_file)
    species_key = args.species.replace(' ', '_')
    return default_dict.get(species_key, default_dict.get('default', {}))


def process_flanking_gene(gene, flanking_genes_dir, species, receptor_chromosomes):
    download_flanking_genes(gene, flanking_genes_dir, species)
    gene_file = flanking_genes_dir / f"{gene}.fna"
    records = SeqIO.to_dict(SeqIO.parse(gene_file, 'fasta'))
    first_key = next(iter(records))
    first_record = records[first_key]
    chromosome_info = first_record.description.split('=')[-1].strip(']')
    number, trailing = extract_chromosome_number_and_trailing(chromosome_info)
    receptor_chromosomes.add(number + trailing)


def extract_chromosome_number_and_trailing(chromosome_info):
    number = ''.join(filter(str.isdigit, chromosome_info))
    trailing = ''.join(filter(str.isalpha, chromosome_info))
    return number, trailing


def create_config(cwd: Path, args):
    global CONFIG
    initialize_config(args)
    rss_config = load_config(cwd / '_config' / 'rss.yaml')
    deep_merge(CONFIG, rss_config.get(args.receptor_type, {}))
    config_file = cwd / 'config' / 'config.yaml'
    make_dir(config_file.parent)
    save_config_to_file(config_file)


def initialize_config(args):
    global CONFIG
    CONFIG.setdefault('SPECIES', {
        'name': args.species,
        'genome': args.reference,
        'cell': args.receptor_type
    })
    CONFIG.setdefault("ASSEMBLY_CHROMOSOMES", sorted(
        args.chromosomes, key=lambda x: key_sort(x)))
    CONFIG.setdefault("FLANKING_GENES", args.flanking_genes)
    CONFIG.setdefault("BUSCO_DATASET", "primates_odb10")
    CONFIG.setdefault("HAPLOTYPES", [1, 2])


def save_config_to_file(file_name):
    global CONFIG
    with open(file_name, 'w') as f:
        yaml.dump(CONFIG, f, default_flow_style=False)


def run_snakemake(args, snakefile="Snakefile"):
    try:
        result = subprocess.run(
            f"snakemake -s {snakefile} --cores {args.threads} --use-conda -prn",
            shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )
        logger.info("Snakemake check ran successfully.")
        # logger.info("Running the complete pipeline.")
        # subprocess.run(
        #     f"snakemake -s {snakefile} --cores {args.threads} --use-conda -pr",
        #     shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        # )
    except subprocess.CalledProcessError as e:
        logger.error(f"Snakemake failed with return code {e.returncode}.")
        logger.error(f"Snakemake output: {e.stdout.decode()}")
        logger.error(f"Application is shutting down.")
        sys.exit(e.returncode)


def main():
    cwd = Path.cwd()
    args = argparser_setup()
    logger.info("Starting main process")
    filter_and_move_files(cwd, args.nanopore, args.pacbio)
    if args.reference and Path(args.reference).suffix in {'.fasta', '.fna'}:
        fasta_path = cwd / args.reference
        split_chromosomes(cwd, fasta_path)
    else:
        genome_dir = cwd / 'reference' / 'genome'
        genome = download_reference_genome(args.reference, genome_dir)
        split_chromosomes(cwd, genome)
    if args.default:
        cwd = Path.cwd()
        loop_flanking_genes(cwd, args)
    create_config(cwd, args)
    run_snakemake(args)
    logger.info("Main process completed")


if __name__ == '__main__':
    main()
