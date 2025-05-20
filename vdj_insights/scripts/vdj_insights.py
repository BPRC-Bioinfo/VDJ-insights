"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

import argparse
import os
import re
import shutil
import sys
import threading
import webbrowser
from time import sleep
from pathlib import Path
import subprocess
import json

import yaml
import pandas as pd
from Bio import SeqIO

from .logger import console_logger, file_logger
from .annotation import main as annotation_main
from .IMGT_scrape import main as imgt_main
from .env_manager import create_and_activate_env, deactivate_env
from .util import make_dir, validate_file, validate_directory, validate_input, load_config, unzip_file


console_log = console_logger(__name__)
file_log = file_logger(__name__)


def cwd_setup(output_dir):
    settings_dir = Path(__file__).resolve().parent.parent
    output_dir = Path(output_dir).resolve()
    make_dir(output_dir)
    copy_flask(output_dir / 'flask', settings_dir)
    os.chdir(str(output_dir))
    return settings_dir, output_dir


def copy_flask(output_dir, reset=False):
    settings_dir = Path(__file__).resolve().parent.parent
    if not output_dir.is_dir() or reset:
        shutil.copytree(str(settings_dir / "flask"), str(output_dir), dirs_exist_ok=True)


def validate_read_files(file_path):
    """
    Validates a file path to ensure it points to a valid .fastq.gz file.
    This function checks that the file exists and has the correct extension.

    Args:
        file_path (str): Path to the read file.

    Returns:
        Path: Validated Path object pointing to the read file.

    Raises:
        argparse.ArgumentTypeError: If the file is not a .fastq.gz file or does not exist.
    """
    data_path = Path(file_path).resolve()
    validate_files(file_path)
    if data_path.suffixes != ['.fastq', '.gz']:
        file_log.error(f"Invalid file type for {file_path}. Must be a '.fastq.gz' file.")
        raise argparse.ArgumentTypeError(
            f"The file {file_path} must be a '.fastq.gz' file.")
    file_log.info(f"Validated read file: {file_path}")
    return data_path


def validate_files(file_path):
    """
    Validates that the given file path points to an existing file.
    If the path is invalid or points to a directory, an error is logged
    and an exception is raised.

    Args:
        file_path (str): Path to the file that needs validation.

    Returns:
        Path: Validated Path object pointing to the file.

    Raises:
        argparse.ArgumentTypeError: If the path is invalid or points to a directory.
    """
    data_path = Path(file_path)
    if not data_path.is_file():
        file_log.error(f"File does not exist or is a directory: {file_path}")
        raise argparse.ArgumentTypeError(
            f"The file {file_path} does not exist or is a directory.")
    file_log.info(f"Validated file: {file_path}")
    return data_path


def validate_reference(value):
    """
    Validates a reference genome input, ensuring it is either a valid
    .fasta/.fna file or a valid accession code. The function checks the
    file extension against known FASTA extensions and validates the
    presence of the file. If the input is an accession code, it is checked
    against known accession patterns.

    Args:
        value (str): The input reference genome file path or accession code.

    Returns:
        str: The validated reference genome path or accession code.

    Raises:
        argparse.ArgumentTypeError: If the input does not match the required format.
    """
    fasta_extensions = {".fasta", ".fna", ".fa"}
    if Path(value).suffix in fasta_extensions:
        validate_files(value)
        file_log.info(f"Validated reference genome file: {value}")
        return value

    accession_patterns = [
        r"^[A-Z]{1,2}\d{5,6}$",
        r"^[A-Z]{2}_\d{6,}\.\d+$",
        r"^[A-Z]{3}_\d{9}\.\d+$",
        r"^[A-Z]{2}\d{6,}\.\d+$",
    ]

    if any(re.match(pattern, value) for pattern in accession_patterns):
        file_log.info(f"Validated reference genome accession: {value}")
        return value

    file_log.error(
        f"Invalid reference genome input: '{value}'. Must be a .fasta or .fna file, or a valid accession code.")
    raise argparse.ArgumentTypeError(
        f"Invalid reference genome input: '{value}'. Must be a .fasta or .fna file, or a valid accession code.")


def validate_flanking_genes(value):
    """
    Validates a comma-separated list of flanking genes, ensuring that the
    list is composed of uppercase gene names and contains an even number
    of elements. If the list is odd, an error is logged and an exception
    is raised.

    Args:
        value (str): Comma-separated list of flanking genes.

    Returns:
        list: A list of validated flanking gene names.

    Raises:
        argparse.ArgumentTypeError: If the number of genes is not even.
    """
    flanking_genes = [gene.strip().upper() if gene.strip() != '-' else '-' for gene in value.split(',')]
    if len(flanking_genes) % 2 == 1:
        file_log.error(f"The specified flanking genes: {flanking_genes} should be even numbers.")
        raise argparse.ArgumentTypeError(f"The specified flanking genes: {flanking_genes} should be even numbers (e.g., 2, 4, 6, 8) rather than odd (e.g., 1, 3, 5).")
    file_log.info(f"Validated flanking genes: {flanking_genes}")
    return value


def validate_chromosome(value):
    """
    Validates a comma-separated list of chromosome identifiers, ensuring
    that each identifier is a valid chromosome number between 1-22, or
    'X' or 'Y'. If any identifier is invalid, an error is logged and an
    exception is raised.

    Args:
        value (str): Comma-separated list of chromosome identifiers.

    Returns:
        list: A list of validated chromosome identifiers.

    Raises:
        argparse.ArgumentTypeError: If any chromosome identifier is invalid.
    """
    chromosomes = [chromosome.strip() for chromosome in value.split(',')]
    valid_chromosomes = set(map(str, range(1, 23))) | {"X", "Y"}

    if all(chromosome in valid_chromosomes for chromosome in chromosomes):
        file_log.info(f"Validated chromosomes: {chromosomes}")
        return chromosomes
    else:
        file_log.error(f"Invalid chromosome list: '{value}'. All values must be integers between 1-22, or 'X', 'Y'.")
        raise argparse.ArgumentTypeError(f"Invalid chromosome list: '{value}'. All values must be integers between 1-22, or 'X', 'Y'.")


def setup_pipeline_args(subparsers):
    """
    Configures the argument parser for the 'pipeline' command. This command
    processes sequencing data, and the function sets up various argument groups
    including reads, reference genome, and analysis settings. It also sets up
    validation to ensure that the '--default' option is mutually exclusive with
    '--flanking-genes' and '--chromosomes'.

    Args:
        subparsers (argparse._SubParsersAction): The subparsers object to add
        the 'pipeline' command to.

    Raises:
        argparse.ArgumentError: If invalid argument combinations are provided.
    """
    parser_pipeline = subparsers.add_parser('pipeline', help='Run the pipeline for sequencing data processing.')

    reads_group = parser_pipeline.add_argument_group('reads', 'Arguments related to read files')
    reads_group.add_argument('-ont', '--nanopore', required=True, type=validate_read_files, help='Path to the Oxford Nanopore Technologies reads file.')
    reads_group.add_argument('-pb', '--pacbio', required=True, type=validate_read_files, help='Path to the Pacific Biosciences reads file.')

    reference_group = parser_pipeline.add_argument_group('reference', 'Arguments related to the reference genome')
    reference_group.add_argument('-ref', '--reference', type=validate_reference, required=False, help='Path to the reference genome file if one is already present (optional).')

    analysis_group = parser_pipeline.add_argument_group('analysis', 'Arguments related to analysis settings')

    exclusive_group = analysis_group.add_mutually_exclusive_group(required=False)
    exclusive_group.add_argument('--default', action='store_true', help='Use default settings. Cannot be used with --flanking-genes or --chromosomes.')

    analysis_group.add_argument('-M', '--metadata', type=validate_file, help='Directory containing the metadata file relevant to the analysis. (.XLMX)')
    analysis_group.add_argument('-f', '--flanking-genes', type=validate_flanking_genes, help='Comma-separated list of flanking genes, e.g., MGAM2,EPHB6. Add them as pairs.')
    analysis_group.add_argument('-c', '--chromosomes', type=validate_chromosome, help='List of chromosomes where TR or IG is located.')
    analysis_group.add_argument('-s', '--species', type=str.capitalize, required=True, help='Species name, e.g., Homo sapiens.')
    analysis_group.add_argument('-r', '--receptor-type', required=True, type=str.upper, choices=['TR', 'IG'], help='Type of receptor to analyze: TR (T-cell receptor) or IG (Immunoglobulin).')
    analysis_group.add_argument('-t', '--threads', type=int, required=False, default=8, help='Amount of threads to run the analysis.')
    analysis_group.add_argument('-o', '--output', type=str, default=str(Path.cwd() / 'annotation_results'), help='Output directory for the results.')

    parser_pipeline.set_defaults(func=run_pipeline)

    def validate_pipeline_args(args):
        """
        Validates the arguments provided to the pipeline command, ensuring
        that '--default' is not used together with '--flanking-genes' or
        '--chromosomes', and that both '--flanking-genes' and '--chromosomes'
        are provided together unless '--default' is used.

        Args:
            args (argparse.Namespace): The parsed arguments.

        Raises:
            argparse.ArgumentError: If invalid argument combinations are provided.
        """
        if args.default:
            if args.flanking_genes or args.chromosomes:
                parser_pipeline.error("--default cannot be used with --flanking-genes or --chromosomes.")
        else:
            if args.flanking_genes is None or args.chromosomes is None:
                parser_pipeline.error("Both --flanking-genes and --chromosomes must be provided together unless --default is used.")

    parser_pipeline.set_defaults(validate_pipeline_args=validate_pipeline_args)


def setup_annotation_args(subparsers):
    """
    Configures the argument parser for the 'annotation' command. This command
    is used for VDJ segment analysis, and the function sets up various argument
    groups including the library, receptor type, and species. It also sets up
    validation to ensure that the '--default' option is mutually exclusive with
    '--flanking-genes'.

    Args:
        subparsers (argparse._SubParsersAction): The subparsers object to add
        the 'annotation' command to.
    """
    parser_annotation = subparsers.add_parser('annotation', help='Run the annotation tool for VDJ segment analysis.')

    group = parser_annotation.add_mutually_exclusive_group()
    group.add_argument('--default', action='store_true', help='Use default settings. Cannot be used with --flanking-genes.')
    group.add_argument('-f', '--flanking-genes', type=validate_flanking_genes, help='Comma-separated list of flanking genes, e.g., MGAM2,EPHB6. Add them as pairs.')

    parser_annotation.add_argument('-l', '--library', type=validate_file, help='Path to the library file. Expected to be in FASTA format.')
    parser_annotation.add_argument('-r', '--receptor-type', required=True, type=str.upper, choices=['TR', 'IG'], help='Type of receptor to analyze: TR (T-cell receptor) or IG (Immunoglobulin).')

    data_choice = parser_annotation.add_mutually_exclusive_group(required=True)
    data_choice.add_argument('-i', '--input', type=validate_input, help='Directory containing the extracted sequence regions in FASTA format, where VDJ segments can be found. Cannot be used with -f/--flanking-genes.')
    data_choice.add_argument('-a', '--assembly', type=validate_input, help='Directory containing the assembly FASTA files. Must be used with -f/--flanking-genes and -s/--species.')

    parser_annotation.add_argument('-S', '--scaffolding',required=False, type=validate_file, help='Path to the reference genome (FASTA) containing the chromosomes of interest for the selected species')

    parser_annotation.add_argument('-M', '--metadata',required=False, type=validate_file, help='Directory containing the metadata file relevant to the analysis. (.XLMX)')
    parser_annotation.add_argument('-s', '--species', type=str, help='Species name, e.g., Homo sapiens. Required with -a/--assembly.')
    parser_annotation.add_argument('-o', '--output', type=str, default=str(Path.cwd() / 'annotation_results'), help='Output directory for the results.')
    mapping_options = ['minimap2', 'bowtie', 'bowtie2']
    parser_annotation.add_argument('-m', '--mapping-tool', nargs='*', choices=mapping_options, default=mapping_options, help='Mapping tool(s) to use. Choose from: minimap2, bowtie, bowtie2. Defaults to all.')
    parser_annotation.add_argument('-t', '--threads', type=int, required=False, default=8, help='Amount of threads to run the analysis.')

    parser_annotation.set_defaults(func=run_annotation)


def validate_html(value):
    input_path = Path(value).resolve()
    if input_path.name == 'flask':
        validate_directory(str(input_path))
    else:
        input_path = input_path / 'flask'
        validate_directory(str(input_path))
    return input_path


def setup_html(subparsers):
    parser_html = subparsers.add_parser('html', help='Display the generated HTML report.')

    parser_html.add_argument('--show', action='store_true', required=True, help='Open the HTML report in your default web browser.')
    parser_html.add_argument('-i', '--input', required=True, type=validate_html, help='Path to the directory containing the HTML report.')
    parser_html.add_argument('--reset-flask', required=False, default=False, action='store_true', help="Reset the flask directory.")
    parser_html.add_argument('--dev-mode', required=False, default=False, action='store_true', help='Enable developer mode to see the flask output.')
    parser_html.set_defaults(func=run_html)


def run_pipeline(args):
    """
    Executes the main pipeline process for sequencing data processing. This
    function orchestrates various steps such as filtering and moving files,
    splitting chromosomes, running Snakemake for the pipeline, and generating
    HTML reports.

    Args:
        args (argparse.Namespace): The parsed arguments for the pipeline command.
    """

    try:
        settings_dir, output_dir = cwd_setup(args.output)
        os.chdir(output_dir)
        cwd = Path.cwd()
        config = {}
        final_output = cwd / 'final'
        file_log.info("Starting main process")

        filter_and_move_files(cwd, args.nanopore, args.pacbio, config)
        if args.reference and Path(args.reference).suffix in {'.fasta', '.fna', 'fa'}:
            fasta_path = cwd / args.reference
            split_chromosomes(cwd, fasta_path, config)
        else:
            genome_dir = cwd / 'reference' / 'genome'
            genome = download_reference_genome(args.reference, genome_dir)
            split_chromosomes(cwd, genome, config)
        if args.default:
            loop_flanking_genes(settings_dir, output_dir, args)
        create_config(output_dir, settings_dir, args)
        if not final_output.is_dir():
            snakefile = Path(settings_dir / 'Snakefile')
            run_snakemake(args, output_dir, str(snakefile))
    except Exception as e:
        file_log.error(f"Pipeline execution failed: {str(e)}")
    finally:
        cleanup(config)
        file_log.info("Cleanup completed")

    file_log.info("Main process completed successfully")


def run_annotation(args):
    """
    Executes the annotation process for VDJ segment analysis. This function
    handles tasks such as generating or validating the library file, creating
    configuration files, running the main annotation process, and generating
    HTML reports.

    Args:
        args (argparse.Namespace): The parsed arguments for the annotation command.
    """
    settings_dir, output_dir = cwd_setup(args.output)
    cwd = Path.cwd()
    file_log.info('Running the annotation program')
    library = cwd / 'library' / f'library.fasta'

    if not args.library:
        if not library.is_file():
            file_log.info(
                'No library specified, generating it with the IMGT scraper.')
            try:
                create_and_activate_env(settings_dir / 'envs' / 'IMGT.yaml')
                file_log.info(f"Starting scrape for species: {args.species}, type: {args.receptor_type}")
                imgt_main(species=args.species, immune_type=args.receptor_type)
                file_log.info(f"Downloaded the library from the IMGT for {args.species}")
            except subprocess.CalledProcessError as e:
                log_subprocess_error(e)
                exit()
            finally:
                deactivate_env()
                args.library = cwd / 'library' / library
        else:
            args.library = cwd / 'library' / library
    else:
        library_path = Path(args.library)
        destination = cwd / 'library' / "library.fasta"
        make_dir(cwd / 'library')
        shutil.copy(args.library, destination)
        args.library = cwd / 'library' / library

    create_config(output_dir, settings_dir, args)
    try:
        create_and_activate_env(settings_dir / 'envs' / 'scripts.yaml')
        annotation_main(args)
    except Exception as e:
        file_log.error(f"Annotation failed with error: {str(e)}")
    finally:
        deactivate_env()



def open_browser():
    webbrowser.open_new('http://127.0.0.1:5000/')


def get_python_executable():
    """
    Checks whether 'python' or 'python3' is available and returns the appropriate executable.

    Returns:
        str: The name of the Python executable ('python' or 'python3').

    Raises:
        RuntimeError: If neither 'python' nor 'python3' is found in the PATH.
    """
    if shutil.which("python"):
        return "python"
    elif shutil.which("python3"):
        return "python3"
    else:
        raise RuntimeError(
            "No suitable Python executable found. Ensure 'python' or 'python3' is available in PATH."
        )


def run_html(args):
    console_log.info(
        "Running the HTML report, which should automatically open in your browser.\n"
        "If it doesn't, try entering this address manually: http://127.0.0.1:5000.\n"
        "If that doesn't work, try http://localhost:8080.\n"
        "If neither address works, you may need to forward the port using SSH with the following command:\n"
        "'ssh -L 8080:localhost:5000 username@your_server_ip'.\n"
        "For further instructions, please refer to the README."
    )
    output_dir = args.input
    copy_flask(output_dir, args.reset_flask)
    os.chdir(output_dir.parent)
    try:
        if not args.dev_mode:
            threading.Timer(1, open_browser).start()
        python_executable = get_python_executable()
        process = subprocess.Popen(
            [python_executable, str(output_dir / 'app.py')],
            #stdout=False if args.dev_mode else subprocess.DEVNULL,
            #stderr=False if args.dev_mode else subprocess.DEVNULL,
        )
        process.communicate()
    except KeyboardInterrupt:
        console_log.info(
            "Process interrupted by user (Ctrl + C), HTML report closed!")


def key_sort(s):
    """
    Sorts strings with numbers in a natural order (e.g., "chr1", "chr2",
    "chr10" instead of "chr1", "chr10", "chr2").

    Args:
        s (str): String to sort.

    Returns:
        list: A list of integers and strings, sorted naturally.
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def split_chromosomes(cwd, fasta_path, config):
    """
    Splits the chromosomes from a reference genome FASTA file into
    individual files for each chromosome. The function creates an
    output directory for the chromosome files and a new reference
    genome file, and logs the process.

    Args:
        cwd (Path): The current working directory.
        fasta_path (Path): Path to the reference genome FASTA file.
    """
    output_dir = cwd / 'reference' / 'chromosomes'
    genome_dir = cwd / 'reference' / 'genome'
    new_genome_file = genome_dir / 'new_reference.fasta'
    if not new_genome_file.is_file():
        [make_dir(i) for i in [output_dir, genome_dir]]
        write_chromosomes(fasta_path, output_dir, new_genome_file, config)
    else:
        for chromosome in sorted(output_dir.glob("reference_chr*.fasta"), key=lambda x: key_sort(x.stem)):
            chromosome_number = chromosome.stem.replace('reference_chr', '')
            config.setdefault('ALL_CHROMOSOMES', []).append(chromosome_number)


def write_chromosomes(fasta_path, output_dir, new_genome_file, config):
    """
    Writes each chromosome from a reference genome FASTA file to
    individual files and a combined reference genome file. It logs
    the file paths and saves each chromosome as a separate FASTA file.

    Args:
        fasta_path (Path): Path to the reference genome FASTA file.
        output_dir (Path): Directory to save individual chromosome FASTA files.
        new_genome_file (Path): Path to the new combined reference genome file.
    """
    files = {}
    files["all"] = new_genome_file.open('a')
    for record in SeqIO.parse(fasta_path, 'fasta'):
        process_record(record, output_dir, files, config)

    for file in files.values():
        file.close()

    file_log.info(f"Files have been written to: {output_dir}")


def process_record(record, output_dir, files, config):
    """
    Processes a single record from a FASTA file, extracting chromosome
    information and writing it to the appropriate chromosome file. The
    function updates the config dictionary with chromosome numbers.

    Args:
        record (SeqRecord): A single record from the FASTA file.
        output_dir (Path): Directory to save the chromosome file.
        files (dict): A dictionary of open file handles for writing.
    """
    splitted = record.description.split(' ')
    try:
        chromosome, chromosome_number = extract_chromosome_info(splitted)
        if chromosome not in files:
            config.setdefault('ALL_CHROMOSOMES', []).append(
                chromosome_number.rstrip(","))
            file_path = output_dir / f"{chromosome}.fasta"
            file_path_txt = output_dir / f"{chromosome}.txt"
            files[chromosome] = file_path.open('a')
            files[f"{chromosome}_txt"] = file_path_txt.open('a')
        SeqIO.write(record, files[chromosome], 'fasta')
        files[f"{chromosome}_txt"].write(f"{record.id}\n")
        SeqIO.write(record, files['all'], 'fasta')
    except (StopIteration, IndexError) as e:
        file_log.warning(
            f"Error processing record: {record.description}. Error: {e}")


def extract_chromosome_info(splitted):
    """
    Extracts chromosome information from a split description string
    of a FASTA record. The function identifies the chromosome number
    and creates a chromosome identifier.

    Args:
        splitted (list): A list of strings split from the record description.

    Returns:
        str, str: The chromosome identifier and chromosome number.
    """
    chromosome, chromosome_number = next(
        (f"reference_chr{splitted[i+1].rstrip(',')}", splitted[i+1])
        for i, x in enumerate(splitted)
        if x == "chromosome"
    )
    return chromosome, chromosome_number


def rename_files(cwd: Path, file: Path, read_type: str):
    """
    Renames a sequencing read file to a standard format based on the
    read type (e.g., nanopore or pacbio) and moves it to the downloads
    directory. The function returns the sample name and the new file path.

    Args:
        cwd (Path): The current working directory.
        file (Path): The original file path of the sequencing read.
        read_type (str): The type of sequencing read (e.g., nanopore, pacbio).

    Returns:
        str, Path: The sample name and the new file path.
    """
    moved_dir = cwd / 'downloads'
    make_dir(moved_dir)
    new_file = get_new_file_path(file, read_type, moved_dir)
    file_log.info(f"Renamed file {file} to {new_file}")
    return file.stem.split('_')[0].upper(), new_file


def get_new_file_path(file, read_type, moved_dir):
    """
    Constructs a new file path for a sequencing read file based on the
    read type and sample name. The file is moved to the downloads directory.

    Args:
        file (Path): The original file path.
        read_type (str): The type of sequencing read (e.g., nanopore, pacbio).
        moved_dir (Path): The directory to move the renamed file to.

    Returns:
        Path: The new file path in the downloads directory.
    """
    conversion = {"ONT": "nanopore",
                  "NANOPORE": "nanopore",
                  "PACBIO": "pacbio",
                  "PB": "pacbio"
                  }
    stripped_extensions = file.with_suffix('').with_suffix('')
    sample = stripped_extensions.stem.split('_')[0].upper()
    new_sample_name = sample if sample not in conversion else "SAMPLE"
    return moved_dir / f"{new_sample_name}_{read_type}.fastq.gz"


def filter_and_move_files(cwd: Path, original_nanopore: Path, original_pacbio: Path, config):
    """
    Filters and moves sequencing read files (nanopore and pacbio) to
    standard locations in the working directory. It updates the config
    dictionary with the original and moved file paths.

    Args:
        cwd (Path): The current working directory.
        original_nanopore (Path): Path to the original nanopore read file.
        original_pacbio (Path): Path to the original pacbio read file.
    """
    ont_sample, moved_nanopore = rename_files(
        cwd, original_nanopore, 'nanopore')
    pb_sample, moved_pacbio = rename_files(cwd, original_pacbio, 'pacbio')
    config.setdefault('DATA', {
        'ORIGINAL': {
            'pacbio': str(original_pacbio),
            'nanopore': str(original_nanopore)
        },
        'MOVED': {
            'nanopore': str(moved_nanopore),
            'pacbio': str(moved_pacbio)
        }
    })
    if not moved_nanopore.is_file() and not moved_pacbio.is_file():
        move_files(original_nanopore, original_pacbio, moved_nanopore,
                   moved_pacbio, ont_sample, pb_sample)


def move_files(original_nanopore, original_pacbio, moved_nanopore, moved_pacbio, ont_sample, pb_sample):
    """
    Moves the original nanopore and pacbio read files to their new
    locations in the downloads directory. If the sample names differ,
    the files are renamed with a standard naming convention.

    Args:
        original_nanopore (Path): Path to the original nanopore read file.
        original_pacbio (Path): Path to the original pacbio read file.
        moved_nanopore (Path): Path to the new nanopore read file location.
        moved_pacbio (Path): Path to the new pacbio read file location.
        ont_sample (str): The sample name extracted from the nanopore read file.
        pb_sample (str): The sample name extracted from the pacbio read file.
    """
    if ont_sample != pb_sample:
        moved_nanopore = moved_nanopore.with_name(f"SAMPLE_ont.fastq.gz")
        moved_pacbio = moved_pacbio.with_name(f"SAMPLE_pacbio.fastq.gz")
    shutil.move(original_nanopore, moved_nanopore)
    shutil.move(original_pacbio, moved_pacbio)
    file_log.info(
        f"Moved {original_nanopore} to {moved_nanopore} and {original_pacbio} to {moved_pacbio}")


def download_flanking_genes(gene, dir: Path, species):
    """
    Downloads the flanking genes for a given gene and species using
    the NCBI datasets command-line tool. The downloaded files are
    extracted and processed into a FASTA file.

    Args:
        gene (str): The name of the gene for which flanking genes are to be downloaded.
        dir (Path): The directory to save the downloaded files.
        species (str): The species name for which the genes are being downloaded.
    """
    dir = Path(dir)
    output_zip = dir / f"{gene}.zip"
    output_fna = dir / f"{gene}.fna"
    if not output_fna.is_file():
        run_download_command(gene, species, output_zip, output_fna, dir)


def run_download_command(gene, species, output_zip, output_fna, dir):
    """
    Executes the command to download flanking genes for a given gene
    and species. The function logs the command and processes the
    downloaded files.

    Args:
        gene (str): The name of the gene for which flanking genes are to be downloaded.
        species (str): The species name for which the genes are being downloaded.
        output_zip (Path): The path to save the downloaded zip file.
        output_fna (Path): The path to save the processed FASTA file.
        dir (Path): The directory to save the downloaded and processed files.

    Raises:
        subprocess.CalledProcessError: If the command execution fails.
        Exception: If any unexpected error occurs during file processing.
    """
    try:
        command = f'datasets download gene symbol {gene} --taxon "{species.capitalize()}" --include gene --filename {output_zip}'
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process_downloaded_files(output_zip, dir, output_fna)
        file_log.info(f"Downloaded and processed flanking genes for {gene}")
    except subprocess.CalledProcessError as e:
        log_subprocess_error(e)
    except Exception as e:
        file_log.error(f"An unexpected error occurred: {e}")


def process_downloaded_files(output_zip, dir, output_fna):
    """
    Processes the downloaded zip file by extracting its contents and
    renaming the gene FASTA file. The function then cleans up the
    downloaded files and logs the process.

    Args:
        output_zip (Path): The path to the downloaded zip file.
        dir (Path): The directory where the files are extracted.
        output_fna (Path): The path to save the processed FASTA file.
    """
    unzip_file(output_zip, dir)
    ncbi_fna_path = dir / "ncbi_dataset" / "data" / "gene.fna"
    ncbi_fna_path.rename(output_fna)
    cleanup_downloaded_files(output_zip, dir)


def cleanup_downloaded_files(output_zip, dir):
    """
    Cleans up the downloaded files by deleting the zip file, README file,
    and the extracted directory. The function logs the cleanup process.

    Args:
        output_zip (Path): The path to the downloaded zip file.
        dir (Path): The directory where the files were extracted.
    """
    output_zip.unlink()
    readme_path = dir / "README.md"
    if readme_path.exists():
        readme_path.unlink()
    shutil.rmtree(dir / "ncbi_dataset")
    sleep(2)


def log_subprocess_error(e):
    """
    Logs the error information from a failed subprocess command. The
    function captures the standard output and error streams and logs
    them for debugging.

    Args:
        e (subprocess.CalledProcessError): The exception object from a failed subprocess command.
    """
    file_log.error(f"Error occurred while running command: {e}")
    #file_log.error(e.stdout.decode())
    #file_log.error(e.stderr.decode())


def download_reference_genome(genome_code, reference_dir: Path):
    """
    Downloads a reference genome using the NCBI datasets command-line tool,
    extracts the genome FASTA file, and saves it to the specified directory.

    Args:
        genome_code (str): The accession code for the reference genome.
        reference_dir (Path): The directory to save the downloaded genome.

    Returns:
        Path: The path to the downloaded and processed reference genome FASTA file.
    """
    reference_dir = Path(reference_dir)
    make_dir(reference_dir)
    output_zip = reference_dir / f"{genome_code}.zip"
    output_fna = reference_dir / f"reference.fna"
    if not output_fna.is_file():
        run_reference_download_command(genome_code, output_zip, reference_dir, output_fna)
    return output_fna


def run_reference_download_command(genome_code, output_zip, reference_dir, output_fna):
    """
    Executes the command to download a reference genome based on the
    accession code. The function logs the command and processes the
    downloaded files.

    Args:
        genome_code (str): The accession code for the reference genome.
        output_zip (Path): The path to save the downloaded zip file.
        reference_dir (Path): The directory to save the extracted genome files.
        output_fna (Path): The path to save the processed reference genome FASTA file.

    Raises:
        subprocess.CalledProcessError: If the command execution fails.
        Exception: If any unexpected error occurs during file processing.
    """
    try:
        command = f'datasets download genome accession {genome_code} --include genome --filename {output_zip}'
        file_log.info(f"Running command: {command}")
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process_reference_files(output_zip, reference_dir,output_fna, genome_code)
    except subprocess.CalledProcessError as e:
        log_subprocess_error(e)
    except Exception as e:
        file_log.error(f"An unexpected error occurred: {e}")


def process_reference_files(output_zip, reference_dir, output_fna, genome_code):
    """
    Processes the downloaded reference genome zip file by extracting its
    contents and renaming the genome FASTA file. The function then cleans
    up the downloaded files and logs the process.

    Args:
        output_zip (Path): The path to the downloaded zip file.
        reference_dir (Path): The directory where the files were extracted.
        output_fna (Path): The path to save the processed reference genome FASTA file.
        genome_code (str): The accession code for the reference genome.
    """
    unzip_file(output_zip, reference_dir)
    ncbi_fna_path = (reference_dir / "ncbi_dataset" /
                     "data" / genome_code).glob('*.fna')
    ncbi_fna_path = list(ncbi_fna_path)[0]
    ncbi_fna_path.rename(output_fna)
    cleanup_downloaded_files(output_zip, reference_dir)
    file_log.info(f"Downloaded the reference genome: {genome_code}")


def deep_merge(d1, d2):
    """
    Deep merges two dictionaries. If there are nested dictionaries, the
    function recursively merges them. If there are conflicting keys, the
    values from the second dictionary overwrite those in the first.

    Args:
        d1 (dict): The first dictionary to merge.
        d2 (dict): The second dictionary to merge.
    """
    for key, value in d2.items():
        if key in d1:
            if isinstance(d1[key], dict) and isinstance(value, dict):
                deep_merge(d1[key], value)
            else:
                d1[key] = value
        else:
            d1[key] = value


def loop_flanking_genes(settings_dir, output_dir, args):
    """
    Processes flanking genes by downloading them and determining their
    associated chromosomes. The function updates the argument namespace
    with the downloaded flanking genes and their chromosomes.

    Args:
        cwd (Path): The current working directory.
        args (argparse.Namespace): The parsed arguments for the pipeline command.
    """
    receptor_chromosomes = set()
    flanking_genes_dir = output_dir / "flanking_genes"
    make_dir(flanking_genes_dir)
    species_dict = get_species_dict(settings_dir, args)
    flanking_genes = species_dict.get(
        args.receptor_type, {}).get('FLANKING_GENES', [])
    args.flanking_genes = flanking_genes
    for gene in args.flanking_genes:
        process_flanking_gene(gene, flanking_genes_dir,
                              args.species, receptor_chromosomes)
    args.chromosomes = list(receptor_chromosomes)


def get_species_dict(settings_dir, args):
    """
    Retrieves the species-specific configuration dictionary from the
    species YAML file. If a species-specific configuration is not found,
    the default configuration is returned.

    Args:
        cwd (Path): The current working directory.
        args (argparse.Namespace): The parsed arguments for the pipeline command.

    Returns:
        dict: The species-specific or default configuration dictionary.
    """
    default_dict = load_config(settings_dir / '_config' / 'species.yaml')
    species_key = args.species.replace(' ', '_')
    return default_dict.get(species_key, default_dict.get('default', {}))


def process_flanking_gene(gene, flanking_genes_dir, species, receptor_chromosomes):
    """
    Downloads the flanking gene for a given species and adds the
    associated chromosome information to the receptor_chromosomes set.

    Args:
        gene (str): The name of the gene for which flanking genes are to be downloaded.
        flanking_genes_dir (Path): The directory to save the downloaded files.
        species (str): The species name for which the genes are being downloaded.
        receptor_chromosomes (set): A set to store the associated chromosomes.
    """
    download_flanking_genes(gene, flanking_genes_dir, species)
    gene_file = flanking_genes_dir / f"{gene}.fna"
    records = SeqIO.to_dict(SeqIO.parse(gene_file, 'fasta'))
    first_key = next(iter(records))
    first_record = records[first_key]
    chromosome_info = first_record.description.split('=')[-1].strip(']')
    number, trailing = extract_chromosome_number_and_trailing(chromosome_info)
    receptor_chromosomes.add(number + trailing)


def extract_chromosome_number_and_trailing(chromosome_info):
    """
    Extracts the chromosome number and any trailing characters
    (e.g., 'X', 'Y') from a chromosome description string.

    Args:
        chromosome_info (str): The chromosome description string.

    Returns:
        str, str: The chromosome number and trailing characters.
    """
    number = ''.join(filter(str.isdigit, chromosome_info))
    trailing = ''.join(filter(str.isalpha, chromosome_info))
    return number, trailing


def create_config(output_dir, settings_dir, args):
    """
    Creates a configuration file based on the provided arguments
    and species-specific settings. The configuration is saved as a
    YAML file in the config directory.

    Args:
        cwd (Path): The current working directory.
        args (argparse.Namespace): The parsed arguments for the pipeline or annotation command.
    """
    config = {}
    species_config = load_config(settings_dir / '_config' / 'species.yaml')
    initialize_config(args, species_config, settings_dir, config)
    rss_config = load_config(settings_dir / '_config' / 'rss.yaml')
    deep_merge(config, rss_config.get(args.receptor_type, {}))
    config_file = output_dir / 'tmp/config' / 'config.yaml'
    make_dir(config_file.parent)
    save_config_to_file(config_file, config)


def initialize_config(args, species_config, settings_dir, config):
    """
    Initializes the global config dictionary based on the provided
    arguments and species-specific settings. It populates the
    config dictionary with species name, receptor type, flanking genes,
    and other relevant information.

    Args:
        args (argparse.Namespace): The parsed arguments for the pipeline or annotation command.
        species_config (dict): The species-specific configuration dictionary.
    """
    config.setdefault('SPECIES', {
        'name': args.species,
        'cell': args.receptor_type
    })
    config["SETTINGS"] = str(settings_dir)
    if hasattr(args, 'reference') and args.reference is not None:
        config['SPECIES']['genome'] = args.reference
        config.setdefault("BUSCO_DATASET", "primates_odb10")
        config.setdefault("HAPLOTYPES", [1, 2])
    if hasattr(args, 'chromosomes') and args.chromosomes is not None:
        config.setdefault("ASSEMBLY_CHROMOSOMES", sorted(
            args.chromosomes, key=lambda x: key_sort(x)))
    if hasattr(args, 'flanking_genes') and args.flanking_genes is not None:
        config.setdefault("FLANKING_GENES", args.flanking_genes)

    # Handling flanking genes based on default settings
    flanking_genes = getattr(args, 'flanking_genes', None)
    if not flanking_genes and args.default:
        species = args.species.capitalize().replace(" ", "_")
        receptor = args.receptor_type

        species_data = species_config.get(species, species_config.get("default", {}))
        receptor_data = species_data.get(receptor, species_config["default"].get(receptor, {}))
        flanking_genes = receptor_data.get("FLANKING_GENES", {})

        #flanking_genes = species_config.get(
        #    args.species.capitalize().replace(" ", "_"), "default"
        #).get(args.receptor_type, {}).get("FLANKING_GENES")

    if flanking_genes:
        config.setdefault("FLANKING_GENES", flanking_genes)
        args.flanking_genes = flanking_genes

    if hasattr(args, 'assembly') and args.assembly is not None:
        config.setdefault('DATA', {})['library'] = str(args.library)
        config.setdefault('DATA', {})['assembly'] = str(args.assembly)


def save_config_to_file(file_name, config):
    """
    Saves the global config dictionary to a YAML configuration file.

    Args:
        file_name (Path): The path to the YAML configuration file.
    """
    with open(file_name, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)


def run_snakemake(args, output_dir, config, snakefile="Snakefile"):
    """
    Runs the Snakemake workflow using the provided number of threads.
    The function first performs a dry-run to check the workflow and
    then runs the complete pipeline if the dry-run is successful.

    Args:
        args (argparse.Namespace): The parsed arguments for the pipeline command.
        snakefile (str): The path to the Snakefile. Default is "Snakefile".

    Raises:
        SystemExit: If the Snakemake workflow fails.
    """
    config_file = Path(output_dir / 'config' / 'config.yaml')
    try:
        result = subprocess.run(
            f"snakemake -s {snakefile} --cores {args.threads} --use-conda --configfile {config_file} -prn",
            shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        file_log.info("Snakemake check ran successfully.")
        file_log.info("Running the complete pipeline.")
        subprocess.run(
            f"snakemake -s {snakefile} --cores {args.threads} --use-conda --configfile {config_file} -pr",
            shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        file_log.error(f"Snakemake failed with return code {e.returncode}.")
        file_log.error(f"Snakemake output: {e.stdout.decode()}")
        file_log.error(f"Application is shutting down.")
        cleanup(config)
        sys.exit(e.returncode)


def cleanup(config):
    """
    Cleans up by moving the original sequencing read files back
    to their original locations from the downloads directory,
    and removes the downloads directory.

    """
    pacbio_original = config['DATA']['ORIGINAL'].get('pacbio')
    nanopore_original = config['DATA']['ORIGINAL'].get('nanopore')
    nanopore_moved = config['DATA']['MOVED'].get('nanopore')
    pacbio_moved = config['DATA']['MOVED'].get('pacbio')
    shutil.move(nanopore_moved, nanopore_original)
    shutil.move(pacbio_moved, pacbio_original)
    Path('downloads').rmdir()


def main():
    """
    The main entry point of the script. It sets up the argument
    parsers for the 'pipeline' and 'annotation' commands, parses
    the command-line arguments, and invokes the corresponding
    function based on the command provided.

    Raises:
        SystemExit: If no command is provided or if invalid arguments are provided.
    """
    console_log.info("Starting the VDJ Insights pipeline.")

    parser = argparse.ArgumentParser(description="Tool for sequencing data processing and VDJ annotation")
    subparsers = parser.add_subparsers(dest='command', help='Available commands', required=True)

    setup_pipeline_args(subparsers)
    setup_annotation_args(subparsers)
    setup_html(subparsers)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if args.command == 'pipeline':
        args.validate_pipeline_args(args)
    if not args.command:
        parser.print_help()
    else:
        args.func(args)


if __name__ == '__main__':
    main()
