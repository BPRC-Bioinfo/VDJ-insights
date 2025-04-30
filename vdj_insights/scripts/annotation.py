"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""
import ast

from tqdm import tqdm

from .scaffolding import scaffolding_main

from .mapping import mapping_main
from .report import report_main

from .functionality import main_functionality
from .RSS import main_rss
from .CDR import main_cdr

from pathlib import Path
import pandas as pd
import argparse
from .blast import blast_main
from .map_genes import map_main
from .extract_region import region_main

from .report import make_bed, make_gtf

from .figures.barplot import main as barplot_main
from .figures.boxplot import main as boxplot_main
from .figures.broken_regions import main as broken_regions_main
from .figures.heatmap import main as heatmap_main
from .figures.sub_families import main as sub_families_main
from .figures.venn_diagram import main as venn_diagram_main

from .util import make_dir, validate_file, validate_input, validate_metadata_coverage
from .property import log_error

from .logger import console_logger, file_logger


console_log = console_logger(__name__)
file_log = file_logger(__name__)


@log_error()
def validate_flanking_genes(value: str) -> list:
    """
    Validates and processes a list of flanking genes, ensuring they are uppercase and handles empty values.
    Splits a comma-separated string of flanking genes, strips whitespace, converts them to uppercase, and
    handles any empty values (represented by '-'). Ensures that the number of flanking genes is even; if not,
    raises an `argparse.ArgumentTypeError`.

    Args:
        value (str): Comma-separated string of flanking genes.

    Returns:
        list: List of processed flanking genes.

    Raises:
        argparse.ArgumentTypeError: If the number of flanking genes is odd.
    """
    flanking_genes = [gene.strip().upper() for gene in value.split(',')]
    if len(flanking_genes) % 2 == 1:
        raise argparse.ArgumentTypeError(f"""The specified flanking genes: {flanking_genes} should be even numbers (e.g., 2, 4, 6, 8) rather than odd (e.g., 1, 3, 5).""")
    return value


def argparser_setup(include_help: bool = True) -> argparse.ArgumentParser:
    """
    Configures the argument parser for the annotation tool. Sets up the command-line argument parser, defining
    the required and optional arguments, mutually exclusive groups, and their validation logic. Supports customization
    of the help option.

    Args:
        include_help (bool): Whether to include the help argument. Defaults to True.

    Returns:
        argparse.ArgumentParser: Configured argument parser ready for parsing input arguments.
    """
    parser = argparse.ArgumentParser(description='A tool for finding known and novel VDJ segments in certain data, with a library of choice.',formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=include_help)

    required_group = parser.add_argument_group('Required Options')
    required_group.add_argument('-l', '--library', required=True, type=validate_file, help='Path to the library file. Expected to be in FASTA format.')
    required_group.add_argument('-r', '--receptor-type', required=True, type=str.upper, choices=['TR', 'IG'], help='Type of receptor to analyze: TR (T-cell receptor) or IG (Immunoglobulin).')

    regions_or_assembly_group = parser.add_argument_group('Data Source Options','Select the data source: regions or assembly. These options are mutually exclusive.')
    data_choice = regions_or_assembly_group.add_mutually_exclusive_group(required=True)
    data_choice.add_argument('-i', '--input', type=validate_input, help='Directory containing the extracted sequence regions in FASTA format, where VDJ segments can be found. Cannot be used with -f/--flanking-genes or -s/--species.')
    data_choice.add_argument('-a', '--assembly', type=validate_input, help='Directory containing the assembly FASTA files. Must be used with -f/--flanking-genes and -s/--species.')

    assembly_options = parser.add_argument_group('Assembly-Specific Options','These options are required if -a/--assembly is chosen:')
    assembly_options.add_argument('-M', '--metadata',required=False, type=validate_file, help='Directory containing the metadata file relevant to the analysis. (.XLMX)')
    assembly_options.add_argument('-s', '--species', type=str.capitalize, help='Species name, e.g., Homo sapiens. Required with -a/--assembly.')
    assembly_options.add_argument('-S', '--scaffolding',required=False, type=validate_file, help='Path to the reference genome (FASTA) containing the chromosomes of interest for the selected species')

    exclusive_group = parser.add_argument_group('Exclusive Options')
    exclusive_mutually_exclusive = exclusive_group.add_mutually_exclusive_group()
    exclusive_mutually_exclusive.add_argument('-f', '--flanking-genes', type=validate_flanking_genes, help='Comma-separated list of flanking genes, e.g., MGAM2,EPHB6. Add them as pairs. Required with -a/--assembly.')
    exclusive_mutually_exclusive.add_argument('--default', action='store_true', help='Use default settings. Cannot be used with -f/--flanking-genes or -c/--chromosomes.')

    optional_group = parser.add_argument_group('Optional Options')
    optional_group.add_argument('-o', '--output', type=str, default='annotation', help='Output directory for the results.')
    mapping_options = ['minimap2', 'bowtie', 'bowtie2']
    optional_group.add_argument('-m', '--mapping-tool', nargs='*', choices=mapping_options, default=mapping_options, help='Mapping tool(s) to use. Choose from: minimap2, bowtie, bowtie2. Defaults to all.')
    optional_group.add_argument('-t', '--threads', type=int,required=False, default=8, help='Number of threads to run the analysis.')

    return parser

import time


@log_error()
def main(args=None):
    """
    Main function for the annotation tool. Performs the following steps:

        1. Parses command-line arguments and sets up paths for input and output files.
        2. Validates input parameters and configures the environment based on the selected options.
        3. Creates necessary directories and initializes the annotation process.
        4. Executes the mapping and region extraction processes if assembly mode is selected:
           - Calls `map_main` to map flanking genes against the assembly.
           - Calls `region_main` to extract regions based on mapped genes.
        5. Retrieves or creates a DataFrame from existing or new mapping results:
           - Uses `get_or_create` to check for existing reports or generate new ones.
           - Calls `combine_df` to merge mapping results and remove duplicates.
        6. Runs BLAST operations to align sequences and generate results:
           - Calls `blast_main` to perform BLAST alignment on the combined DataFrame.
        7. Generates the final annotation report and RSS file:
           - Calls `report_main` to generate an Excel report summarizing the findings.
           - Calls `RSS_main` to produce an RSS feed if required.

    Args:
        args (list or None): Command-line arguments to parse. Defaults to None.

    Returns:
        None

    Raises:
        ValueError: If invalid arguments are passed to the main function.
        OSError: If any file or directory operation fails.
        subprocess.CalledProcessError: If a BLAST command or database creation fails.
    """
    cwd = Path.cwd()

    console_log.info(f"Initialise pipeline")
    update_args = argparser_setup()

    if args.assembly and args.metadata:
        if not validate_metadata_coverage(args.assembly, args.metadata):
            console_log.error("Shutting down script. Metadata does not cover all assembly files.")
            exit()

    if args is None:
        args = update_args.parse_args()
    elif isinstance(args, list):
        args = update_args.parse_args(args)
    elif not isinstance(args, argparse.Namespace):
        raise ValueError("Invalid arguments passed to the main function")
    if args.assembly:
        if not args.flanking_genes or not args.species:
            update_args.error('-a/--assembly requires -f/--flanking-genes and -s/--species.')
    if args.scaffolding:
        if not args.assembly:
            update_args.error('-S/--scaffolding requires -a/--assembly.')
    if args.input:
        if not args.receptor_type or not args.species:
            update_args.error('-i/--input requires -r/--region and -s/--species.')
    args.species = args.species.capitalize() if args.species else None

    region_dir = "tmp/region"
    if args.input:
        region_dir = args.input

    if args.scaffolding:
        scaffolding_dir = "tmp/scaffold_assemblies"
        args.assembly = scaffolding_main(args.scaffolding, args.assembly, scaffolding_dir,  args.threads)

    annotation_folder = cwd / 'annotation'
    make_dir(annotation_folder)

    timing_results = []

    #mapping flanking genes and region extraction
    flanking_genes_dict = ast.literal_eval(str(args.flanking_genes))
    if args.assembly:
        start = time.time()
        map_main(flanking_genes_dict, args.assembly, args.species, args.threads)
        end = time.time()
        timing_results.append(["Mapping flanking genes", round(end - start, 2)])

        start = time.time()
        region_main(flanking_genes_dict, args.assembly, args.threads)
        end = time.time()
        timing_results.append(["Region of intrest extraction", round(end - start, 2)])

    #mapping library and creating report
    report = annotation_folder / "tmp/report.csv"
    if not report.exists():
        report_df = pd.DataFrame()
        start = time.time()
        for tool in args.mapping_tool:
            file_log.info(f"Processing tool: {tool}")
            mapping_df = mapping_main(tool, args.receptor_type, region_dir, args.library, args.threads)
            report_df = pd.concat([report_df, mapping_df])
        report_df = report_df.drop_duplicates(subset=["reference", "start", "stop", "name"]).reset_index(drop=True)
        make_dir(annotation_folder / "tmp/")
        report_df.to_csv(report, index=False)
        end = time.time()
        timing_results.append([f"Mapping library", round(end - start, 2)])

    else:
        report_df = pd.read_csv(report)

    #reevaluate mapping genes of library with blast
    blast_file = annotation_folder / "tmp/blast_results.csv"
    if not blast_file.exists():
        start = time.time()
        blast_main(report_df, blast_file, args.library, args.threads)
        end = time.time()
        timing_results.append(["BLAST revaluation", round(end - start, 2)])

    #create report and rss
    start = time.time()
    report_main(annotation_folder, blast_file, args.receptor_type, args.library, args.assembly, args.metadata)
    end = time.time()
    timing_results.append(["Report and filtering", round(end - start, 2)])

    start = time.time()
    main_functionality(args.receptor_type, args.species, args.threads)
    end = time.time()
    timing_results.append(["Predict functionality", round(end - start, 2)])

    start = time.time()
    main_rss(args.threads)
    end = time.time()
    timing_results.append(["RSS extraction and validation", round(end - start, 2)])

    start = time.time()
    main_cdr(args.species, args.receptor_type, args.threads)
    end = time.time()
    timing_results.append(["CDR identification", round(end - start, 2)])

    #create BED and GTF files
    data = pd.read_excel(annotation_folder / "annotation_report_all.xlsx")
    make_bed(data, annotation_folder / "BED")
    make_gtf(data, annotation_folder / "GTF")

    #create figures
    if args.metadata:
        functions = [barplot_main, boxplot_main, broken_regions_main, sub_families_main, venn_diagram_main, heatmap_main]
        args_list = [(annotation_folder,), (annotation_folder,), (cwd,), (annotation_folder,), (annotation_folder, args.receptor_type), (annotation_folder,)]
        with tqdm(total=len(functions), desc="Creating plots", unit="Plot") as pbar:
            for func, args in zip(functions, args_list):
                start = time.time()
                func(*args)
                end = time.time()
                timing_results.append([func.__name__, round(end - start, 2)])

    excel_path = annotation_folder / "execution_times.xlsx"
    df = pd.DataFrame(timing_results, columns=["Process", "Execution Time (s)"])
    df.to_excel(excel_path, index=False)

    file_log.info(f"Annotation process completed. Results are available in {annotation_folder}")
    console_log.info(f"Annotation process completed. Results are available in {annotation_folder}")

