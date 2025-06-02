"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""
import argparse
import ast
import os
import shutil
import subprocess
import sys
import threading
import time
import webbrowser
from pathlib import Path

import pandas as pd
import yaml
from tqdm import tqdm

from .CDR import main_cdr
from .IMGT_scrape import main as imgt_main
from .RSS import main_rss
from .blast_v2 import blast_main
#from .blast import blast_main

from .functionality import main_functionality
from .logger import console_logger, file_logger
from .mapping import mapping_main
from .map_genes import map_main
from .extract_region import region_main
from .property import log_error
from .report import make_bed, make_gtf, report_main
from .scaffolding import scaffolding_main
from .util import load_config, make_dir,validate_directory, validate_file, validate_input, validate_metadata_coverage
from .env_manager import create_and_activate_env, deactivate_env

from .figures.barplot import main as barplot_main
from .figures.boxplot import main as boxplot_main
from .figures.broken_regions import main as broken_regions_main
from .figures.heatmap import main as heatmap_main
from .figures.sub_families import main as sub_families_main
from .figures.venn_diagram import main as venn_diagram_main

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def cwd_setup(output_dir):
    """
    Sets up the current working directory for the pipeline.

    Args:
        output_dir (str): Path to the output directory.

    Returns:
        tuple: Paths to the settings directory and output directory.
    """
    settings_dir = Path(__file__).resolve().parent.parent
    output_dir = Path(output_dir).resolve()
    make_dir(output_dir)
    copy_flask(output_dir / 'flask', settings_dir)
    os.chdir(str(output_dir))
    return settings_dir, output_dir


def copy_flask(output_dir, reset=False):
    """
    Copies the Flask directory to the specified output directory.

    Args:
        output_dir (Path): Path to the output directory.
        reset (bool): Whether to overwrite the existing Flask directory.
    """
    settings_dir = Path(__file__).resolve().parent.parent
    if not output_dir.is_dir() or reset:
        shutil.copytree(
            settings_dir / "flask",
            output_dir,
            dirs_exist_ok=True
        )


def validate_flanking_genes(value):
    """
    Validates the flanking genes argument.

    Args:
        value (str): Comma-separated list of flanking genes.

    Returns:
        str: Validated flanking genes.

    Raises:
        argparse.ArgumentTypeError: If the number of genes is not even.
    """
    genes = [
        g.strip().upper() if g.strip() != '-' else '-'
        for g in value.split(',')
    ]
    if len(genes) % 2:
        file_log.error(f"Flanking genes should be an even number: {genes}")
        raise argparse.ArgumentTypeError(
            f"Specify an even number of genes (e.g., 2, 4, 6): {genes}"
        )
    file_log.info(f"Validated flanking genes: {genes}")
    return value


def setup_annotation_args(subparsers):
    """
    Configures the command-line arguments for the annotation command.

    Args:
        subparsers (argparse._SubParsersAction): Subparsers object for adding commands.
    """
    p = subparsers.add_parser(
        'annotation',
        help='Run VDJ segment annotation.'
    )
    group = p.add_mutually_exclusive_group()
    group.add_argument(
        '--default',
        action='store_true',
        help='Use default settings (cannot be used with --flanking-genes).'
    )
    group.add_argument(
        '-f', '--flanking-genes',
        type=validate_flanking_genes,
        help='Comma-separated list of flanking genes (must be even count).'
    )

    p.add_argument(
        '-l', '--library',
        type=validate_file,
        help='Path to library FASTA.'
    )
    p.add_argument(
        '-r', '--receptor-type',
        required=True,
        type=str.upper,
        choices=['TR', 'IG'],
        help='TR or IG.'
    )

    data_choice = p.add_mutually_exclusive_group(required=True)
    data_choice.add_argument(
        '-i', '--input',
        type=validate_input,
        help='Directory with input FASTA regions.'
    )
    data_choice.add_argument(
        '-a', '--assembly',
        type=validate_input,
        help='Directory with assembly FASTA (requires --flanking-genes and --species).'
    )

    p.add_argument(
        '-S', '--scaffolding',
        type=validate_file,
        help='Reference genome FASTA.'
    )
    p.add_argument(
        '-M', '--metadata',
        type=validate_file,
        help='Metadata file directory.'
    )
    p.add_argument(
        '-s', '--species',
        type=str,
        help='Species name (required with --assembly).'
    )
    p.add_argument(
        '-o', '--output',
        type=str,
        default=str(Path.cwd() / 'annotation_results'),
        help='Output directory.'
    )
    p.add_argument(
        '-m', '--mapping-tool',
        nargs='*',
        choices=['minimap2', 'bowtie', 'bowtie2'],
        default=['minimap2', 'bowtie', 'bowtie2'],
        help='Mapping tools to use.'
    )
    p.add_argument(
        '-t', '--threads',
        type=int,
        default=8,
        help='Number of threads.'
    )
    p.set_defaults(func=run_annotation, parser=p)


def setup_html(subparsers):
    """
    Configures the command-line arguments for the HTML report command.

    Args:
        subparsers (argparse._SubParsersAction): Subparsers object for adding commands.
    """
    p = subparsers.add_parser(
        'html',
        help='Display the generated HTML report.'
    )
    p.add_argument(
        '-i', '--input',
        required=True,
        type=validate_html,
        help='Path to HTML report directory.'
    )
    p.add_argument(
        '--reset-flask',
        action='store_true',
        help='Re-copy the flask directory.'
    )
    p.add_argument(
        '--dev-mode',
        action='store_true',
        help='Show flask output in console.'
    )
    p.set_defaults(func=run_html)


def validate_html(value):
    """
    Validates the HTML report directory.

    Args:
        value (str): Path to the HTML report directory.

    Returns:
        Path: Validated path to the Flask directory.
    """
    p = Path(value).resolve()
    if p.name == 'flask':
        validate_directory(str(p))
    else:
        validate_directory(str(p / 'flask'))
        p = p / 'flask'
    return p


def run_annotation(args):
    """
    Executes the annotation pipeline.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    settings_dir, output_dir = cwd_setup(args.output)
    cwd = Path.cwd()
    file_log.info('Starting annotation')
    lib_dest = cwd / 'library' / 'library.fasta'

    if not args.library:
        if not lib_dest.is_file():
            create_and_activate_env(settings_dir / 'envs' / 'IMGT.yaml')
            imgt_main(species=args.species, immune_type=args.receptor_type)
            deactivate_env()
            args.library = lib_dest
        else:
            args.library = lib_dest
    else:
        make_dir(cwd / 'library')
        shutil.copy(args.library, lib_dest)
        args.library = lib_dest

    create_config(output_dir, settings_dir, args)
    create_and_activate_env(settings_dir / 'envs' / 'scripts.yaml')
    annotation_main(args)
    deactivate_env()


def open_browser():
    """
    Opens the default web browser to display the Flask application.
    """
    webbrowser.open_new('http://127.0.0.1:5000/')


def get_python_executable():
    """
    Finds the Python executable in the system PATH.

    Returns:
        str: Path to the Python executable.

    Raises:
        RuntimeError: If no Python executable is found.
    """
    for cmd in ("python", "python3"):
        if shutil.which(cmd):
            return cmd
    raise RuntimeError("No python executable found in PATH.")


def run_html(args):
    """
    Launches the HTML report using Flask.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    console_log.info("Launching HTML report...")
    output_dir = args.input
    copy_flask(output_dir, args.reset_flask)
    os.chdir(output_dir.parent)
    if not args.dev_mode:
        threading.Timer(1, open_browser).start()
    cmd = [get_python_executable(), str(output_dir / 'app.py')]
    subprocess.Popen(cmd).communicate()


def log_subprocess_error(e):
    """
    Logs errors from subprocess execution.

    Args:
        e (Exception): Exception object containing error details.
    """
    file_log.error(f"Subprocess failed: {e}")
    file_log.error(e.stdout.decode())


def deep_merge(d1, d2):
    """
    Recursively merges two dictionaries.

    Args:
        d1 (dict): Target dictionary.
        d2 (dict): Source dictionary.
    """
    for k, v in d2.items():
        if k in d1 and isinstance(d1[k], dict) and isinstance(v, dict):
            deep_merge(d1[k], v)
        else:
            d1[k] = v


def create_config(output_dir, settings_dir, args):
    """
    Creates the configuration file for the pipeline.

    Args:
        output_dir (Path): Path to the output directory.
        settings_dir (Path): Path to the settings directory.
        args (argparse.Namespace): Parsed command-line arguments.
    """
    config = {}
    species_cfg = load_config(settings_dir / '_config' / 'species.yaml')
    initialize_config(args, species_cfg, settings_dir, config)
    rss_cfg = load_config(settings_dir / '_config' / 'rss.yaml')
    deep_merge(config, rss_cfg.get(args.receptor_type, {}))
    cfg_path = output_dir / 'tmp' / 'config' / 'config.yaml'
    make_dir(cfg_path.parent)
    with open(cfg_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)


def initialize_config(args, species_cfg, settings_dir, config):
    """
    Initializes the configuration dictionary based on command-line arguments.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        species_cfg (dict): Species configuration data.
        settings_dir (Path): Path to the settings directory.
        config (dict): Configuration dictionary to populate.
    """
    parser = args.parser

    if args.default and not args.species:
        parser.error("--default requires --species/-s to be set")
    if args.assembly and not args.species:
        parser.error("--assembly requires --species/-s to be set")

    config['SPECIES'] = {'name': args.species, 'cell': args.receptor_type}
    config['SETTINGS'] = str(settings_dir)

    if getattr(args, 'flanking_genes', None):
        config['FLANKING_GENES'] = args.flanking_genes

    elif getattr(args, 'default', False):
        species_key = args.species.capitalize().replace(" ", "_")
        species_data = species_cfg.get(species_key, species_cfg.get("default", {}))
        receptor_data = species_data.get(
            args.receptor_type,
            species_cfg.get("default", {}).get(args.receptor_type, {})
        )
        default_flanking = receptor_data.get("FLANKING_GENES")
        if default_flanking:
            config['FLANKING_GENES'] = default_flanking
            args.flanking_genes = default_flanking

    if getattr(args, 'assembly', None):
        data = config.setdefault('DATA', {})
        data['library'] = str(args.library)
        data['assembly'] = str(args.assembly)


@log_error()
def annotation_main(args: argparse.Namespace):
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
    parser = args.parser

    cwd = Path.cwd()
    console_log.info(f"Initialise pipeline")

    if args.assembly and args.metadata:
        if not validate_metadata_coverage(args.assembly, args.metadata):
            console_log.error("Shutting down script. Metadata does not cover all assembly files.")
            exit()

    if args is None:
        args = args.parse_args()

    elif isinstance(args, list):
        args = args.parse_args(args)

    elif not isinstance(args, argparse.Namespace):
        raise ValueError("Invalid arguments passed to the main function")

    if args.assembly:
        if not args.flanking_genes or not args.species:
            parser.error('-a/--assembly requires -f/--flanking-genes and -s/--species.')

    if args.scaffolding:
        if not args.assembly:
            parser.error('-S/--scaffolding requires -a/--assembly.')

    if args.input:
        if not args.receptor_type or not args.species:
            parser.error('-i/--input requires -r/--region and -s/--species.')

    args.species = args.species.capitalize() if args.species else None

    region_dir = "tmp/region"
    if args.input:
        region_dir = args.input

    if args.scaffolding:
        scaffolding_dir = "tmp/scaffold_assemblies"
        args.assembly = scaffolding_main(args.scaffolding, args.assembly, scaffolding_dir,  args.threads)

    annotation_folder = cwd / 'annotation'
    make_dir(annotation_folder)

    #mapping flanking genes and region extraction
    flanking_genes_dict = ast.literal_eval(str(args.flanking_genes))
    if args.assembly:
        map_main(flanking_genes_dict, args.assembly, args.species, args.threads)
        region_main(flanking_genes_dict, args.assembly, args.threads)

    #mapping library and creating report
    report = annotation_folder / "tmp/report.csv"
    if not report.exists():
        report_df = pd.DataFrame()
        for tool in args.mapping_tool:
            file_log.info(f"Processing tool: {tool}")
            mapping_df = mapping_main(tool, region_dir, args.library, args.threads)
            report_df = pd.concat([report_df, mapping_df])
        report_df = report_df.drop_duplicates(subset=["reference", "start", "stop", "name"]).reset_index(drop=True)
        make_dir(annotation_folder / "tmp/")
        report_df.to_csv(report, index=False)

    else:
        report_df = pd.read_csv(report)

    #reevaluate mapping genes of library with blast
    blast_file = annotation_folder / "tmp/blast_results.csv"
    if not blast_file.exists():
        blast_main(report_df, blast_file, args.library, args.threads)

    #create report and rss
    report_main(annotation_folder, blast_file, args.receptor_type, args.library, args.assembly, args.metadata)
    main_functionality(args.receptor_type, args.species, args.threads)
    main_rss(args.threads)
    main_cdr(args.species, args.receptor_type, args.threads)

    #create BED and GTF files
    data = pd.read_excel(annotation_folder / "annotation_report_all.xlsx")
    make_bed(data, annotation_folder / "BED")
    make_gtf(data, annotation_folder / "GTF")

    #create figures
    if args.metadata:
        functions = [barplot_main, boxplot_main, sub_families_main, venn_diagram_main, heatmap_main]
        args_list = [(annotation_folder,), (annotation_folder,), (annotation_folder,), (annotation_folder, args.receptor_type), (annotation_folder,)]

        with tqdm(total=len(functions), desc="Creating plots", unit="Plot") as pbar:
            for func, args in zip(functions, args_list):
                func(*args)

    file_log.info(f"Annotation process completed. Results are available in {annotation_folder}")
    console_log.info(f"Annotation process completed. Results are available in {annotation_folder}")


def main():
    """
    Entry point for the VDJ Insights pipeline. Parses command-line arguments and executes the appropriate command.
    """
    console_log.info("Starting VDJ Insights pipeline.")
    parser = argparse.ArgumentParser(description="VDJ-Insights")
    subs = parser.add_subparsers(dest='command', required=True)

    setup_annotation_args(subs)
    setup_html(subs)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()

