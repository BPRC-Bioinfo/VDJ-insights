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
import webbrowser
from pathlib import Path
import json
import time
from datetime import datetime

import pandas as pd
import yaml
from tqdm import tqdm

from .CDR import main_cdr
from .IMGT_scrape import main as imgt_main, set_release
from .RSS import main_rss
from .blast import blast_main

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


def setup_annotation_args2(subparsers):
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
    p.add_argument(
        '--verbose',
        help='Enable verbose logging.',
        action='store_true'
    )
    p.set_defaults(func=run_annotation, parser=p)



def setup_annotation_args(subparsers):
    """
    Configures the command-line arguments for the annotation command.
    """

    class IndentedHelpFormatter(argparse.HelpFormatter):
        def __init__(self, prog):
            super().__init__(prog,
                             indent_increment=2,
                             max_help_position=50,
                             width=None)

    p = subparsers.add_parser(
        'annotation',
        usage=(
            "vdj_insights.py annotation [OPTIONS]\n\n"
        ),
        formatter_class=IndentedHelpFormatter,
        add_help=False,
    )

    p.add_argument(
        '-h', '--help', action='help',
        help='Show this help message and exit\n'
    )

    group = p.add_mutually_exclusive_group()
    group.add_argument(
        '-f', '--flanking-genes', metavar='<str>',
        type=validate_flanking_genes,
        help='Comma-separated list of flanking genes (must be even count).'
    )

    # Required arguments
    req = p.add_argument_group('required arguments')
    req.add_argument(
        '-r', '--receptor-type', metavar='<TR|IG>',
        required=True, type=str.upper, choices=['TR', 'IG'],
        help='Receptor type: TR (T-cell receptor) or IG (immunoglobulin).'
    )
    data_choice = req.add_mutually_exclusive_group(required=True)
    data_choice.add_argument(
        '-i', '--input', metavar='<dir>',
        type=validate_input,
        help='Directory with input FASTA regions.'
    )
    data_choice.add_argument(
        '-a', '--assembly', metavar='<dir>',
        type=validate_input,
        help='Directory with genome assembly FASTA (requires --flanking-genes and --species).'
    )

    opt = p.add_argument_group('optional arguments')
    opt.add_argument(
        '-l', '--library', metavar='<file>',
        type=validate_file,
        help='Path to library file.'
    )
    opt.add_argument(
        '-S', '--scaffolding', metavar='<file>',
        type=validate_file,
        help='Path to reference genome file.'
    )
    opt.add_argument(
        '-M', '--metadata', metavar='<file>',
        type=validate_file,
        help='Path to metadata file.'
    )
    opt.add_argument(
        '-s', '--species', metavar='<str>',
        type=str,
        help='Species name (required with --assembly).'
    )
    opt.add_argument(
        '-o', '--output', metavar='<path>',
        type=str,
        default=str(Path.cwd() / 'annotation_results'),
        help='Path to output directory. [default: <current working directory>/annotation_results]'
    )
    opt.add_argument(
        '-m', '--mapping-tool', metavar='<tool>',
        nargs='*', choices=['minimap2', 'bowtie', 'bowtie2'],
        default=['minimap2', 'bowtie', 'bowtie2'],
        help='Mapping tools to use. [default: minimap2, bowtie, bowtie2]'
    )
    opt.add_argument(
        '-t', '--threads', metavar='<int>',
        type=int, default=8,
        help='Number of threads. [default: 8]'
    )
    opt.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging.'
    )

    p.set_defaults(func=run_annotation, parser=p)


def setup_html(subparsers):
    """
    Configures the command-line arguments for the HTML report viewer.

    Args:
        subparsers (argparse._SubParsersAction): The subparsers object to which the 'html' command will be added.
    """
    p = subparsers.add_parser(
        'html',
        help='Serve and view the annotation report in a browser.'
    )
    p.add_argument(
        '-i', '--input',
        required=True,
        type=validate_html,
        help='Path to the directory containing annotation data.'
    )
    p.add_argument(
        '--port',
        default=5002,
        type=int,
        help='Port to run the local server on (default: 5002).'
    )
    p.add_argument(
        '--host',
        default="0.0.0.0",
        type=str,
        help='Host address to bind the server (default: 0.0.0.0).'
    )
    p.add_argument(
        '-d', '--dev-mode',
        action='store_true',
        help='Run in development mode and display Flask logs in the console.'
    )

    p.set_defaults(func=run_html)



def setup_scrape_args(subparsers):
    """
    Configures the command-line arguments for the scrape command.

    Args:
        subparsers (argparse._SubParsersAction): Subparsers object for adding commands.
    """
    p = subparsers.add_parser(
        'scrape',
        help='Scrape immune gene segments from IMGT.'
    )
    p.add_argument(
        '-r', '--receptor-type',
        required=True,
        type=str.upper,
        choices=['TR', 'IG'],
        help='TR or IG.'
    )
    p.add_argument(
        '-s', '--species',
        required=True,
        type=str,
        help='Species name.'
    )
    p.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        help='Output directory for scraped data.'
    )
    p.set_defaults(func=run_scrape, parser=p)


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
    create_and_activate_env(settings_dir / 'envs' / 'vdj-insights_env.yaml', args.verbose)
    create_config(output_dir, settings_dir, args)
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
    os.chdir(output_dir.parent)
    if not args.dev_mode:
        threading.Timer(1, open_browser).start()
    cmd = [
        "flask", "--app",
        str(output_dir / 'app.py'),
        "run",
        "--host", str(args.host),
        "--port", str(args.port)
    ]
    if args.dev_mode:
        cmd.append("--debug")
    subprocess.Popen(cmd).communicate()


def run_scrape(args):
    """
    Runs the IMGT scraping tool.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    console_log.info(f"Starting scrape for species: {args.species} ({args.receptor_type})")
    imgt_main(species=args.species, immune_type=args.receptor_type, output_dir=args.output)


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


def save_args(args):
    cwd = Path.cwd()
    args_path = cwd / "used_commando.json"

    skip_keys = {"func", "parser"}
    def is_json_serializable(v):
        return isinstance(v, (str, int, float, bool, list, dict, type(None)))

    serializable_args = {
        k: (v if is_json_serializable(v) else str(v))
        for k, v in vars(args).items()
        if k not in skip_keys
    }

    serializable_args["command line"] = " ".join(sys.argv)
    with open(args_path, "w") as f:
        json.dump(serializable_args, f, indent=4)


def initialize_config(args, species_cfg, settings_dir, config):
    parser = args.parser

    if getattr(args, 'assembly', None) and not args.species:
        parser.error("--assembly requires --species/-s to be set")

    if getattr(args, 'flanking_genes', None):
        flanking_genes = args.flanking_genes
    else:
        species_key = args.species.capitalize().replace(" ", "_")

        species_data = species_cfg.get(species_key, species_cfg['default'])
        receptor_data = species_data.get(args.receptor_type, species_cfg['default'][args.receptor_type])
        flanking_genes = receptor_data.get("FLANKING_GENES")

        if species_key not in species_cfg:
            console_log.warning(f"Species '{args.species}' not found in configuration. falling back to default flanking genes. {flanking_genes}")

    args.flanking_genes = flanking_genes


@log_error()
def annotation_main(args: argparse.Namespace):
    """
    Main function for the annotation tool. Performs the following steps:

    1. Parses command-line arguments and sets up paths for input and output files.
    2. Validates input parameters and configures the environment based on the selected options.
    3. Creates necessary directories and initializes the annotation process.
    4. Executes the mapping and region extraction processes if assembly mode is selected:
      - Calls map_main to map flanking genes against the assembly.
      - Calls region_main to extract regions based on mapped genes.
    5. Retrieves or creates a DataFrame from existing or new mapping results:
      - Uses get_or_create to check for existing reports or generate new ones.
      - Calls combine_df to merge mapping results and remove duplicates.
    6. Runs BLAST operations to align sequences and generate results:
      - Calls blast_main to perform BLAST alignment on the combined DataFrame.
    7. Generates the final annotation report and RSS file:
      - Calls report_main to generate an Excel report summarizing the findings.
      - Calls RSS_main to produce an RSS feed if required.

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

    timings = {}
    start_total = time.time()

    # Scaffold
    if args.scaffolding:
        t0 = time.time()
        scaffolding_dir = "tmp/scaffold_assemblies"
        args.assembly = scaffolding_main(args.scaffolding, args.assembly, scaffolding_dir, args.threads, args.verbose)
        timings["scaffolding"] = round(time.time() - t0, 2)

    # Library
    console_log.info(f"IMGT scrape: {args.species} ({args.receptor_type})")
    lib_dest = cwd / 'library' / 'library.fasta'
    if not args.library:
        if not lib_dest.is_file():
            t0 = time.time()
            imgt_main(species=args.species, immune_type=args.receptor_type)
            args.library = lib_dest
            timings["imgt_scraping"] = round(time.time() - t0, 2)
        else:
            json_dest = cwd / 'library' / "library_info.json"
            with open(json_dest, "r") as f:
                config = json.load(f)

            old_release = config["set_release"]
            new_release = set_release()
            if old_release != new_release:
                console_log.warning(f"New release detected of library: {old_release} → {new_release}")
            args.library = lib_dest
    else:
        make_dir(cwd / 'library')
        shutil.copy(args.library, lib_dest)
        args.library = lib_dest

    annotation_folder = cwd / 'annotation'
    make_dir(annotation_folder)

    # Mapping en regio-extractie
    flanking_genes_dict = ast.literal_eval(str(args.flanking_genes))
    if args.assembly:
        t0 = time.time()
        map_main(flanking_genes_dict, args.assembly, args.species, args.threads, args.verbose)
        timings["mapping flanking genes"] = round(time.time() - t0, 2)

        t0 = time.time()
        region_main(flanking_genes_dict, args.assembly, args.threads, args.verbose)
        timings["extract regions"] = round(time.time() - t0, 2)

    # Mapping library
    t0 = time.time()
    region_dir = args.input if args.input else "tmp/region"
    report = annotation_folder / "tmp/report.csv"
    if not report.exists():
        report_df = pd.DataFrame()
        for tool in args.mapping_tool:
            file_log.info(f"Processing tool: {tool}")
            mapping_df = mapping_main(tool, region_dir, args.library, args.threads, args.verbose)
            report_df = pd.concat([report_df, mapping_df])
        report_df = report_df.drop_duplicates(subset=["reference", "start", "stop", "name"]).reset_index(drop=True)
        make_dir(annotation_folder / "tmp/")
        report_df.to_csv(report, index=False)
    else:
        report_df = pd.read_csv(report)
    timings["library_mapping"] = round(time.time() - t0, 2)

    # BLAST
    t0 = time.time()
    blast_file = annotation_folder / "tmp/blast_results.csv"
    if not blast_file.exists():
        blast_main(report_df, blast_file, args.library, args.threads, args.verbose)
    timings["blast"] = round(time.time() - t0, 2)

    # Report + RSS + CDR
    t0 = time.time()
    report_main(annotation_folder, blast_file, args.receptor_type, args.library, args.assembly, args.metadata)
    timings["report"] = round(time.time() - t0, 2)

    t0 = time.time()
    main_functionality(args.receptor_type, args.species, args.threads, args.verbose)
    timings["functionality"] = round(time.time() - t0, 2)

    t0 = time.time()
    main_rss(args.threads, args.verbose)
    timings["rss"] = round(time.time() - t0, 2)

    t0 = time.time()
    main_cdr(args.species, args.receptor_type, args.threads, args.verbose)
    timings["cdr"] = round(time.time() - t0, 2)

    # BED + GTF
    data = pd.read_excel(annotation_folder / "annotation_report_all.xlsx")
    make_bed(data, annotation_folder / "BED")
    make_gtf(data, annotation_folder / "GTF")
    save_args(args)

    # Figures
    t0 = time.time()
    if args.metadata:
        functions = [barplot_main, boxplot_main, sub_families_main, venn_diagram_main, heatmap_main]
        args_list = [(annotation_folder,), (annotation_folder,), (annotation_folder,), (annotation_folder, args.receptor_type), (annotation_folder,)]
        with tqdm(total=len(functions), desc="Creating plots", unit="Plot") as pbar:
            for func, args_ in zip(functions, args_list):
                func(*args_)
                pbar.update()
    timings["figures"] = round(time.time() - t0, 2)

    timings["total_time"] = round(time.time() - start_total, 2)
    timings["date"] = datetime.today().strftime("%Y-%m-%d")

    with open("timing.json", "w") as f:
        json.dump(timings, f, indent=4)


    file_log.info(f"Annotation process completed. Results are available in {annotation_folder}")
    console_log.info(f"Annotation process completed. Results are available in {annotation_folder}")
    console_log.info(f"Showing report:  vdj-insights html - i {cwd}")


def main():
    """
    Entry point for the VDJ Insights pipeline. Parses command-line arguments and executes the appropriate command.
    """
    parser = argparse.ArgumentParser(description="VDJ-Insights")
    subs = parser.add_subparsers(dest='command', required=True)

    setup_annotation_args(subs)
    setup_html(subs)
    setup_scrape_args(subs)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    args.func(args)
    console_log.info("Starting VDJ Insights pipeline.")



if __name__ == '__main__':
    main()

