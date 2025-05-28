import argparse
import os
import shutil
import sys
import threading
import webbrowser
import subprocess
from pathlib import Path

import yaml

from .logger import console_logger, file_logger
from .annotation import main as annotation_main
from .IMGT_scrape import main as imgt_main
from .env_manager import create_and_activate_env, deactivate_env
from .util import make_dir, validate_file, validate_directory, validate_input, load_config

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
        shutil.copytree(
            settings_dir / "flask",
            output_dir,
            dirs_exist_ok=True
        )


def validate_flanking_genes(value):
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
    p.set_defaults(func=run_annotation)


def setup_html(subparsers):
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
    p = Path(value).resolve()
    if p.name == 'flask':
        validate_directory(str(p))
    else:
        validate_directory(str(p / 'flask'))
        p = p / 'flask'
    return p


def run_annotation(args):
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
    webbrowser.open_new('http://127.0.0.1:5000/')


def get_python_executable():
    for cmd in ("python", "python3"):
        if shutil.which(cmd):
            return cmd
    raise RuntimeError("No python executable found in PATH.")


def run_html(args):
    console_log.info("Launching HTML report...")
    output_dir = args.input
    copy_flask(output_dir, args.reset_flask)
    os.chdir(output_dir.parent)
    if not args.dev_mode:
        threading.Timer(1, open_browser).start()
    cmd = [get_python_executable(), str(output_dir / 'app.py')]
    subprocess.Popen(cmd).communicate()


def log_subprocess_error(e):
    file_log.error(f"Subprocess failed: {e}")
    file_log.error(e.stdout.decode())


def deep_merge(d1, d2):
    for k, v in d2.items():
        if k in d1 and isinstance(d1[k], dict) and isinstance(v, dict):
            deep_merge(d1[k], v)
        else:
            d1[k] = v


def create_config(output_dir, settings_dir, args):
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



def main():
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
