from mapping import mapping_main
from RSS import RSS_main
from report import report_main
from pathlib import Path
import pandas as pd
import argparse
from logger import custom_logger
from blast import blast_main
from map_genes import map_main
from extract_region import region_main

"""
Used python packages:
    1. pandas
    2. openpyxl
"""
# Method for logging current states of the program.
logger = custom_logger(__name__)


def make_dir(dir):
    """
    Create a directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """
    Path(dir).mkdir(parents=True, exist_ok=True)
    return dir


def combine_df(mapping_tools, cell_type, input_dir, library, threads):
    """
    Calling the mapping script annotation.py. This script is called
    with either "bowtie", "bowtie2" or "minimap2" as input value.
    These all return a df.
    Then all df's are concatenated to form a new complete one.
    Finally dropping all duplicates in the df based on a
    subset of ["start", "stop"], resetting the index of this df
    and return this new df.

    Returns:
        unique_combinations (DataFrame): A df containing all the unique
        mapping entries from bowtie(2) and minimap2.
    """
    df = pd.DataFrame()
    for tool in mapping_tools:
        mapping_df = mapping_main(tool, cell_type, input_dir, library, threads)
        df = pd.concat([df, mapping_df])
        df["haplotype"] = df["file"].str.extract(r'_([^_]+)\.')[0]
    unique_combinations = df.drop_duplicates(
        subset=["start", "stop", "haplotype"])
    return unique_combinations.reset_index(drop=True)


def write_report(df, report):
    """
    Creates an excel (xlsx) file of the df and saves it as "report.xlsx".

    Args:
        df (DataFrame): df to be saved.
        report (Path): Path of an excel to be saved to.
    """
    df.to_excel(report, index=False)


def get_or_create(cell_type, annotation_folder, mapping_tool, input_dir, library, threads):
    """
    Verifies if the report.xlsx is present.
    If present it returns the content of the file as a df. Otherwise
    call the combine_df() function to create the df;
    then save it as report.xlsx and return the df.

    Args:
        annotation_folder (Path): Path of the annotation folder.

    Returns:
        df (DataFrame): df containing information from report.xlsx.
    """
    report = annotation_folder / "report.xlsx"
    if not report.exists():
        logger.info("The report.xlsx file does not exist! Creating it!")
        df = combine_df(mapping_tool, cell_type, input_dir, library, threads)
        write_report(df, report)
        return df
    else:
        return pd.read_excel(report)


def validate_directory(directory_path):
    """Check if the specified directory exists."""
    if not Path(directory_path).is_dir():
        raise argparse.ArgumentTypeError(
            f"The directory {directory_path} does not exist. \
                Try another directory!")
    return directory_path


def validate_file(file_path):
    """Check if the specified file exists."""
    if not Path(file_path).is_file():
        raise argparse.ArgumentTypeError(
            f"The file {file_path} does not exist. \
                Try another file please!")
    return file_path


def validate_input(input_path):
    input_path = Path(input_path)
    validate_directory(input_path)
    if not any(entry.is_file() for ext in ["*.fasta", "*.fa", "*.fna"] for entry in input_path.glob(ext)):
        raise argparse.ArgumentTypeError(
            f"The directory {input_path} is empty or does not contain any fasta files!")
    return input_path


def validate_flanking_genes(value):
    """
    Validates and processes a list of flanking genes to ensure they are uppercase and handle empty values.

    Args:
        Comma-separated string of flanking genes.
    Returns:
        List of processed flanking genes.
    """
    flanking_genes = [gene.strip().upper() if gene.strip() !=
                      '-' else '' for gene in value.split(',')]
    if len(flanking_genes) % 2 == 1:
        raise argparse.ArgumentTypeError(
            f"The specified flanking genes: {flanking_genes} should be even numbers (e.g., 2, 4, 6, 8) rather than odd (e.g., 1, 3, 5).")
    return flanking_genes


def argparser_setup(include_help=True):
    """
    Set up the argument parser for the annotation tool.

    Args:
        include_help (bool): Whether to include the help argument. Defaults to True.
    Returns:
        ArgumentParser: Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description='A tool for finding known and novel VDJ segments in certain data, with a library of choice.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=include_help  # Control whether to add the help argument
    )

    required_group = parser.add_argument_group('Required Options')
    required_group.add_argument('-l', '--library', required=True, type=validate_file,
                                help='Path to the library file. Expected to be in FASTA format.')
    required_group.add_argument(
        '-r', '--receptor-type', required=True, type=str.upper, choices=['TR', 'IG'],
        help='Type of receptor to analyze: TR (T-cell receptor) or IG (Immunoglobulin).')

    regions_or_assembly_group = parser.add_argument_group('Data Source Options',
                                                          'Select the data source: regions or assembly. These options are mutually exclusive.')
    data_choice = regions_or_assembly_group.add_mutually_exclusive_group(
        required=True)
    data_choice.add_argument('-i', '--input', type=validate_input,
                             help='Directory containing the extracted sequence regions in FASTA format, where VDJ segments can be found. Cannot be used with -f/--flanking-genes or -s/--species.')
    data_choice.add_argument('-a', '--assembly', type=validate_input,
                             help='Directory containing the assembly FASTA files. Must be used with -f/--flanking-genes and -s/--species.')

    assembly_options = parser.add_argument_group('Assembly-Specific Options',
                                                 'These options are required if -a/--assembly is chosen:')
    assembly_options.add_argument('-f', '--flanking-genes', type=validate_flanking_genes,
                                  help='Comma-separated list of flanking genes, e.g., MGAM2,EPHB6. Add them as pairs. Required with -a/--assembly.')
    assembly_options.add_argument('-s', '--species', type=str,
                                  help='Species name, e.g., Homo sapiens. Required with -a/--assembly.')

    optional_group = parser.add_argument_group('Optional Options')
    optional_group.add_argument('-o', '--output', type=str,
                                default='annotation',
                                help='Output directory for the results.')
    mapping_options = ['minimap2', 'bowtie', 'bowtie2']
    optional_group.add_argument('-m', '--mapping-tool', nargs='*',
                                choices=mapping_options, default=mapping_options,
                                help='Mapping tool(s) to use. Choose from: minimap2, bowtie, bowtie2. Defaults to all.')

    optional_group.add_argument('-t', '--threads', type=int,
                                required=False, default=8, help='Amount of threads to run the analysis.')
    return parser


def main(args=None):
    """
    Main function for this annotation.py script. 
    It first sets a cwd (current directory the user is in) object 
    based on the current location. Based on this cwd some 
    input and output directories and files are set. It fetches the 
    needed blast db and the DataFrame (df). 
    Then it checks if the "blast_results.xlsx" is created. If this is not
    the case the needed df is made. Then the query cov is converted to numeric
    and checked if it is equal to 100%. The start and stop coordinates 
    are filtered out of the query column and also stored in the df.
    This df is saved to "blast_results.xlsx".
    The write_annotation_report() function is called.  
    """
    update_args = argparser_setup()
    region_dir = "region"
    if args is None:
        args = update_args.parse_args()
    elif isinstance(args, list):
        args = update_args.parse_args(args)
    elif not isinstance(args, argparse.Namespace):
        raise ValueError("Invalid arguments passed to the main function")
    if args.assembly:
        if not args.flanking_genes or not args.species:
            update_args.error(
                '-a/--assembly requires -f/--flanking-genes and -s/--species.')
        args.species = args.species.capitalize() if args.species else None

    if args.input and (args.flanking_genes or args.species):
        update_args.error(
            '-i/--input cannot be used with -f/--flanking-genes or -s/--species.')
    if args.input:
        region_dir = args.input
    cwd = Path.cwd()
    annotation_folder = cwd / args.output
    make_dir(annotation_folder)
    if args.assembly:
        map_main(args.flanking_genes, args.assembly, args.species)
        region_main(args.flanking_genes, args.assembly)
    df = get_or_create(args.receptor_type, annotation_folder, args.mapping_tool,
                       region_dir, args.library, args.threads)
    blast_file = annotation_folder / "blast_results.xlsx"
    if not blast_file.exists():
        blast_main(df, blast_file, args.library)
    report_main(annotation_folder, blast_file, args.receptor_type, args.library)
    RSS_main()
    logger.info(
        f"Annotation process completed. Results are available in {args.output}.")


if __name__ == '__main__':
    main()
