import time
import requests
from bs4 import BeautifulSoup
import logging
import argparse
from pathlib import Path
from urllib.parse import urlencode
from Bio import SeqIO

COLORS = {
    'WARNING': '\033[93m',
    'INFO': '\033[92m',
    'DEBUG': '\033[94m',
    'CRITICAL': '\033[91m',
    'ERROR': '\033[91m',
    'ENDC': '\033[0m',
}


class CustomFormatter(logging.Formatter):
    def format(self, record):
        levelname = record.levelname
        message = logging.Formatter.format(self, record)
        return f"{COLORS.get(levelname, '')}{message}{COLORS['ENDC']}"


# Configure root logger
logging.basicConfig(level=logging.DEBUG)

# Get the root logger
logger = logging.getLogger()

# Create console handler with a higher log level
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)

# Create formatter and add it to the handler
formatter = CustomFormatter(
    '%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Remove all handlers associated with the root logger
logger.handlers = []

# Add the custom handler to the root logger
logger.addHandler(handler)


def make_dir(location):
    Path(location).mkdir(parents=True, exist_ok=True)


def cleanup(folder):
    for file in folder.glob("*.fasta"):
        Path(file).unlink()
    Path(folder).rmdir()
    logging.info(
        f"Deleting folder: {folder.name}, because --cleanup was selected")


def create_library(folder: Path):
    library = folder.parent / "library"
    make_dir(library)
    with open(library / "library.fasta", 'w') as w:
        for file in folder.glob("*.fasta"):
            for record in SeqIO.parse(file, "fasta"):
                w.write(f">{record.description}\n{record.seq.upper()}\n")


def construct_url(segment, species, frame):
    base_url = "https://www.imgt.org/genedb/GENElect"
    params = {
        'query': f"{frame} {segment}",
        'species': species
    }
    return f"{base_url}?{urlencode(params)}"


def scrape(response):
    soup = BeautifulSoup(response.text, 'html.parser')
    paragraphs = soup.find_all('pre')
    for p in paragraphs:
        seq = p.text.strip()
        if ">" in seq:
            logging.info(f"Succeeded to retrieve sequences from IMGT")
            return seq


def write_sequence(name, folder, sequence):
    file = folder / f"{name}.fasta"
    logging.info(f"Writing sequences from {name} to {file.name}")
    with open(file, 'w') as f:
        f.write(str(sequence) + "\n")


def fetch_sequence(segment, folder, species, frame):
    url = construct_url(segment, species, frame)
    response = requests.get(url)
    if response.status_code == 200:
        sequence = scrape(response)
        if sequence:
            write_sequence(segment, folder, sequence)
        else:
            logging.warning(f"No sequences found for {segment} of {species}.")
    else:
        logging.warning(
            f"Failed to fetch data for {segment} of {species}. Status code: {response.status_code}")
    logging.info("Waiting 2 seconds to avoid overloading the IMGT server!")
    time.sleep(2)


def scrape_IMGT(species, imune_type, folder, frame):
    segments = {"TCR": ["TRBV", "TRBJ", "TRBD", "TRAV", "TRAJ", "TRDD", "TRDJ", "TRDV", "TRGV", "TRGJ"],
                "IG": ['IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ']}
    if not folder.exists():
        logging.info(f"Folder {folder.name} does not exist, creating it!")
        make_dir(folder)
    for segment in segments[imune_type]:
        logging.info(
            f"Retrieving sequences from IMGT for the {segment} of {species}")
        if not Path(folder / f"{segment}.fasta").exists():
            fetch_sequence(segment, folder, species, frame)
        else:
            logging.info(f"File {segment}.fasta already exists skipping!")
            time.sleep(2)


def to_capital(string: str):
    return string.capitalize()


def convert_frame(frame):
    options = {"all": "7.2", "in-frame": "7.5", "in-frame-gaps": "7.1"}
    return options[frame or "all"]


def validate_directory(path):
    path = Path(path)
    if not path.exists():
        make_dir(path)
    return path


def argparser_setup():
    latin_names = [
        'Homo sapiens', 'Mus', 'Bos taurus',
        'Ovis aries', 'Danio rerio', 'Canis lupus familiaris',
        'Macaca mulatta', 'Tursiops truncatus', 'Oryctolagus cuniculus',
        'Heterocephalus glaber'
    ]

    parser = argparse.ArgumentParser(
        description='Scrape IMGT for TCR and IG segment sequences of a given species.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Grouping arguments
    required_group = parser.add_argument_group('Required Options')
    optional_group = parser.add_argument_group('Optional Options')

    # Species options
    required_group.add_argument('-S', "--species", type=to_capital, choices=latin_names, required=True,
                                help='Name of the species to scrape for (e.g., "Homo sapiens"). Capitalization is handled automatically.')
    required_group.add_argument('-T', '--type', type=str.upper, choices=['TCR', 'IG'], required=True,
                                help='Type of sequence to scrape: TCR (T-cell receptor) or IG (Immunoglobulin).')

    # Output options
    optional_group.add_argument('-O', '--output', type=validate_directory,
                                help='Output directory where the results will be saved. The directory will be created if it does not exist.')
    optional_group.add_argument('-f', '--frame-selection', type=str, choices=['all', 'in-frame', 'in-frame-gaps'],
                                help='ORF frame analysis type. Choices are "all" for F+ORF+all P, "in-frame" for F+ORF+in-frame P, or "in-frame-gaps" for F+ORF+in-frame P with IMGT gaps.')
    optional_group.add_argument('--create-library', action='store_true',
                                help='Create a library from the IMGT files if specified.')
    optional_group.add_argument('--cleanup', action='store_true',
                                help='Clean up leftover IMGT files after processing.')

    args = parser.parse_args()
    return args


def main():
    args = argparser_setup()
    logging.info(
        f"Starting scrape for species: {args.species}, type: {args.type}")

    output_dir = Path(args.output) if args.output else Path.cwd(
    ) / args.species.replace(" ", "_").lower()
    frame_selection = convert_frame(args.frame_selection)
    make_dir(output_dir)
    scrape_IMGT(args.species, args.type, output_dir, frame_selection)
    if args.create_library:
        logging.info("Creating a library from generated files.")
        create_library(output_dir)

    if args.cleanup:
        cleanup(output_dir)

    logging.info("Scrape completed successfully.")


if __name__ == '__main__':
    main()
