import time
import requests
from bs4 import BeautifulSoup
import logging
import argparse
from pathlib import Path
from urllib.parse import urlencode
from Bio import SeqIO

# Custom formatter for the logging module

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


def make_dir(dir):
    """
    Create an directory when not existing.

    Args:
        location (str): Path of the directory to create.
    """

    Path(dir).mkdir(parents=True, exist_ok=True)


def cleanup(directory):
    """
    Gets a path of a directory as input and loops over the fasta files 
    inside. It removes every found fasta file. After the removal of the 
    fasta files it will remove the directory itself.

    Args:
        directory (Path): Path to the directory.
    """
    for file in directory.glob("*.fasta"):
        Path(file).unlink()
    Path(directory).rmdir()
    logging.info(
        f"Deleting folder: {directory.name}, because --cleanup was selected")


def edit_header(header):
    """
    Simplifies the header of the current record of that is being parsed.
    It keeps only the ID, region, segment and species.

    Args:
        header (str): Old header of the record.

    Returns:
        str: Parsed header of the record.
    """
    new_header = "_".join(header.split("|")[0:3])
    return new_header.replace(" ", "_")


def create_library(directory: Path, simple_headers):
    """
    Gets a path of a directory that contain the fasta files needed to 
    create the library.fasta file. It first creates the library directory.
    Then loops over the fasta files and opens it with biopython SeqIO and 
    append the contents of the current fasta file to the library.fasta file.
    In the end there is logged that the library is created.


    Args:
        directory (Path): Path to the directory.
    """
    library = directory.parent / "library"
    make_dir(library)
    with open(library / "library.fasta", 'w') as w:
        for file in directory.glob("*.fasta"):
            for record in SeqIO.parse(file, "fasta"):
                if simple_headers:
                    record.description = edit_header(record.description)
                w.write(f">{record.description}\n{record.seq.upper()}\n")
    logging.info("Creating a library from generated files.")


def construct_url(segment, species, frame):
    """
    Constructs the url that is needed to fetch the segments from the
    IMGT web server based on the segment type, species and frame. 
    It takes the base url and adds the right values for segment, 
    species and frame, which are stored in the params dict. 
    With urlencode the dict is converted to a url.

    Args:
        segment (str): The segment type itself.
        species (str): The chosen species, which is choses from the
        argparse list for species.
        frame (str): the chosen frame, which is choses from the argparse
        list for regarding frame.

    Returns:
        str: constructed url for data retrieval from the IMGT web server.
    """
    base_url = "https://www.imgt.org/genedb/GENElect"
    params = {
        'query': f"{frame} {segment}",
        'species': species
    }
    return f"{base_url}?{urlencode(params)}"


def scrape(response):
    """
    Uses the response and parses it with BeautifulSoup. 
    It first create a BeautifulSoup object from the response and searches
    for all the "pre" paragraphs in the "soup" object. Then it checks if ">" 
    is present in the paragraph, if this is the case the paragraph gets
    returned and gets logged that VDJ sequences are found.

    Args:
        response (requests.models.Response): Response containing
        information of the IMGT html page.

    Returns:
        str: All the found VDJ sequences in this response/paragraph. 
    """
    soup = BeautifulSoup(response.text, 'html.parser')
    paragraphs = soup.find_all('pre')
    for p in paragraphs:
        seq = p.text.strip()
        if ">" in seq:
            logging.info(f"Succeeded to retrieve sequences from IMGT")
            return seq


def write_sequence(name, directory, sequence):
    """
    Write the found VDJ segment sequenecs. It first establish a base fasta file. 
    Then it writes the VDJ sequences to the fasta file and logs that 
    sequence are written to the file.

    Args:
        name (str): Name of current VDJ segment.
        directory (Path): Path to the directory.
        sequence (str): Large string containing all the different sequences
        for a current VDJ segment.
    """
    file = directory / f"{name}.fasta"
    logging.info(f"Writing sequences from {name} to {file.name}")
    with open(file, 'w') as f:
        f.write(str(sequence) + "\n")


def fetch_sequence(segment, directory, species, frame):
    """
    Uses the constrcted url and fetches the response based on the url. 
    It checks the response code. If it is equal 200 the response is
    scraped otherwise it is logged that it failed to retrieve data from IMGT server.
    it also checks if something is found after scraping. 
    Otherwise it logs that it could not scrape any data from the response.
    In the end a 2 second sleep it called to avoid overloading the IMGT server.


    Args:
        segment (str): The segment type itself.
        directory (Path): Path to the directory.
        species (str): The chosen species, which is choses from the
        argparse list for species.
        frame (str): the chosen frame, which is choses from the argparse
        list for regarding frame.
    """
    url = construct_url(segment, species, frame)
    response = requests.get(url)
    if response.status_code == 200:
        sequence = scrape(response)
        if sequence:
            write_sequence(segment, directory, sequence)
        else:
            logging.warning(f"No sequences found for {segment} of {species}.")
    else:
        logging.warning(
            f"Failed to fetch data for {segment} of {species}. Status code: {response.status_code}")
    logging.info("Waiting 2 seconds to avoid overloading the IMGT server!")
    time.sleep(2)


def scrape_IMGT(species, immune_type, directory, frame):
    """
    Uses the argparse values to scrape from the IMGT server. First a dictionary 
    is created with the different VDJ segments that are known.
    First it checks if the output folder exists, if not it creates 
    it otherwise it uses it. Then loop over the right VDJ segments list,
    chosen from the dict based on the immune_type. If fasta file exists based 
    on the current segment it is not fetched, logged that it already exists 
    and a 2 second sleep is called. 
    Otherwise the segment sequences are fetched and logged.

    Args:

        species (str): The chosen species, which is choses from the
        argparse list for species.
        immune_type (str): The type of receptor that is being used. 
        Either TR or IG.
        directory (Path): Path to the directory.
        frame (str): the chosen frame, which is choses from the argparse
        list for regarding frame.
    """
    segments = {
        "TR": [
            "TRBV", "TRBJ", "TRBD", "TRAV", "TRAJ",
            "TRDD", "TRDJ", "TRDV", "TRGV", "TRGJ"
        ],
        "IG": [
            'IGHV', 'IGHD', 'IGHJ', 'IGKV',
            'IGKJ', 'IGLV', 'IGLJ'
        ]
    }
    if not directory.exists():
        logging.info(f"Folder {directory.name} does not exist, creating it!")
        make_dir(directory)
    for segment in segments[immune_type]:
        logging.info(
            f"Retrieving sequences from IMGT for the {segment} of {species}")
        if not Path(directory / f"{segment}.fasta").exists():
            fetch_sequence(segment, directory, species, frame)
        else:
            logging.info(f"File {segment}.fasta already exists skipping!")
            time.sleep(2)


def to_capital(species):
    """
    Capitalize the input given for the species argument with argparse.

    Args:
        species (str): The chosen species, which is choses from the
        argparse list for species.

    Returns:
        str: Capitalized species string.
    """
    return species.capitalize()


def convert_frame(frame):
    """
    Transform the given frame to its secondary value. It creates a dict 
    for this and converts the value based on this dict. 
    If the frame has the None type, it is set to all.

    Args:
        frame (str/None): The chosen frame. If not set its value is None.

    Returns:
        str: decimal number indicating the which type of frame is needed.
    """
    options = {
        "all": "7.2", "in-frame": "7.5", "in-frame-gaps": "7.1"
    }
    return options[frame or "all"]


def validate_directory(path):
    """
    Validate the given path. First the path is converted to a Path object
    and than there is checked if the path to a directory exists. If it does the
    path is returned, otherwise the directory is created and returned.

    Args:
        path (str): Path to a directory.

    Returns:
        Path: Path to a directory, converted to a Path object.
    """
    path = Path(path)
    if not path.exists():
        make_dir(path)
    return path


def argparser_setup():
    """
    Establish a list with all the different species the user can choose 
    from to fetch. Then argparse itself it initiated. A base description is set,
    as well as the two groups for required parameters and optional parameters.
    Then all the parameters are set and their extra function if validation 
    or transformation is needed. In the end a parser (args) is returned 
    with all the given parameters.

    Returns:
        args (argparse.Namespace): Object with all the given parameters.
    """
    latin_names = [
        'Homo sapiens', 'Mus', 'Bos taurus',
        'Ovis aries', 'Danio rerio', 'Canis lupus familiaris',
        'Macaca mulatta', 'Tursiops truncatus', 'Oryctolagus cuniculus',
        'Heterocephalus glaber', 'Macaca fascicularis',
    ]

    parser = argparse.ArgumentParser(
        description='Scrape IMGT for TR and IG segment sequences of a given species.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Grouping arguments
    required_group = parser.add_argument_group('Required Options')
    optional_group = parser.add_argument_group('Optional Options')

    # Species options
    required_group.add_argument('-S', "--species", type=to_capital, choices=latin_names, required=True,
                                help='Name of the species to scrape for (e.g., "Homo sapiens"). Capitalization is handled automatically.')
    required_group.add_argument('-T', '--type', type=str.upper, choices=['TR', 'IG'], required=True,
                                help='Type of sequence to scrape: TR (T-cell receptor) or IG (Immunoglobulin).')

    # Output options
    optional_group.add_argument('-O', '--output', type=validate_directory,
                                help='Output directory where the results will be saved. The directory will be created if it does not exist.')
    optional_group.add_argument('-f', '--frame-selection', type=str, choices=['all', 'in-frame', 'in-frame-gaps'],
                                help='ORF frame analysis type. Choices are "all" for F+ORF+all P, "in-frame" for F+ORF+in-frame P, or "in-frame-gaps" for F+ORF+in-frame P with IMGT gaps.')
    optional_group.add_argument('--create-library', action='store_true',
                                help='Create a library from the IMGT files if specified.')
    optional_group.add_argument('--cleanup', action='store_true',
                                help='Clean up leftover IMGT files after processing.')
    optional_group.add_argument('--simple-headers', action='store_true',
                                help='Create simplified headers to improve readability.')

    args = parser.parse_args()
    return args


def main():
    """
    Main function of the IMGT_scrape script. 
    It first retrieves all the given parameters from argparse. Then logs that 
    is starting the scrape for a given species and immune type. It determines 
    the output. If not given a Path object is created based on the species name.
    Otherwise a Path object is created based on the output parameter. 
    If the output directory is not created yet it is created. Next the scraping
    itself is done. After scraping there is checked if the library needs to be 
    created and/or that the individual fasta files of the segments
    need to be removed. Lastly there is logged that the scrape is finished.
    """
    args = argparser_setup()
    logging.info(
        f"Starting scrape for species: {args.species}, type: {args.type}")

    output_dir = Path(args.output) if args.output else Path.cwd(
    ) / args.species.replace(" ", "_").lower()
    frame_selection = convert_frame(args.frame_selection)
    make_dir(output_dir)
    scrape_IMGT(args.species, args.type, output_dir, frame_selection)
    if args.create_library:
        create_library(output_dir, args.simple_headers)

    if args.cleanup:
        cleanup(output_dir)

    logging.info("Scrape completed successfully.")


if __name__ == '__main__':
    main()
