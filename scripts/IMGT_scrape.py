import re
import time
import requests
from bs4 import BeautifulSoup
import argparse
from pathlib import Path
from urllib.parse import urlencode
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from logger import custom_logger

logger = custom_logger(__name__)
log_info = []  # List to collect log information
fasta_files_info = []  # List to collect fasta file information


def make_dir(dir):
    """
    Create a directory when not existing.

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
    log_info.append(
        {"message": f"Deleted folder: {directory.name}, because --cleanup was selected", "level": "warning"})
    logger.info(
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
    log_info.append(
        {"message": "Library created from generated files.", "level": "success"})
    logger.info("Creating a library from generated files.")


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
            log_info.append(
                {"message": "Succeeded to retrieve sequences from IMGT.", "level": "success"})
            logger.info("Succeeded to retrieve sequences from IMGT")
            return seq


def file_empty(file_path):
    """
    Checks if a file is empty by evaluating its size using the pathlib module. 
    This function assesses whether the file at the specified path exists and, if it does, 
    whether its size is zero bytes. It uses the `Path` object from `pathlib` to interact with the 
    file system, making it more intuitive and readable compared to older os-based methods.

    Args:
    file_path (str or Path): The path to the file to check. Can be a string or a Path object.

    Returns:
    bool: True if the file exists and is empty (size is 0 bytes), False otherwise.
    """
    path = Path(file_path)
    return path.exists() and path.stat().st_size == 0


def set_release():
    response = requests.get("https://www.imgt.org/vquest/refseqh.html")
    soup = BeautifulSoup(response.text, 'html.parser')
    release = soup.find(
        'h2', string=lambda text: text and 'IMGT/V-QUEST reference directory' in text)
    return str(re.findall(r'\(.*?\)', release.text.strip())[0][1:-1])


def write_sequence(name, directory, sequence):
    """
    Write the found VDJ segment sequences. It first establishes a base fasta file. 
    Then it writes the VDJ sequences to the fasta file and logs that 
    sequence are written to the file.

    Args:
        name (str): Name of current VDJ segment.
        directory (Path): Path to the directory.
        sequence (str): Large string containing all the different sequences
        for a current VDJ segment.
    """
    file = directory / f"{name}.fasta"
    log_info.append(
        {"message": f"Writing sequences from {name} to {file.name}.", "level": "success"})
    logger.info(f"Writing sequences from {name} to {file.name}")
    with open(file, 'w') as f:
        f.write(str(sequence) + "\n")
    # Count the number of entries in the fasta file
    count = sum(1 for _ in SeqIO.parse(file, "fasta"))
    fasta_files_info.append({"name": file.name, "entries": count})
    return file_empty(file)

    # Count the number of entries in the fasta file
    count = sum(1 for _ in SeqIO.parse(file, "fasta"))
    fasta_files_info.append({"name": file.name, "entries": count})


def fetch_sequence(segment, directory, species, frame, retry_limit=3):
    """
    Fetches sequence data from IMGT server using a constructed URL and handles failures with specified retry logic.
    Implements differentiated wait times based on the cause of retry need.

    Args:
        segment (str): The segment type.
        directory (Path): Directory path where sequence data is to be stored.
        species (str): Species name.
        frame (str): Frame specification.
        retry_limit (int): Maximum number of retries for fetching data.
    """
    for attempt in range(retry_limit):
        url = construct_url(segment, species, frame)
        response = requests.get(url)
        sleep_time = 2
        if response.status_code == 200:
            sequence = scrape(response)
            if sequence:
                empty_value = write_sequence(segment, directory, sequence)
                if not empty_value:
                    logger.info("Sequence successfully retrieved and written.")
                    break
                else:
                    sleep_time = 30
                    logger.warning(
                        "Retrieved sequence was empty. Will retry after extended wait.")
            else:
                logger.warning(
                    f"No sequences found for {segment} of {species}.")
        else:
            logger.warning(
                f"Failed to fetch data for {segment} of {species} with status code: {response.status_code}")
        logger.info(
            f"Waiting {sleep_time} seconds to avoid overloading the IMGT server.")
        time.sleep(sleep_time)
    if attempt == retry_limit - 1 and response.status_code == 200 and not sequence:
        logger.error(
            f"Max retries reached for {segment} of {species} without successful data retrieval.")


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
        log_info.append(
            {"message": f"Folder created: {directory.name}.", "level": "success"})
        logger.info(f"Folder {directory.name} does not exist, creating it!")
        make_dir(directory)
    for segment in segments[immune_type]:
        segment_file = directory / f"{segment}.fasta"
        log_info.append(
            {"message": f"Retrieving sequences from IMGT for the {segment} of {species}.", "level": "info"})
        logger.info(
            f"Retrieving sequences from IMGT for the {segment} of {species}")
        if not Path(segment_file).exists():
            fetch_sequence(segment, directory, species, frame)
        else:
            # Count the number of entries in the fasta file
            count = sum(1 for _ in SeqIO.parse(segment_file, "fasta"))
            fasta_files_info.append(
                {"name": segment_file.name, "entries": count})
            log_info.append(
                {"message": f"File exists: {segment}.fasta already exists, skipping!", "level": "info"})
            logger.info(f"File {segment}.fasta already exists skipping!")
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


def save_log_to_html(release, fasta_files_info, file_path, args):
    """
    Save the collected log information to an HTML file using Jinja2 template.

    Args:
        log_info (list): List of log strings to save.
        fasta_files_info (list): List of dictionaries containing fasta file names and number of entries.
        file_path (str): Path to the HTML file.
        args (argparse.Namespace): Parsed command line arguments.
    """
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template('source/template/template.html')

    date_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_content = template.render(
        date_time=date_time,
        set_release=release,
        species=args.species,
        type=args.type,
        output=args.output,
        frame_selection=args.frame_selection,
        create_library=args.create_library,
        cleanup=args.cleanup,
        simple_headers=args.simple_headers,
        fasta_files=fasta_files_info
    )

    with open(file_path, 'w') as html_file:
        html_file.write(html_content)


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
    global args  # Make args global to access it in save_log_to_html function
    args = argparser_setup()
    log_info.append(
        {"message": f"Starting scrape for species: {args.species}, type: {args.type}.", "level": "info"})
    logger.info(
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

    log_info.append(
        {"message": "Scrape completed successfully.", "level": "success"})
    logger.info("Scrape completed successfully.")

    # Save log to HTML
    save_log_to_html(set_release(), fasta_files_info,
                     "source/html/scrape_report.html", args)


if __name__ == '__main__':
    main()
