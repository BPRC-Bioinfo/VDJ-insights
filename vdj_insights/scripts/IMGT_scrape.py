"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

import json
import re
import time
import requests
from typing import Union

from bs4 import BeautifulSoup
import argparse
from pathlib import Path
from urllib.parse import urlencode
from Bio import SeqIO

from .util import make_dir
from .logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)

fasta_files_info = []
fasta_file_urls = []

def cleanup(directory: Union[str, Path]):
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
    file_log.info(f"Deleting folder: {directory.name}, because --cleanup was selected")


def create_library(directory: Union[str, Path], receptor: str, simple_headers):
    """
    Gets a path of a directory that contain the fasta files needed to 
    create the library.fasta file. It first creates the library directory.
    Then loops over the fasta files and opens it with biopython SeqIO,get
    the ID, region, segment and species and append the contents of the
    current fasta file to the library.fasta file. In the end there is
    logged that the library is created.

    Args:
        directory (Path): Path to the directory.
    """
    library = directory.parent / "library"
    make_dir(library)
    with open(library / f"library.fasta", 'w') as w:
        for file in directory.glob("*.fasta"):
            for record in SeqIO.parse(file, "fasta"):
                if simple_headers:
                    new_header = "_".join(record.description.split("|")[0:3])
                    record.description = new_header.replace(" ", "_")
                w.write(f">{record.description}\n{record.seq.upper()}\n")
    file_log.info("Creating a library from generated files.")
    return library


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
            file_log.info("Succeeded to retrieve sequences from IMGT")
            return seq



def set_release():
    response = requests.get("https://www.imgt.org/vquest/refseqh.html")
    soup = BeautifulSoup(response.text, 'html.parser')
    release = soup.find('h2', string=lambda text: text and 'IMGT/V-QUEST reference directory' in text)
    return str(re.findall(r'\(.*?\)', release.text.strip())[0][1:-1])


def write_sequence(name: str, directory: Union[str, Path], sequence: str):
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
    path = Path(directory / f"{name}.fasta")
    file_log.info(f"Writing sequences from {name} to {path.name}")

    with open(path, 'w') as f:
        f.write(str(sequence) + "\n")

    count = sum(1 for _ in SeqIO.parse(path, "fasta"))
    fasta_files_info.append({"name": path.name, "entries": count})

    return path.exists() and path.stat().st_size == 0



def fetch_sequence(segment: str, directory: Union[str, Path], species: str, frame: str, retry_limit=3):
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
        url = f"https://www.imgt.org/genedb/GENElect?{urlencode({'query': f'{frame} {segment}', 'species': species})}"
        response = requests.get(url)
        sleep_time = 2
        if response.status_code == 200:
            sequence = scrape(response)
            if sequence:
                empty_value = write_sequence(segment, directory, sequence)
                if not empty_value:
                    file_log.info("Sequence successfully retrieved and written.")
                    break
                else:
                    sleep_time = 30
                    file_log.warning("Retrieved sequence was empty. Will retry after extended wait.")
            else:
                file_log.warning(f"No sequences found for {segment} of {species}.")
        else:
            file_log.warning(f"Failed to fetch data for {segment} of {species} with status code: {response.status_code}")
        file_log.info(f"Waiting {sleep_time} seconds to avoid overloading the IMGT server.")
        time.sleep(sleep_time)
    if attempt == retry_limit - 1 and response.status_code == 200 and not sequence:
        file_log.error(f"Max retries reached for {segment} of {species} without successful data retrieval.")
        console_log.error(f"Max retries reached for {segment} of {species} without successful data retrieval.")
        raise
    fasta_file_urls.append(url)

def scrape_IMGT(species: str, immune_type: str, directory: Union[str, Path], frame: str):
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
    """
    segments = {
        "TR": [
            "TRBV", "TRBJ", "TRBD", "TRBC",
            "TRAV", "TRAJ", "TRAC",
            "TRDD", "TRDJ", "TRDC", "TRDV",
            "TRGV", "TRGJ"
        ],
        "IG": [
            "IGHV", "IGHD", "IGHJ", "IGHC",
            "IGKV", "IGKJ", "IGKC",
            "IGLV", "IGLJ", "IGLC"
        ]
    }
    """
    make_dir(directory)
    for segment in segments[immune_type]:
        segment_file = directory / f"{segment}.fasta"
        file_log.info(
            f"Retrieving sequences from IMGT for the {segment} of {species}")
        if not Path(segment_file).exists():
            fetch_sequence(segment, directory, species, frame)
        else:
            count = sum(1 for _ in SeqIO.parse(segment_file, "fasta"))
            fasta_files_info.append(
                {"name": segment_file.name, "entries": count})
            file_log.info(f"File {segment}.fasta already exists skipping!")
            time.sleep(2)


def convert_frame(frame: Union[str, None]):
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


def save_json(release: str,
              fasta_files_info: list[dict],
              fasta_files_url: list[str],
              file_path: Union[str, Path],
              species: str,
              immune_type: str,
              output_dir: Union[str, Path],
              frame_selection: str,
              create_library_bool: bool,
              cleanup_bool: bool,
              simple_headers_bool: bool):

    for info, url in zip(fasta_files_info, fasta_files_url):
        info['url'] = url

    json_dict = {
        "set_release": release,
        "species": species,
        "type": immune_type,
        "output": str(output_dir),
        "frame_selection": frame_selection,
        "create_library": create_library_bool,
        "cleanup": cleanup_bool,
        "simple_headers": simple_headers_bool,
        "fasta_files": fasta_files_info,
    }

    with open(file_path, 'w') as w:
        json.dump(json_dict, w, indent=4)


def main(species: str,
         immune_type: str,
         output_dir: Union[str, Path] = None,
         frame_selection: str = "all",
         create_library_bool: bool = True,
         cleanup_bool: bool = True,
         simple_headers_bool: bool = True
         ):
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
    file_log.info(f"Starting scrape for species: {species}, type: {immune_type}")

    cwd = Path.cwd()
    if output_dir:
        path = cwd / output_dir
    else:
        path = cwd
    frame_selection = convert_frame(frame_selection)
    imgt_dir = path / species.replace(" ", "_").lower()
    make_dir(imgt_dir)
    
    scrape_IMGT(species, immune_type, imgt_dir, frame_selection)
    if create_library_bool:
        library_dir = create_library(imgt_dir, immune_type, simple_headers_bool)

    if cleanup_bool:
        cleanup(imgt_dir)

    file_log.info("Scrape completed successfully.")
    json_file = library_dir / "library_info.json"

    save_json(
        release=set_release(),
        fasta_files_info=fasta_files_info,
        fasta_files_url=fasta_file_urls,
        file_path=json_file,
        species=species,
        immune_type=immune_type,
        output_dir=path,
        frame_selection=frame_selection,
        create_library_bool=create_library_bool,
        cleanup_bool=cleanup_bool,
        simple_headers_bool=simple_headers_bool
    )


if __name__ == '__main__':
    main("Homo sapiens", "IG" )
