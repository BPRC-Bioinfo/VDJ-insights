import os
import subprocess
import questionary
from questionary import Style
import json

current_cwd = os.getcwd()

custom_style = Style([
    ('separator', 'fg:#6C6C6C'),
    ('qmark', 'fg:#FF9D00 bold'),  # question mark at the start of the prompt
    ('question', 'fg:#9ad4f5'),  # question text
    ('selected', 'fg:#5F819D'),  # for the selected item
    ('pointer', 'fg:#fff bold'),  # pointer used in select and checkbox prompts
    ('answer', 'fg:#fe2e2e bold'),  # for answer
    ('highlighted', 'fg:#a6e0ff bold'),
])


def get_ids(path):
    """
    creates a dictionary with as key the SRA identifier and as value
    a ENA link.

    Parameters:
        path (str): path to the input.txt file to be used.

    Returns:
        id_dict (dict): dictionary with SRA identifier and ENA link.
    """
    with open(path) as f:
        input_content = f.readlines()
        id_dict = {}
        for id in input_content:
            id_list = id.strip().split("/")
            id_dict[id_list[-2]] = id.strip()
        return id_dict


def fetch_chromosome():
    chromosomes = {}
    for chromosome in ["chr3", "chr7"]:
        assembly_file = f"{current_cwd}/downloads/reports/mmul10_assembly_report.txt"
        egrep_cmd = f"cat {assembly_file} | egrep '{chromosome}' | cut -f 5"
        result = subprocess.getoutput(egrep_cmd).split("\n")
        chromosomes[chromosome] = result
    with open(f"{current_cwd}/input/chromosome_conversion.json", "w") as f:
        json.dump(chromosomes, f, indent=4)


def cal_chr_length():
    mmul10 = {}
    with open("downloads/mmul10.fna", "r") as f:
        start = None
        for line in f:
            if line.startswith(">"):
                line = line.split()
                start = f"{line[0][1::]}"
                mmul10[start] = 0
            else:
                mmul10[start] += len(line.strip())
    with open(f"{current_cwd}/input/chromosome_lengths.json", "w") as f:
        json.dump(mmul10, f, indent=4)


def fetchall_args_sra_download():
    """
    Visual prompt for asking which download type you want to use as
    download type: wget, king or sra.
    Returns:
        selected_option (str): chosen option as string.
    """
    selected_option = questionary.select(
                "Select a option for downloading SRA files",
                choices=["wget", "kingfisher", "sra"],
                style=custom_style,
                pointer="❯",
                use_jk_keys=True,
                show_selected=True,
                qmark=""
            ).ask()
    return selected_option


def fetchall_args_input_file():
    """
    Visual prompt for asking which input file you want to use.
    Returns:
        selected_option (str): chosen input file as string.
    """
    input_files = [f for f in os.listdir(f"{current_cwd}/input/") if f.endswith(".txt")]
    selected_input = questionary.select(
                "Select a option as input file",
                choices=input_files,
                style=custom_style,
                pointer="❯",
                use_jk_keys=True,
                show_selected=True,
                qmark=""
            ).ask()
    return selected_input
