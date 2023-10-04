import os
import questionary
from questionary import Style,confirm

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
    with open(path) as f:
        input_content = f.readlines()
        id_dict = {}
        for id in input_content:
            id_list = id.strip().split("/")
            id_dict[id_list[-2]] = id.strip()
        return id_dict



def fetchall_args_sra_download():
    selected_option = questionary.select(
                "Select a option for downloading SRA files",
                choices=["wget", "king", "sra"],
                style=custom_style,
                pointer="❯",
                use_jk_keys=True,
                show_selected=True,
                qmark=""
            ).ask()
    return selected_option

def fetchall_args_input_file():
    input_files = [ f for f in os.listdir(f"{current_cwd}/input/") if f.endswith(".txt")]
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


if __name__ == "__main__":
    print(get_ids("input/input.txt"))