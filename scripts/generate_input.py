import os
import json
import subprocess
from random import choice

current_cwd = os.getcwd()


def generate_input_file(file):
    with open(f"{current_cwd}/input/{file}") as f:
        input_content = [file_id["fastq_ftp"] for file_id in json.loads(f.read())]
        top_five = [choice(input_content) for _ in range(5)]
        with open(f"{current_cwd}/input/input.txt", "w") as f:
            f.write("\n".join(top_five))
        


for file in os.listdir(f"{current_cwd}/input"):
    if file.endswith("_json.txt"):
        json_file = f"{current_cwd}/input"
        new_name = {json_file}/{file.split('_json')[0]}.json
        subprocess.call(f"mv {json_file}/{file} {new_name}", shell=True)
        id_dict = generate_input_file(file) 
    if file.endswith(".json"):
        id_dict = generate_input_file(file)


# id_dict = get_ids()

