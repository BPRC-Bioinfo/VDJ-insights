import os
current_cwd = os.getcwd()


def get_ids():
    with open(f"{current_cwd}/input/input.txt") as f:
        input_content = f.readlines()
        id_dict = {}
        for id in input_content:
            id_list = id.strip().split("/")[-2::]
            id_dict[id_list[0]] = id_list[1].split(".")[0]
        return id_dict