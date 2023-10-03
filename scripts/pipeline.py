import os
current_cwd = os.getcwd()


def get_ids(path):
    with open(path) as f:
        input_content = f.readlines()
        id_dict = {}
        for id in input_content:
            id_list = id.strip().split("/")
            id_dict[id_list[-2]] = id.strip()
        return id_dict

if __name__ == "__main__":
    print(get_ids("input/input.txt"))