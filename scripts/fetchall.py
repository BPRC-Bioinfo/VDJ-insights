import subprocess
import os
import shutil
import argparse
from argparse import RawTextHelpFormatter
from pipeline import get_ids

current_pwd = os.getcwd()
output_dir = "reads"


def parser_args():
    parser = argparse.ArgumentParser(
        description=f"Choosing a downloading style", formatter_class=RawTextHelpFormatter)
    group1 = parser.add_argument_group('required arguments')
    group1.add_argument("-t", "--type", nargs="+", type=str, required=True, metavar="type", choices=["sra", "kingfisher", "wget"], default="wget",
                        help="choose a downloading style.\nOptions include sra (sra), kingfisher (king) or wget")
    group1.add_argument("-i", "--input", nargs="+", type=str, required=True, metavar="input", 
                        help="choose a input file")
    # group2 = parser.add_argument_group('optional arguments')

    args = parser.parse_args()
    return args, parser

def sra():
    print(f"Fetching ids with sra!")
    sra_get_prefetch()
    sra_get_fastq_files()
    

def sra_get_prefetch(id):
    prefetch = f"prefetch {id}"
    print(f"Prefetching {id}!")
    subprocess.call(prefetch, shell=True)


def sra_get_fastq_files():
    for i in os.listdir('.'):
        if i.startswith("SRR"):
            sra_file = f"{current_pwd}/{i}/{i}.sra"
            print (f"Generating fastq for: {i}")
            fastq_dump = f"fasterq-dump  {sra_file}"
            print ("The command used was: " + fastq_dump)
            subprocess.call(fastq_dump, shell=True)

def kingfisher(id):
    print("Fetching files with kingfisher!")
    kingfisher_fetch = f"kingfisher get -r {id} -m ena-ftp"
    print(f"Fetching {id} with kingfisher!")
    subprocess.call(kingfisher_fetch, shell=True)
    

def wget(id):
    print("Fetching files with wget!")
    wget_fetch = f"wget {id}"
    subprocess.call(wget_fetch, shell=True)

# def move_files():
#     fastq_files = [f for f in os.listdir('.') if f.endswith('.fastq')]
#     if not os.path.exists(f"{current_pwd}/{output_dir}"):
#         os.mkdir(f"{current_pwd}/{output_dir}")
#     for fastq_file in fastq_files:
#         shutil.move(f"{current_pwd}/{fastq_file}", f"{current_pwd}/{output_dir}")

# def remove_prefetch():
#     for dir in os.listdir('.'):
#         if os.path.isdir(f"{current_pwd}/{dir}") and dir != output_dir:
#             print(f"Removing {dir}!")
#             shutil.rmtree(f"{current_pwd}/{dir}")

if __name__ == "__main__":
    args, parser = parser_args() 
    ids = get_ids(path=args.input[0])
    if args.type[0] == "sra":
        for key, value in ids.items():
            sra(key)
    elif args.type[0] == "kingfisher":
        for key, value in ids.items():
            kingfisher(key)
    elif args.type[0] == "wget":
        with open(args.input[0], "r") as f:
            for line in f:
                wget(line.strip())
    else:
        print(f"{args.type[0]} is not a valid option!")
        
    # move_files()
    # remove_prefetch()
    
    