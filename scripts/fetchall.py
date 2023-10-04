import subprocess
import os
import shutil
import argparse
from argparse import RawTextHelpFormatter
import gzip

current_pwd = os.getcwd()
output_dir = "downloads"


def parser_args():
    parser = argparse.ArgumentParser(
        description=f"Choosing a downloading style", formatter_class=RawTextHelpFormatter)
    group1 = parser.add_argument_group('required arguments')
    group1.add_argument("-t", "--type", nargs="+", type=str, required=True, metavar="type", choices=["sra", "kingfisher", "wget"], default="wget",
                        help="choose a downloading style.\nOptions include sra (sra), kingfisher (king) or wget")
    group1.add_argument("-i", "--input", nargs="+", type=str, required=True, metavar="input", 
                        help="choose a input file")
    group1.add_argument("-o", "--output", nargs="+", type=str, required=True, metavar="output", 
                        help="choose a output file location")
     
    # group2 = parser.add_argument_group('optional arguments')

    args = parser.parse_args()
    return args, parser

def gunzip_file(file, new_file):
    with open(file, 'rb') as f_in:
        with gzip.open(f"{new_file}.fastq.gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def clean(fastq_file: str):
    return fastq_file.split("/")[-1].split("_")[0]

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

def move_files(file, location):
    if file:
        current_file = [i for i in os.listdir(".") if file in i][0]
        if current_file.endswith(".fastq.gz"):
            print(f"Renaming and moving {current_file} to {location}!")
            shutil.move(f"{current_pwd}/{current_file}", f"{current_pwd}/{location}")
        else:
            gunzip_file(current_file, current_file)
    else:
        print("No fastq files found!")
            

def remove_prefetch():
    for dir in os.listdir('.'):
        if os.path.isdir(f"{current_pwd}/{dir}") and dir != output_dir:
            print(f"Removing {dir}!")
            shutil.rmtree(f"{current_pwd}/{dir}")

def run():
    print("Running...")
    args, parser = parser_args() 
    print(args)
    if args.type[0] == "sra":
        sra(clean(args.input[0]))
    elif args.type[0] == "kingfisher":
        kingfisher(clean(args.input[0]))
    elif args.type[0] == "wget":
        wget(args.input[0])
    else:
        print(f"{args.type[0]} is not a valid option!")  
    move_files(clean(args.input[0]), args.output[0])
    # remove_prefetch()
    

if __name__ == "__main__":
    run()