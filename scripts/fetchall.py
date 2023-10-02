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
    # group2.add_argument("-sep", "--separator", type=str, choices=[",", "t", "|", ";", " "], default=",", metavar="SEPARATOR",
    #                     help="Choose a separator for the output file.\nOptions include comma (,), tab (t), pipe (|), semicolon (;), or space ( ). Default is comma ','.")
    # group2.add_argument("-l", "--list_datasets", action="store_true",
    #                     help="Shows a list of all populations which have a haplotype dataset")

    args = parser.parse_args()
    return args, parser

def sra():
    print(f"Fetching ids with sra!")
    # for id in get_ids().keys():
    #     sra_get_prefetch(id)
    # sra_get_fastq_files()
    

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
    print("Fetching ids with kingfisher!")
    # kingfisher_fetch = f"kingfisher get -r {id} -m ena-ftp"
    # print(f"Fetching {id} with kingfisher!")
    # subprocess.call(kingfisher_fetch, shell=True)
    

def wget(id):
    print("Fetching ids with wget!")
    # wget_fetch = f"wget -r {id} -m ena-ftp"
    # print(f"Fetching {id} with wget!")
    # subprocess.call(wget_fetch, shell=True)


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
    
    if args.type[0] == "sra":
        sra()
    elif args.type[0] == "kingfisher":
        kingfisher()
    elif args.type[0] == "wget":
        wget()
        
    # move_files()
    # remove_prefetch()
    
    