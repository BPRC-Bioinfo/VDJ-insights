import subprocess
import os
import shutil
import argparse
from scripts.pipeline import get_ids

current_pwd = os.getcwd()
output_dir = "reads"


def sra(id):
    for id in get_ids().keys():
        sra_get_prefetch(id)
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


def kingfisher():
    for id in get_ids().keys():
        kingfisher_fetch = f"kingfisher get -r {id} -m ena-ftp"
        print(f"Fetching {id} with kingfisher!")
        subprocess.call(kingfisher_fetch, shell=True)


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
    pass
    # move_files()
    # remove_prefetch()
    
    