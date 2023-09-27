import subprocess
import os
import shutil


current_pwd = os.getcwd()
output_dir = "reads"

def get_prefetch():
    id = input("Please enter the id of the run or of a file you want to generate fastq files for: ")
    id = id.strip()
    prefetch = f"prefetch {id}"
    print(f"Prefetching {id}!")
    subprocess.call(prefetch, shell=True)


def get_fastq_files():
    for i in os.listdir('.'):
        if i.startswith("SRR"):
            sra_file = f"{current_pwd}/{i}/{i}.sra"
            print (f"Generating fastq for: {i}")
            fastq_dump = f"fasterq-dump  {sra_file}"
            print ("The command used was: " + fastq_dump)
            subprocess.call(fastq_dump, shell=True)


def move_files():
    fastq_files = [f for f in os.listdir('.') if f.endswith('.fastq')]
    if not os.path.exists(f"{current_pwd}/{output_dir}"):
        os.mkdir(f"{current_pwd}/{output_dir}")
    for fastq_file in fastq_files:
        shutil.move(f"{current_pwd}/{fastq_file}", f"{current_pwd}/{output_dir}")

def remove_prefetch():
    for dir in os.listdir('.'):
        if os.path.isdir(f"{current_pwd}/{dir}") and dir != output_dir:
            print(f"Removing {dir}!")
            shutil.rmtree(f"{current_pwd}/{dir}")

if __name__ == "__main__":
    get_prefetch()
    get_fastq_files()
    move_files()
    remove_prefetch()
    
    