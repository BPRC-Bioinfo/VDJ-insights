import os
import tempfile
import subprocess
from Bio.Seq import Seq
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

cwd = os.getcwd()
accession = "EAW"
folder = "contig"
start, stop = 100, 85
positive = {}
negative = {}

def write_temp_file(bed_content, original_file):
    base_fasta = "_".join(original_file.split(".")[0].split("_")[1:]) + ".fasta"
    ref_fasta = os.path.join(cwd, folder, base_fasta)
    extension = '.bed'
    extension2 = '.fasta'
    with tempfile.NamedTemporaryFile(suffix=extension, mode="w+", delete=True) as temp:
        temp.write("\t".join(bed_content))
        temp.flush()
        temp.seek(0)
        temp_file_path = temp.name

        with tempfile.NamedTemporaryFile(suffix=extension2, mode="w+", delete=True) as temp2:
            temp2_file_path = temp2.name
            command = f"bedtools getfasta -fi {ref_fasta} -bed {temp_file_path} -fo {temp2_file_path}"
            subprocess.call(command, shell=True)

            temp2.seek(0)
            sequence = temp2.readlines()[-1].strip()
            return sequence


def add_to_dict(name, dictionary, rss):
    if name not in dictionary:
        dictionary[name] = [rss]
    else:
        if rss not in dictionary[name]:
            dictionary[name].append(rss)
 

def parse_bedfile(folder, filename):
    bedfile = os.path.join(cwd, folder, filename)
    with open(bedfile, 'r') as f:
        for line in f:
            sline = line.split()
            ref, start, stop, name, strand = (sline[0], sline[1], 
                                              sline[2], sline[3], 
                                              sline[-1])
            if strand == '+':
                rss = write_temp_file([ref, stop, str(int(stop)+39)], filename)
                add_to_dict(name, positive, rss)
            elif strand == '-':
                rss = write_temp_file([ref, str(int(start)-39), start], filename)
                rss = Seq(rss)
                rss = str(rss.reverse_complement())
                add_to_dict(name, negative, rss)
                

def write_fasta_file(dictionary, folder, filename):
    ffile = os.path.join(cwd, folder, filename)
    with open(ffile, 'w') as f:
        for key, value in dictionary.items():
            for x, rss in enumerate(value):
                f.write(f">{key}-{x+1}\n{rss}\n")
        
    
    
def write_json_file():
    rss_folder = os.path.join(cwd, "RSS")
    if not os.path.exists(rss_folder):
        os.mkdir(rss_folder)
    logging.info("Writing RSS sequences from the positive strand to positive.json!")
    with open(os.path.join(rss_folder, "positive.json"), "w") as f:
        f.write(json.dumps(positive, indent=4))
    logging.info("Writing RSS sequences from the negative strand to negative.json!")
    with open(os.path.join(rss_folder, "negative.json"), "w") as f:
        f.write(json.dumps(negative, indent=4))
    complete = {**positive, **negative}
    write_fasta_file(complete, rss_folder, "complete.fasta")
    
def main():
    logging.info(f"Fetching RSS sequences from files in the range of {start}%acc to {stop}%acc.")
    for accuracy in range(start, stop - 1, -1):
        logging.info(f"Retrieving RSS sequences from the {accuracy}%acc bedfile!")
        folder = os.path.join(cwd, "segments_mapping", f"{accuracy}%acc", "secondary_mapping")
        for bfile in os.listdir(folder):
            if accession in bfile and bfile.endswith("bed"):
                parse_bedfile(folder, bfile)
    write_json_file()

if __name__ == '__main__':
    main()
        