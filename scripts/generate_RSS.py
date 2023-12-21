import re
import logging
import tempfile
import subprocess
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
import json

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

cwd = Path.cwd()
accession = "EAW"
folder = "contig"
start, stop = 100, 85
all_info = []
positive, negative, separated_segments = {}, {}, {}


def write_fasta_file(dictionary, folder):
    for key, value in dictionary.items():
        for region, segments in value.items():
            ffile = cwd / folder / f"{key}{region}.fasta"
            with open(ffile, 'w') as f:
                for header, segment in segments.items():
                    for x, rss in enumerate(segment):
                        f.write(f">{header}-{x+1}\n{rss}\n")


def create_segments_dict():
    Path(cwd / "RSS").mkdir(parents=True, exist_ok=True)
    merged_data = {**positive, **negative}
    for key, value in merged_data.items():
        prefix = [i for i in key.split("_") if i.startswith(("TR", "LOC"))][0]
        prefix = re.sub(r"[0-9]", "", prefix)
        if prefix.startswith("TR"):
            region, segment = prefix[0:3], prefix[3]
            separated_segments.setdefault(region, {}).setdefault(segment, {})[key] = value


def add_to_dict(name, dictionary, rss):
    dictionary.setdefault(name, []).append(rss)


def create_sequence(bed_content, original_file):
    base_fasta = "_".join(original_file.stem.split("_")[1:]) + ".fasta"
    ref_fasta = cwd / folder / base_fasta
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


def create_rss(combination):
    name, reference, begin, end, strand, filename = combination
    TR_seg = {
        "V": {"+": [end, str(int(end) + 39)], "-": [str(int(begin) - 39), begin]},
        "D": {"+": [str(int(begin) - 28), str(int(end) + 28)], "-": [str(int(begin) - 28), str(int(end) + 28)]},
        "J": {"+": [str(int(begin) - 39), begin], "-": [end, str(int(end) + 39)]}
    }
    prefix = [i for i in name.split("_") if i.startswith(("TR", "LOC"))][0]
    prefix = re.sub(r"[0-9]", "", prefix)
    if prefix.startswith("TR"):
        segment = prefix[3]
        rss_bed = TR_seg[segment][strand]
        bed_content = [reference, rss_bed[0], rss_bed[-1]]
        rss = create_sequence(bed_content, Path(filename))
        if strand == '+':
            add_to_dict(name, positive, rss)
        elif strand == '-':
            rss = Seq(rss)
            rss = str(rss.reverse_complement())
            add_to_dict(name, negative, rss)


def make_df():
    headers = ["name", "reference", "start", "stop", "strand", "file"]
    df = pd.DataFrame(all_info, columns=headers)
    unique_combinations = df.drop_duplicates(subset=['name', 'start', 'stop'])
    return unique_combinations.values.tolist()


def parse_bedfile(folder, filename):
    bedfile = cwd / folder / filename
    with open(bedfile, 'r') as f:
        for line in f:
            sline = line.split()
            ref, start, stop, name, strand = (sline[0], sline[1],
                                              sline[2], sline[3],
                                              sline[-1])
            all_info.append([name, ref, start, stop, strand, str(filename)])


def gather_list():
    logging.info(f"Fetching RSS sequences from files in the range of {start}%acc to {stop}%acc.")
    for accuracy in range(start, stop - 1, -1):
        logging.info(f"Retrieving RSS sequences from the {accuracy}%acc bedfile!")
        minimap2_folder = cwd / "segments_mapping" / f"{accuracy}%acc" / "secondary_mapping"
        for minibfile in minimap2_folder.glob(f"*_{accession}*.bed"):
            parse_bedfile(minimap2_folder, minibfile)
    bowtie2_folder = cwd / "bowtie2" / f"100%acc"
    for bowbfile in bowtie2_folder.glob(f"*_{accession}*.bed"):
        parse_bedfile(bowtie2_folder, bowbfile)
    return make_df()


def main():
    unique = gather_list()
    for i in unique:
        create_rss(i)
    create_segments_dict()
    write_fasta_file(separated_segments, "RSS")


if __name__ == '__main__':
    main()
