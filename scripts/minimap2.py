import logging
import re
import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

cwd = Path.cwd()
start, stop = 100, 0
reference = cwd / "library" / "retained.fasta"
region_folder = "contig"


class MappingFiles:
    def __init__(self, prefix, primary_dir, secondary_dir):
        self.sam = primary_dir / f"{prefix}.sam"
        self.bam = primary_dir / f"sorted_{prefix}.bam"
        self.bed = primary_dir / f"{prefix}.bed"
        self.tags = secondary_dir / f"tags_{prefix}.txt"
        self.tags_fasta = secondary_dir / f"mapped_{prefix}.fasta"
        self.tags_sam = secondary_dir / f"mapped_{prefix}.sam"
        self.tags_bam = secondary_dir / f"mapped_{prefix}.bam"
        self.tags_bed = secondary_dir / f"mapped_{prefix}.bed"


def get_sequence(line, fasta):
    reference, start, stop = line[0:3]
    sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    specific_sequence = sequences[reference].seq
    return str(specific_sequence[int(start): int(stop)])


def get_region_and_segment(name):
    prefix = [i for i in name.split("_") if i.startswith(("TR", "LOC"))][0]
    prefix = re.sub(r"[0-9]", "", prefix)
    if prefix.startswith("TR"):
        return prefix[0:3], prefix[3]
    else:
        return "LOC", "-"


def parse_bed(file_path, accuracy, fasta):
    entries = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            region, segment = get_region_and_segment(line[3])
            line.extend([get_sequence(line, fasta), accuracy,
                        str(file_path), "minimap2", region, segment, fasta])
            entries.append(line)
    return entries


def execute(command):
    try:
        logging.info(f"Executing command: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e}")


def first_mapping_commands(files, accuracy, threads, reference, fasta):
    commands = [
        f"minimap2 -a -m {accuracy} -t {threads} {reference} {fasta} > {files.sam}",
        f"samtools sort -@ {threads} -o {files.bam} {files.sam}",
        f"samtools index -@ {threads} {files.bam}",
        f"bedtools bamtobed -i {files.bam} > {files.bed}",
        f'cat {files.sam} | egrep -v "^@" | cut -f3 | sort | uniq > {files.tags}',
    ]
    return commands


def second_mapping_commands(files, accuracy, threads, reference, fasta):
    commands = [
        f"seqtk subseq {reference} {files.tags} > {files.tags_fasta}",
        f"minimap2 -a -m {accuracy} -t {threads} {fasta} {files.tags_fasta} > {files.tags_sam}",
        f"samtools sort -@ {threads} -o {files.tags_bam} {files.tags_sam}",
        f"samtools index -@ {threads} {files.tags_bam}",
        f"bedtools bamtobed -i {files.tags_bam} > {files.tags_bed}",
    ]
    return commands


def make_directory(directory):
    Path(directory).mkdir(parents=True, exist_ok=True)


def process_files_for_accuracy(accuracy, fasta_files):
    initial_directory = cwd / "segments_mapping" / f"{accuracy}%acc"
    primary = initial_directory / "primary_mapping"
    secondary = initial_directory / "secondary_mapping"
    threads = 10
    for dir in [initial_directory, primary, secondary]:
        make_directory(dir)

    for fasta in fasta_files:
        prefix = fasta.stem
        fasta = cwd / region_folder / fasta
        files = MappingFiles(prefix, primary, secondary)
        if not files.bed.exists():
            for command in first_mapping_commands(
                files, accuracy, threads, reference, fasta
            ):
                execute(command)
        if not files.tags_bed.exists():
            for command in second_mapping_commands(
                files, accuracy, threads, reference, fasta
            ):
                execute(command)
        if files.tags_bed.exists():
            entries = parse_bed(files.tags_bed, accuracy, fasta)
            logging.info(f"Parsed {len(entries)} entries from {files.tags_bed}")
            yield entries
        else:
            logging.warning(f"Required file missing: {files.tags_bed}")
            yield list()


def make_df(all_entries):
    headers = [
        "reference",
        "start",
        "stop",
        "name",
        "score",
        "strand",
        "sequence",
        "accuracy",
        "file",
        "tool",
        "region",
        "segment",
        "fasta-file"
    ]
    df = pd.DataFrame(all_entries, columns=headers)
    unique_combinations = df.drop_duplicates(subset=["name", "start", "stop"])
    return unique_combinations.reset_index(drop=True)


def minimap2_main():
    regions = cwd / region_folder
    fasta_files = [f for f in regions.glob("*.fasta")]
    all_entries = []
    for accuracy in range(start, stop - 1, -1):
        for current_entries in process_files_for_accuracy(accuracy, fasta_files):
            if current_entries:
                all_entries.extend(current_entries)
    return make_df(all_entries)


if __name__ == "__main__":
    minimap2_main()
