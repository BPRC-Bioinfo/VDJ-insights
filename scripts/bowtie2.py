import re
import logging
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class MappingFiles:
    def __init__(self, prefix, index, beddir):
        self.bowtie_db = index / f"{prefix}_index"
        self.sam = beddir / f"{prefix}.sam"
        self.bam = beddir / f"{prefix}_sorted.bam"
        self.bed = beddir / f"mapped_{prefix}.bed"


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


def parse_bed(file_path, accuracy, fasta, bowtie_type):
    entries = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            region, segment = get_region_and_segment(line[3])
            line.extend([get_sequence(line, fasta),
                        accuracy, str(file_path), bowtie_type, region, segment, fasta])
            entries.append(line)
    return entries


def run_command(command):
    try:
        subprocess.run(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Command '{command}' failed with exit code {e.returncode}")


def make_bowtie2_command(acc, bowtie_db, rfasta, sam_file):
    N = 1 if acc < 50 else 0
    L = int(15 + (acc / 100) * 5)
    score_min_base = -0.1 + (acc / 100) * 0.08
    score_min = f"L,0,{score_min_base:.2f}"
    command = f"bowtie2 -N {N} -L {L} --score-min {score_min} -f -x {bowtie_db} -U {rfasta} -S {sam_file}"
    return command


def make_bowtie_command(acc, bowtie_db, rfasta, sam_file):
    # Map score to mismatches
    # Higher score (towards 100) means fewer mismatches allowed
    mismatches = 3 if acc <= 33 else (2 if acc <= 66 else 1 if acc < 100 else 0)

    # Constructing the Bowtie1 command
    command = f"bowtie -v {mismatches} -m 1 -f -x {bowtie_db} {rfasta} -S {sam_file}"
    return command


def all_commands(files: MappingFiles, fasta_file, rfasta, acc, bowtie_type):
    if bowtie_type == "bowtie2":
        bowtie_command = make_bowtie2_command(
            acc, files.bowtie_db, rfasta, files.sam)
    else:
        bowtie_command = make_bowtie_command(
            acc, files.bowtie_db, rfasta, files.sam)
    commands = [
        f"{bowtie_type}-build {fasta_file} {files.bowtie_db}",
        bowtie_command,
        f"samtools sort -o {files.bam} {files.sam}",
        f"samtools index {files.bam}",
        f"bedtools bamtobed -i {files.bam} > {files.bed}",
    ]
    return commands


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


def run(cwd, indir, outdir, rfasta, beddir, acc, bowtie_type):
    for fasta in indir.glob("*.fasta"):
        prefix = fasta.stem
        index = outdir / prefix
        create_directory(index)
        files = MappingFiles(prefix, index, beddir)
        if not files.bed.exists():
            for command in all_commands(files, fasta, rfasta, acc, bowtie_type):
                run_command(command)
        if files.bed.exists():
            entries = parse_bed(files.bed, acc, fasta, bowtie_type)
            logging.info(f"Parsed {len(entries)} entries from {files.bed}")
            yield entries
        else:
            logging.warning(f"Required file missing: {files.bed}")
            yield list()


def create_directory(location):
    Path(location).mkdir(parents=True, exist_ok=True)


def bowtie2_main(bowtie_type):
    cwd = Path.cwd()
    outdir = cwd / f"{bowtie_type}_db"
    indir = cwd / "contig"
    rfasta = cwd / "library" / "retained.fasta"
    start, stop = 100, 0
    all_entries = []
    for acc in range(start, stop - 1, -1):
        beddir = cwd / bowtie_type / f"{acc}%acc"
        create_directory(beddir)
        for current_entry in run(cwd, indir, outdir, rfasta, beddir, acc, bowtie_type):
            all_entries.extend(current_entry)
    df = make_df(all_entries)
    return df


if __name__ == "__main__":
    bowtie_type = "bowtie"
    bowtie2_main(bowtie_type)
