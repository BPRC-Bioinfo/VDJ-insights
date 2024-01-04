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


def parse_bed(file_path, accuracy, fasta):
    entries = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            region, segment = get_region_and_segment(line[3])
            line.extend([get_sequence(line, fasta),
                        accuracy, str(file_path), "bowtie2", region, segment, fasta])
            entries.append(line)
    return entries


def run_command(command):
    try:
        subprocess.run(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Command '{command}' failed with exit code {e.returncode}")


def all_commands(files: MappingFiles, fasta_file, rfasta):
    commands = [
        f"bowtie2-build {fasta_file} {files.bowtie_db}",
        f"bowtie2 --very-sensitive -N 0 -L 22 --score-min L,0,-0.02 -f -x {files.bowtie_db} -U {rfasta} -S {files.sam}",
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


def run(cwd, indir, outdir, rfasta, beddir):
    for fasta in indir.glob("*.fasta"):
        prefix = fasta.stem
        index = outdir / prefix
        create_directory(index)
        files = MappingFiles(prefix, index, beddir)
        if not files.bed.exists():
            for command in all_commands(files, fasta, rfasta):
                run_command(command)
        if files.bed.exists():
            entries = parse_bed(files.bed, 100, fasta)
            logging.info(f"Parsed {len(entries)} entries from {files.bed}")
            yield entries
        else:
            logging.warning(f"Required file missing: {files.bed}")
            yield list()


def create_directory(location):
    Path(location).mkdir(parents=True, exist_ok=True)


def bowtie2_main():
    cwd = Path.cwd()
    outdir = cwd / "bowtie2_db"
    indir = cwd / "contig"
    rfasta = cwd / "library" / "retained.fasta"
    beddir = cwd / "bowtie2" / "100%acc"
    all_entries = []
    create_directory(beddir)
    for current_entry in run(cwd, indir, outdir, rfasta, beddir):
        all_entries.extend(current_entry)
    df = make_df(all_entries)
    return df


if __name__ == "__main__":
    bowtie2_main()
