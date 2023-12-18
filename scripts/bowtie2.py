import os
import subprocess
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)

def run_command(command):
    try:
        subprocess.run(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command '{command}' failed with exit code {e.returncode}")


def create_directory(location):
    Path(location).mkdir(parents=True, exist_ok=True)
    

def run(cwd, indir, outdir, rfasta, beddir):
    for fasta_file in indir.glob("*.fasta"):
        prefix = fasta_file.stem
        index = outdir / prefix
        create_directory(index)
        db_file = index / f"{prefix}_index"
        sam_file = beddir / f"{prefix}.sam"
        sorted_bam_file = beddir / f"{prefix}_sorted.bam"
        bed_file = beddir / f"mapped_{prefix}.bed"
        db_command = f"bowtie2-build {fasta_file} {db_file}"
        align_command = f"bowtie2 --very-sensitive -N 0 -L 22 --score-min L,0,-0.02 -f -x {db_file} -U {rfasta} -S {sam_file}"
        sort_command = f"samtools sort -o {sorted_bam_file} {sam_file}"
        index_command = f"samtools index {sorted_bam_file}"
        bed_command = f"bedtools bamtobed -i {sorted_bam_file} > {bed_file}"

        for command in [db_command, align_command, sort_command, index_command, bed_command]:
            run_command(command)

    

def main():
    cwd = Path.cwd()
    outdir = cwd / 'bowtie2_db'
    indir = cwd / "contig"
    rfasta = cwd / "library" / "retained.fasta"
    beddir = cwd / "bowtie2" / "100%acc"
    create_directory(beddir)
    run(cwd, indir, outdir, rfasta, beddir)

if __name__ == '__main__':
    main()