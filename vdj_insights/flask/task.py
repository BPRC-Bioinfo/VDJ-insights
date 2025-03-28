from time import sleep

from celery import shared_task
import subprocess
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from click import command


def get_fasta(file_names):
    sequences = []
    for file_name in file_names:
        for record in SeqIO.parse(file_name, "fasta"):
            sequences.append({"id": record.id, "sequence": str(record.seq)})
    return sequences


@shared_task
def merge_sequences(file_names: list, output_dir: str):
    merged_dir = os.path.join(output_dir, "merged")
    aligned_dir = os.path.join(output_dir, "aligned")
    split_dir = os.path.join(output_dir, "split_regions")

    os.makedirs(merged_dir, exist_ok=True)
    os.makedirs(aligned_dir, exist_ok=True)
    os.makedirs(split_dir, exist_ok=True)

    merged_path = os.path.join(merged_dir, "merged_sequences.fasta")
    aligned_path = os.path.join(aligned_dir, "aligned_sequences.fasta")

    sequences = get_fasta(file_names)
    records = [SeqRecord(Seq(seq["sequence"]), id=seq["id"], description="") for seq in sequences]
    with open(merged_path, "w") as f:
        SeqIO.write(records, f, "fasta")

    command = ["mafft", "--auto", merged_path]
    with open(aligned_path, "w") as output:
        subprocess.run(command, stdout=output, stderr=subprocess.PIPE, check=True)

    aligned_sequences = get_fasta([aligned_path])
    records = [SeqRecord(Seq(seq["sequence"]), id=seq["id"], description="") for seq in aligned_sequences]

    for record in records:
        fasta_file = os.path.join(split_dir, f"{record.id}.fasta")
        with open(fasta_file, "w") as f:
            SeqIO.write(record, f, "fasta")
    sleep(10)
    return split_dir


if __name__ == '__main__':
    paths = ["data1/region/sample_1.fa", "data1/region/sample_2.fa"]
    output_dir = "data1/region/"
    regions = merge_sequences(paths, output_dir)
    print(regions)