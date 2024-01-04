import logging
import subprocess
import pandas as pd
from pathlib import Path


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

cwd = Path.cwd()
start, stop = 100, 0
reference = cwd / "library" / "retained.fasta"
region_folder = "contig"


def execute(command):
    try:
        logging.info(f"Executing command: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e}")


def parse_bed(file_path):
    entries = set()
    with open(file_path, 'r') as file:
        for line in file:
            chrom, start, end = line.strip().split('\t')[:3]
            entries.add((chrom, start, end))
    return entries


def compare_alignments(base_entries, other_entries):
    new_segments = other_entries - base_entries
    return new_segments


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
        sam_file = primary / f"{prefix}.sam"
        bam_file = primary / f"sorted_{prefix}.bam"
        bed_file = primary / f"{prefix}.bed"
        tags_file = secondary / f"tags_{prefix}.txt"
        tags_fasta_file = secondary / f"mapped_{prefix}.fasta"
        tags_sam_file = secondary / f"mapped_{prefix}.sam"
        tags_bam_file = secondary / f"mapped_{prefix}.bam"
        tags_bed_file = secondary / f"mapped_{prefix}.bed"
        if not bed_file.exists():
            sam_command = f"minimap2 -a -m {accuracy} -t {threads} {reference} {fasta} > {sam_file}"
            bam_command = f"samtools sort -@ {threads} -o {bam_file} {sam_file}"
            index_command = f"samtools index -@ {threads} {bam_file}"
            bed_command = f"bedtools bamtobed -i {bam_file} > {bed_file}"
            tags = f'cat {sam_file} | egrep -v "^@" | cut -f3 | sort | uniq > {tags_file}'
            for command in [sam_command, bam_command, index_command,
                            bed_command, tags]:
                execute(command)

        if not tags_bed_file.exists():
            seqtk_command = f"seqtk subseq {reference} {tags_file} > {tags_fasta_file}"
            sam_mapped = f"minimap2 -a -m {accuracy} -t {threads} {fasta} {tags_fasta_file} > {tags_sam_file}"
            bam_mapped = f"samtools sort -@ {threads} -o {tags_bam_file} {tags_sam_file}"
            index_mapped = f"samtools index -@ {threads} {tags_bam_file}"
            bed_mapped = f"bedtools bamtobed -i {tags_bam_file} > {tags_bed_file}"
            for command in [seqtk_command, sam_mapped, bam_mapped,
                            index_mapped, bed_mapped]:
                execute(command)
        if bed_file.exists():
            entries = parse_bed(bed_file)
            logging.info(f"Parsed {len(entries)} entries from {bed_file}")
            yield prefix, entries
        else:
            logging.warning(f"Required file missing: {bed_file}")
            yield prefix, set()


def write_alignments(comparison_results):
    folder_name = "alignments_excel"
    folder = cwd / folder_name
    make_directory(folder)
    for prefix in comparison_results:
        rows = []
        for accuracy, data in comparison_results[prefix].items():
            row = {
                'Accuracy': accuracy,
                'New vs 100%': data["New vs 100%"],
                'New vs Prev': data["New vs Prev"],
                'New Segments vs 100%': '; '.join(['; '.join(map(str, s)) for s in data["New Segments vs 100%"]]),
                'New Segments vs Prev': '; '.join(['; '.join(map(str, s)) for s in data["New Segments vs Prev"]])
            }
            rows.append(row)

        df = pd.DataFrame(rows)
        excel_filename = cwd / folder / f"{prefix}_comparative_table.xlsx"
        df.to_excel(excel_filename, index=False)
        logging.info(f"Comparative table for {prefix} saved to {excel_filename}")


def main():
    regions = cwd / region_folder
    fasta_files = [f for f in regions.glob("*.fasta")]
    base_accuracy_entries = {prefix: entries for prefix, entries in process_files_for_accuracy(100, fasta_files)}

    comparison_results = {prefix: {} for prefix in base_accuracy_entries}
    prev_accuracy_entries = {}

    for accuracy in range(start, stop - 1, -1):
        if accuracy == 100:
            continue

        for prefix, current_entries in process_files_for_accuracy(accuracy, fasta_files):
            new_segments_100 = compare_alignments(base_accuracy_entries[prefix], current_entries)
            new_segments_prev = set()
            if accuracy + 1 in prev_accuracy_entries:
                new_segments_prev = compare_alignments(prev_accuracy_entries[accuracy + 1].get(prefix, set()), current_entries)
            comparison_results[prefix][f"{accuracy}%"] = {
                "New vs 100%": len(new_segments_100),
                "New vs Prev": len(new_segments_prev),
                "New Segments vs 100%": list(new_segments_100),
                "New Segments vs Prev": list(new_segments_prev)
            }
        prev_accuracy_entries[accuracy] = {prefix: current_entries for prefix, current_entries in process_files_for_accuracy(accuracy, fasta_files)}
        write_alignments(comparison_results)


if __name__ == "__main__":
    main()
