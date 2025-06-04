import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

from .logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)

def scrape_imgt(species: str, receptor: str, library_file: Path) -> dict:
    segments = {
        "TR": ["TRBV", "TRAV", "TRDV", "TRGV"],
        "IG": ["IGHV", "IGKV", "IGLV"]
    }
    for segment in segments[receptor]:
        species = species.replace(" ", "_")
        url = f"https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/{species}/{receptor}/{segment}.fasta"
        response = requests.get(url)
        if response.status_code == 200:
            fasta_io = StringIO(response.text)
            sequences = SeqIO.to_dict(SeqIO.parse(fasta_io, "fasta"))
            with open(library_file, 'a') as output_handle:
                SeqIO.write(sequences.values(), output_handle, "fasta")
        else:
            console_log.warning(f"{url} fails ({response.status_code}).")

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def main_cdr(species: str, receptor: str,  threads: int = 12) -> None:
    """
    Hoofdfunctie om de CDR-annotaties te verwerken.
    """
    cwd = Path.cwd()
    output_base = cwd / "tmp/CDR"
    output_base.mkdir(parents=True, exist_ok=True)

    library_path = output_base / "library"
    library_path.mkdir(parents=True, exist_ok=True)
    library_file = library_path / "library.fasta"

    if not Path(library_file).is_file():
        scrape_imgt(species=species, receptor=receptor, library_file=library_file)

    cdr1_library = library_path / "CDR1.fasta"
    cdr2_library = library_path / "CDR2.fasta"

    with (open(cdr1_library, 'w') as cdr1_fasta, open(cdr2_library, 'w') as cdr2_fasta):
        for record in SeqIO.parse(library_file, 'fasta'):
            header = record.description
            seq = record.seq
            allele = header.split("|")[1]

            cdr1_sequence = seq[78:114].upper().replace(".", "")
            cdr2_sequence = seq[165:195].upper().replace(".", "")

            cdr1_fasta.write(f">{allele}_cdr1\n{cdr1_sequence}\n")
            cdr2_fasta.write(f">{allele}_cdr2\n{cdr2_sequence}\n")

    output_blast_result = output_base / "blast"
    output_blast_result.mkdir(exist_ok=True, parents=True)

    inpute_target = output_blast_result / "report"
    inpute_target.mkdir(exist_ok=True, parents=True)

    annotation_report_path = cwd / "annotation" / "annotation_report_all.xlsx"
    annotation_report = pd.read_excel(annotation_report_path)

    target_lengths = {}
    target_sequence = inpute_target / "target_sequence.fasta"
    with open(target_sequence, 'w') as target_fasta:
        for index, row in annotation_report[["Short name", "Target sequence"]].iterrows():
            target_fasta.write(f">{row['Short name']}_{index}\n{row['Target sequence'].replace('-', '')}\n")
            target_lengths[index] = len(row['Target sequence'].replace('-', ''))

    for library in [cdr1_library, cdr2_library]:
        cdr = library.stem
        output_blast_file = output_blast_result / f"{cdr}.txt"
        blast_cmd = f'blastn -task blastn-short -query {library} -subject {target_sequence} -qcov_hsp_perc 100 -out {output_blast_file} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs"'
        subprocess.run(blast_cmd, shell=True, check=True)

        if output_blast_file.exists():
            columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                       "send", "evalue", "bitscore", "qseq", "sseq", "qcovs"]
            blast_results = pd.read_csv(output_blast_file, sep="\t", names=columns)
            if not blast_results.empty:
                grouped_blast_results = blast_results.groupby("sseqid")
                for index, blast_group in grouped_blast_results:
                    s_index = index.split("_")[0]

                    if cdr == "CDR1":
                        sorted_blast_group = blast_group.sort_values(by=["qstart", "qcovs", "length", "sstart"], ascending=[True, False, False, True])
                        if s_index in sorted_blast_group['qseqid'].values:
                            best_hit = sorted_blast_group[sorted_blast_group['qseqid'].str.contains(s_index, regex=False)].iloc[0]
                        else:
                            best_hit = sorted_blast_group.iloc[0]

                    elif cdr == "CDR2":
                        sorted_blast_group = blast_group.sort_values(by=["qstart", "qcovs", "length", "sstart"], ascending=[True, False, False, False])
                        if sorted_blast_group['qseqid'].str.contains(s_index, regex=False).any():
                            best_hit = sorted_blast_group[sorted_blast_group['qseqid'].str.contains(s_index, regex=False)].iloc[0]
                        else:
                            best_hit = sorted_blast_group.iloc[0]

                    index_segment = int(best_hit['sseqid'].split("_")[-1])
                    sstart = int(best_hit["sstart"])
                    send = int(best_hit["send"])

                    if sstart > send:
                        target_length = target_lengths[index_segment]
                        sstart = target_length - sstart + 1
                        send = target_length - send + 1

                    annotation_report.loc[index_segment, f"{cdr}_start"] = sstart - 1
                    annotation_report.loc[index_segment, f"{cdr}_stop"] = send

    known = annotation_report[annotation_report["Status"] == "Known"]
    novel = annotation_report[annotation_report["Status"] == "Novel"]

    if not annotation_report.empty:
        combined_df = annotation_report.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
        combined_df.to_excel(cwd / "annotation" / "annotation_report_all.xlsx", index=False)
    if not known.empty:
        known = known.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
        known.to_excel(cwd / "annotation" / "annotation_report_known.xlsx", index=False)
    if not novel.empty:
        novel = novel.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
        novel.to_excel(cwd / "annotation" / "annotation_report_novel.xlsx", index=False)
