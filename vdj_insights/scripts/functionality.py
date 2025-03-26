import subprocess
from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings
from Bio import BiopythonWarning

from .IMGT_scrape import main as imgt_main

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

warnings.simplefilter('ignore', BiopythonWarning)


def open_files(cwd: Path) -> (pd.DataFrame, pd.DataFrame):
    data = cwd / "annotation" / "annotation_report_all.xlsx"

    data = pd.read_excel(data)

    data_V = data[data["Segment"].isin(["V"])].copy()
    data_D = data[data["Segment"].isin(["D"])].copy()
    data_J = data[data["Segment"].isin(["J"])].copy()
    return data_V, data_D, data_J


def write_l_region_fasta_files(library_file, library_path):
    l_parts = {}

    for record in SeqIO.parse(library_file, 'fasta'):
        header = record.description
        seq = record.seq

        allele = header.split("|")[1]
        header_positions = record.description.split("|")[5]

        cumulative_offset = 0
        for index, positions in enumerate(header_positions.split("+")):
            start = int(positions.split("..")[0])
            end = int(positions.split("..")[1])
            segment_length = end - start + 1

            segment_seq = seq[cumulative_offset: cumulative_offset + segment_length]
            cumulative_offset += segment_length

            part_key = f"L_PART{index + 1}"
            if part_key not in l_parts:
                l_parts[part_key] = {}
            l_parts[part_key][allele] = segment_seq

    for part, allele_dict in l_parts.items():
        records = []
        for allele, seq in allele_dict.items():
            record = SeqRecord(Seq(seq), id=allele, description="")
            records.append(record)

        with open(f"{library_path}/{part}.fasta", "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")


def write_extract_region_from_genome(locus_fasta_path, region, l_part, length, group_locus):
    locus_fasta_file_name = locus_fasta_path / f"{l_part}_{region}.fasta"
    with (open(locus_fasta_file_name, 'w') as l_file):
        for index_segment, row in group_locus.iterrows():

            target_name = row["Target name"]
            start_coord_segment = row["Start coord"]
            end_coord_segment = row["End coord"]
            strand = row["Strand"]
            fasta_path = row["Path"]

            with open(fasta_path, 'r') as fasta_file:
                sequence_region = SeqIO.read(fasta_file, 'fasta')

            if strand == "+":
                v_downstream_region = sequence_region.seq[(start_coord_segment - length):start_coord_segment]
            elif strand == "-":
                v_downstream_region = sequence_region.seq[end_coord_segment:(end_coord_segment + length)]
            if len(v_downstream_region) != 0:
                l_file.write(f">{target_name}_{index_segment}\n{v_downstream_region}\n")
    return locus_fasta_file_name


def get_blast_results(output_base, l_part, region, locus_fasta_file_name, library_path, extract_coords):
    output_blast_result = output_base / "blast"
    output_blast_result.mkdir(exist_ok=True, parents=True)
    output_blast_result = output_blast_result / f"{l_part}_{region}.txt"

    library_file = library_path / f"{l_part}.fasta"

    blast_cmd = f'blastn -task blastn-short -word_size 7 -reward 1  -penalty -2 -gapopen 5 -gapextend 2 -dust no  -query {library_file} -subject {locus_fasta_file_name} -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -out {output_blast_result} -strand both -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qcovs"'
    subprocess.run(blast_cmd, shell=True, check=True)

    if output_blast_result.exists():
        columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq", "qcovs"]
        blast_results = pd.read_csv(output_blast_result, sep="\t", names=columns)
        if not blast_results.empty:
            grouped_blast_results = blast_results.groupby("sseqid")

            length_l_part = {record.id: len(record.seq) for record in SeqIO.parse(library_file, 'fasta')}
            for index, blast_group in grouped_blast_results:
                if index == "":
                    print(blast_group.sort_values(by=["qstart", "qcovs", "length", "mismatch"], ascending=[True, False, False, True]).head())
                best_hit = blast_group.sort_values(by=["qstart", "qcovs", "length", "mismatch"], ascending=[True, False, False, True]).iloc[0]

                sseqid = best_hit["sseqid"]
                if best_hit["qcovs"] == 100 and best_hit["qstart"] == 1:
                    start_coord = best_hit["sstart"]
                    end_coord = best_hit["send"]
                elif best_hit["qcovs"] > 70:
                    start_coord = (best_hit["sstart"] - best_hit["qstart"]) + 1
                    end_coord = best_hit["send"] + (length_l_part[best_hit["qseqid"]] - best_hit["qend"])
                else:
                    start_coord = None
                    end_coord = None

                if start_coord is not None and end_coord is not None :
                    if sseqid not in extract_coords:
                        extract_coords[sseqid] = {}
                    extract_coords[sseqid][f"{l_part}_start"] = start_coord
                    extract_coords[sseqid][f"{l_part}_end"] = end_coord
    return extract_coords


def extraxt_region_from_genome(locus_fasta_path, region, group_locus, extract_coords, off_set):
    locus_fasta_file_name = locus_fasta_path / f"L_PART1_{region}.fasta"
    extracted_region = {record.id: record.seq.upper() for record in SeqIO.parse(locus_fasta_file_name, 'fasta')}

    for seq_id, coords in extract_coords.items():
        L_PART1_start, L_PART1_end = coords.get("L_PART1_start"), coords.get("L_PART1_end")
        L_PART2_start, L_PART2_end = coords.get("L_PART2_start"), coords.get("L_PART2_end")
        if all([L_PART1_start, L_PART1_end, L_PART2_start, L_PART2_end]):
            sequence = extracted_region[seq_id]
            index_segment = int(seq_id.split("_")[-1])

            if L_PART1_start > L_PART1_end or L_PART2_start > L_PART2_end:
                sequence = sequence.reverse_complement()
                length_sequence = len(sequence)

                L_PART1_seq = sequence[(length_sequence - L_PART1_start):(length_sequence - L_PART1_end + 1)]
                L_PART2_seq = sequence[(length_sequence - L_PART2_start):(length_sequence - L_PART2_end + 1)]
                DONOR_SPLICE = sequence[(length_sequence - L_PART1_end + 1):(length_sequence - L_PART1_end + 3)]
                ACCEPTOR_SPLICE = sequence[(length_sequence - L_PART2_start - 2):(length_sequence - L_PART2_start)]

                target_sequence = Seq(group_locus.loc[index_segment, "Target sequence"].replace("-", "")).reverse_complement()
            else:
                L_PART1_seq = sequence[(L_PART1_start - 1):L_PART1_end]
                L_PART2_seq = sequence[(L_PART2_start + off_set - 1):(L_PART2_end + off_set)]
                DONOR_SPLICE = sequence[L_PART1_end:(L_PART1_end + 2)]
                ACCEPTOR_SPLICE = sequence[(L_PART2_start + off_set - 3):(L_PART2_start + off_set - 1)]

                target_sequence = group_locus.loc[index_segment, "Target sequence"].replace("-", "")

            group_locus.loc[index_segment, "L-PART1"] = str(L_PART1_seq)
            group_locus.loc[index_segment, "L-PART2"] = str(L_PART2_seq)
            group_locus.loc[index_segment, "L-PART"] = str(L_PART1_seq + L_PART2_seq)

            group_locus.loc[index_segment, "DONOR-SPLICE"] = str(DONOR_SPLICE)
            group_locus.loc[index_segment, "ACCEPTOR-SPLICE"] = str(ACCEPTOR_SPLICE)

            group_locus.loc[index_segment, "L-PART+V-EXON"] = str(L_PART1_seq + L_PART2_seq + target_sequence)
            group_locus.loc[index_segment, "Protein"] = str(Seq(L_PART1_seq + L_PART2_seq + target_sequence).translate(to_stop=False).upper())

    return group_locus


def check_functional_protein(segment, l_region, target_sequence, strand, protein, donor_splice, acceptor_splice):

    if segment == "V":
        if len(l_region) % 3 != 0:
            return ["pseudo", "Frameshift in L-part region"]

        if protein[0] != 'M':
            return ["pseudo", "No start codon in protein"]

        if "*" in protein[:-1]:
            return ["pseudo", "Stop codon within protein"]

        if "*" not in protein:
            base_status = ["functional", ""]
        elif protein[-1] == "*":
            base_status = ["ORF", "Last protein is stop codon"]
        else:
            base_status = ["pseudo", ""]

        if donor_splice.upper() != "GT":
            return ["ORF", f"Incorrect donor splice site {donor_splice}"]

        if  acceptor_splice.upper() != "AG":
            return ["ORF", f"Incorrect acceptor splice site {acceptor_splice}"]

        if Seq(target_sequence.replace("-", "")).translate(to_stop=False).upper().count("C") < 2:
            return ["ORF", "Missing Cysteine in protein"]

        return base_status

    elif segment == "J":
        if strand == "-":
            target_sequence = Seq(target_sequence.replace("-", "")).reverse_complement()
        J_target_sequence = Seq(target_sequence)
        motief_count = 0
        for frame in range(3):
            J_protein = str(J_target_sequence[frame:].translate())
            if re.compile(r'[FW]G.G').search(J_protein):
                motief_count += 1
                if "*" in J_protein:
                    return ["pseudo", f"Stop codon in protein, frame {frame}"]

        if motief_count == 0:
            return ["ORF", "No motifs ([FW]G.G) found"]

        if donor_splice.upper() != "GT":
            return ["pseudo", f"In correct donor splice site {donor_splice}"]

        return ["functional", ""]


def save_results(combined_results, cwd):
    known = combined_results[combined_results["Status"] == "Known"]
    novel = combined_results[combined_results["Status"] == "Novel"]

    if not combined_results.empty:
        combined_df = combined_results.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
        combined_df.to_excel(cwd / "annotation" / "annotation_report_all.xlsx", index=False)
    if not known.empty:
        known = known.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
        known.to_excel(cwd / "annotation" / "annotation_report_known.xlsx", index=False)
    if not novel.empty:
        novel = novel.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
        novel.to_excel(cwd / "annotation" / "annotation_report_novel.xlsx", index=False)


def main_functionality(immune_type, species, threads: int = 12) -> None:

    cwd = Path.cwd()
    output_base = cwd / "tmp/functionality"
    output_base.mkdir(parents=True, exist_ok=True)

    library_path = output_base / "library"
    library_path.mkdir(parents=True, exist_ok=True)
    library_file = library_path / "library.fasta"

    if not Path(library_file).is_file():
        imgt_main(species=species, immune_type=immune_type, output_dir=output_base, frame_selection="L-PART1+L-PART2", simple_headers_bool=False)

    write_l_region_fasta_files(library_file, library_path)

    data_V, data_D, data_J = open_files(cwd)
    combined_results = data_D

    #V-REGION
    for region, group_locus in data_V.groupby("Region"):
        locus_fasta_path = output_base / "fasta"
        locus_fasta_path.mkdir(exist_ok=True, parents=True)

        extract_coords = {}
        off_set = 500 - 50
        for l_part, length in zip(["L_PART1", "L_PART2"], [500, 50]):
            locus_fasta_file_name = write_extract_region_from_genome(locus_fasta_path, region, l_part, length, group_locus)
            extract_coords = get_blast_results(output_base, l_part, region, locus_fasta_file_name, library_path, extract_coords)
        group_locus = extraxt_region_from_genome(locus_fasta_path, region, group_locus, extract_coords, off_set)
        combined_results = pd.concat([combined_results, group_locus])

    #J-REGION
    for region, group_locus in data_J.groupby("Region"):
        for index_segment, row in group_locus.iterrows():
            start_coord_segment = row["Start coord"]
            end_coord_segment = row["End coord"]
            strand = row["Strand"]
            fasta_path = row["Path"]

            with open(fasta_path, 'r') as fasta_file:
                sequence_region = SeqIO.read(fasta_file, 'fasta')

            if strand == "+":
                j_downstream_region = sequence_region.seq[end_coord_segment:(end_coord_segment + 2)]
            elif strand == "-":
                j_downstream_region = Seq(sequence_region.seq[(start_coord_segment - 2):start_coord_segment]).reverse_complement()
            group_locus.loc[index_segment, "DONOR-SPLICE"] = str(j_downstream_region)
            group_locus.loc[index_segment, "Protein"] = ""

        combined_results = pd.concat([combined_results, group_locus])

    combined_results[["Function", "Function_messenger"]] = combined_results.apply(lambda row: check_functional_protein(row["Segment"], row["L-PART"], row["Target sequence"], row['Strand'], row["Protein"], row["DONOR-SPLICE"], row["ACCEPTOR-SPLICE"]) if pd.notna(row["Protein"]) else ["pseudo", ""], axis=1, result_type="expand")
    mask = combined_results["Segment"].isin(["D"])
    combined_results.loc[mask, ["Function", "Function_messenger"]] = ["functional", ""]
    save_results(combined_results, cwd)


