import json
import os
import re
import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


def open_files(cwd: Path) -> (pd.DataFrame, pd.core.groupby.generic.DataFrameGroupBy):
    file_known = cwd / "annotation" / "annotation_report_known_rss.xlsx"
    file_novel = cwd / "annotation" / "annotation_report_novel_rss.xlsx"

    dataframes = []
    if file_known.exists():
        dataframes.append(pd.read_excel(file_known))
    if file_novel.exists():
        dataframes.append(pd.read_excel(file_novel))
    data = pd.concat(dataframes, ignore_index=True)

    data_other = data[~data["Segment"].isin(["V"])]
    data_v = data[data["Segment"].isin(["V"])].copy()

    data_v["Short name strip"] = data_v["Short name"].str.split(r'[-*]').str[0]
    v_grouped = data_v.groupby("Short name strip")
    return data_other, v_grouped


def scrape_imgt(species: str, path: Path) -> dict:
    segments = {
        "TR": ["TRBV", "TRAV", "TRDV", "TRGV"],
        "IG": ["IGHV", "IGKV", "IGLV"]
    }
    all_sequences = {}
    for immune_type, seg_list in segments.items():
        for segment in seg_list:
            url = f"https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/{species}/{immune_type}/{segment}.fasta"
            response = requests.get(url)
            if response.status_code == 200:
                fasta_io = StringIO(response.text)
                sequences = SeqIO.to_dict(SeqIO.parse(fasta_io, "fasta"))
                all_sequences.update(sequences)
            else:
                print(f"Waarschuwing: kon {url} niet ophalen (status code {response.status_code}).")
        combined_file = path / f"{species}_{immune_type}_combined.fasta"
        os.makedirs(combined_file.parent, exist_ok=True)
        with open(combined_file, 'w') as output_handle:
            SeqIO.write(all_sequences.values(), output_handle, "fasta")
    return all_sequences


def create_dataframe(all_sequences: dict) -> pd.DataFrame:
    data = []
    for seq_id, record in all_sequences.items():
        try:
            gene_info = seq_id.split('|')[1]
            gene_name = re.split(r'[*]', gene_info)[0]
            gene_subgroup = re.split(r'[-*]', gene_info)[0]

            cdr1_sequence = record.seq[78:114].upper().replace(".", "")
            cdr2_sequence = record.seq[165:195].upper().replace(".", "")

            description = record.description.split("|")
            allele = description[1] if len(description) > 1 else ""

            data.append({
                "Gene_subgroup": gene_subgroup,
                "Gene_name": gene_name,
                "Allele": allele,
                "CDR1_sequence": str(cdr1_sequence),
                "CDR2_sequence": str(cdr2_sequence),
            })
        except IndexError:
            print(f"IndexError voor sequentie: {seq_id}. Sla deze over.")
    return pd.DataFrame(data)


def run_meme_analysis(fasta_path: Path, output_dir: Path) -> None:
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    if not sequences:
        return

    lengths = [len(record.seq) for record in sequences]
    minw, maxw = min(lengths), max(lengths)
    multi_command = f"meme {fasta_path} -o {output_dir} -dna -mod zoops -minw {minw} -maxw {maxw}"

    amount = subprocess.run(f"grep -c '^>' {fasta_path}", shell=True, capture_output=True, text=True)
    try:
        if int(amount.stdout.strip()) > 1 and not output_dir.exists():
            subprocess.run(multi_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"{e}")


def run_fimo_analysis(meme_output: Path, fasta_path: Path, output_dir: Path) -> None:
    command = f"fimo --o {output_dir} {meme_output / 'meme.txt'} {fasta_path}"
    amount = subprocess.run(f"grep -c '^>' {fasta_path}", shell=True, capture_output=True, text=True)
    try:
        if int(amount.stdout.strip()) > 1 and not output_dir.exists():
            subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"{e}")


def get_fimo_output(fimo_input: Path) -> pd.DataFrame:
    if fimo_input.is_dir():
        fimo_output_file = fimo_input / "fimo.txt"
        if fimo_output_file.exists():
            df_fimo = pd.read_csv(fimo_output_file, sep='\t')
            df_fimo["index_group_df"] = df_fimo["sequence name"].str.split("_").str[-1]
            df_fimo['index_group_df'] = df_fimo['index_group_df'].astype(int)
            return df_fimo
    return pd.DataFrame()


def process_variant(gene_subgroup: str, group_df: pd.DataFrame, all_sequences_df: pd.DataFrame, output_base: Path,  cwd: Path) -> pd.DataFrame:
    combined_results = pd.DataFrame()

    f_df = all_sequences_df[all_sequences_df["Gene_subgroup"] == gene_subgroup].copy()
    if not f_df.empty:
        safe_gene_subgroup = gene_subgroup.replace("/", "_")
        for cdr_region in [1, 2]:
            cdr_fasta_path = output_base / "fasta" / "library_split" / f"{safe_gene_subgroup}_CDR{cdr_region}.fasta"
            os.makedirs(output_base / "fasta" / "library_split", exist_ok=True)
            with open(cdr_fasta_path, 'w') as out_handle:
                for _, row in f_df.iterrows():
                    out_handle.write(f">{row['Allele']}\n{row[f'CDR{cdr_region}_sequence']}\n")

            cdr_meme_output = output_base / "meme_output" / f"{safe_gene_subgroup}_CDR{cdr_region}"
            os.makedirs(output_base / "meme_output", exist_ok=True)
            run_meme_analysis(fasta_path=cdr_fasta_path, output_dir=cdr_meme_output)

            fasta_known_path = output_base / "fasta" / "samples_fasta" / f"{safe_gene_subgroup}.fasta"
            os.makedirs(output_base / "fasta" / "samples_fasta", exist_ok=True)
            with open(fasta_known_path, 'w') as out_handle:
                for index_segment, row in group_df.iterrows():
                    out_handle.write(f">{row['Target name']}_{index_segment}\n{row['Target sequence']}\n")

            fimo_output = output_base / "fimo_output" / f"{safe_gene_subgroup}_CDR{cdr_region}"
            os.makedirs(output_base / "fimo_output", exist_ok=True)
            run_fimo_analysis(meme_output=cdr_meme_output, fasta_path=fasta_known_path, output_dir=fimo_output)
            df_fimo = get_fimo_output(fimo_input=fimo_output)
            if not df_fimo.empty:
                for index_segment, row in group_df.iterrows():
                    df_match = df_fimo[df_fimo["index_group_df"] == index_segment]
                    if not df_match.empty and df_match["score"].values[0] > 0:
                        group_df.loc[index_segment, f"CDR_{cdr_region}_start"] = df_match["start"].values[0]
                        group_df.loc[index_segment, f"CDR_{cdr_region}_stop"] = df_match["stop"].values[0]
                        matched_seq = df_match["matched sequence"].values[0].upper()
                        aa_seq = str(Seq(matched_seq).translate())
                        group_df.loc[index_segment, f"CDR_{cdr_region}_aa"] = aa_seq
                        group_df.loc[index_segment, f"CDR_{cdr_region}_seq"] = matched_seq

        combined_results = pd.concat([combined_results, group_df])
    return combined_results


def main_cdr(threads: int = 12) -> None:
    """
    Hoofdfunctie om de CDR-annotaties te verwerken.
    """
    cwd = Path.cwd()
    output_base = cwd / "CDR"

    library_path = output_base / "fasta" / "library"
    library_path.mkdir(parents=True, exist_ok=True)

    all_sequences = scrape_imgt(species="Homo_sapiens", path=library_path)
    if all_sequences:
        all_sequences_df = create_dataframe(all_sequences=all_sequences)

        data_other, v_grouped = open_files(cwd)
        combined_df = data_other.copy()

        max_jobs = 2
        with ProcessPoolExecutor(max_workers=max_jobs) as executor:
            futures = [
                executor.submit(process_variant, gene_subgroup, group_locus, all_sequences_df, output_base, cwd)
                for gene_subgroup, group_locus in v_grouped
            ]
            with tqdm(total=v_grouped.ngroups, desc="Processing CDR", unit='task') as pbar:
                for future in as_completed(futures):
                    result = future.result()
                    combined_df = pd.concat([combined_df, result])
                    pbar.update(1)

        known = combined_df[combined_df["Status"] == "Known"]
        novel = combined_df[combined_df["Status"] == "Novel"]

        if not combined_df.empty:
            combined_df = combined_df.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
            combined_df.to_excel(cwd / "annotation" / "annotation_report_all_rss_cdr.xlsx", index=False)
        if not known.empty:
            known = known.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
            known.to_excel(cwd / "annotation" / "annotation_report_known_rss_cdr.xlsx", index=False)
        if not novel.empty:
            novel = novel.sort_values(by=['Sample', 'Region', 'Start coord'], ascending=[True, True, True])
            novel.to_excel(cwd / "annotation" / "annotation_report_novel_rss_cdr.xlsx", index=False)


if __name__ == '__main__':
    main_cdr(threads=12)
