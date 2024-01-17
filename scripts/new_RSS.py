import re
import json
import subprocess
import tempfile
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq


def write_fasta_file(dictionary, folder):
    Path(folder).mkdir(parents=True, exist_ok=True)
    for key, value in dictionary.items():
        for region, segments in value.items():
            ffile = folder / f"{key}{region}.fasta"
            with open(ffile, 'w') as f:
                for header, segment in segments.items():
                    for x, rss in enumerate(segment):
                        f.write(f">{header}-{x+1}\n{rss}\n")


def rss_type(start, end, segment, strand):
    variant_23 = 39
    variant_12 = 28

    TR_seg = {
        "V": {"+": [end, str(int(end) + variant_23)],
              "-": [str(int(start) - variant_23), start]},
        "D": {"+": [str(int(start) - variant_12), str(int(end) + variant_12)],
              "-": [str(int(start) - variant_12), str(int(end) + variant_12)]},
        "J": {"+": [str(int(start) - variant_12), start],
              "-": [end, str(int(end) + variant_12)]}
    }
    return TR_seg[segment][strand]


def fetch_sequence(row):
    query, fasta, start, end, strand, region, segment = row[[
        'Old name-like', 'Path', 'Start coord', 'End coord', 'Strand',
        'Region', 'Segment']]
    with open(fasta, 'r') as fasta_file:
        record = SeqIO.read(fasta_file, 'fasta')
        rss_start, rss_stop = rss_type(start, end, segment, strand)
        rss = record.seq[int(rss_start):int(rss_stop)]
        if strand == '-':
            rss = str(Seq(str(rss)).reverse_complement())
    return query, rss


def add_to_dict(query, dictionary, rss):
    dictionary.setdefault(query, []).append(rss)


def get_mers(segment, rss):
    # ADD D SEGMENT
    mers = {
        "V": [rss[0:7], rss[-9:]],
        "D": [],
        "J": [rss[0:9], rss[-7:]]
    }
    return mers[segment]


def load_header(row, separated_segments):
    region, segment, function = row[["Region", "Segment", "Function"]]
    query, rss_sequence = fetch_sequence(row)
    val1, val2 = get_mers(segment, rss_sequence)
    # ADD D SEGMENT
    if len(val1) == 7:
        row["heptamer"], row["nonamer"] = ''.join(
            val1), ''.join(val2)
    else:
        row["nonamer"], row["heptamer"] = ''.join(
            val1), ''.join(val2)
    if function != "P":
        add_to_dict(query, separated_segments.setdefault(
            region, {}).setdefault(segment, {}), f"{str(rss_sequence)}")
    return row


def run_meme(out, rss_file: Path):
    command = f"meme {rss_file} -o {out} -dna -mod zoops -nmotifs 1 -minw 6 -maxw 50"
    subprocess.run(command, shell=True)


def get_reference_mers(pattern, segment):
    elements = re.findall(r'\[[^\]]*\]|.', pattern)
    ref_mers = {
        "V": [elements[:7], elements[-9:]],
        "D": [elements[:7], elements[-9:]],
        "J": [elements[:9], elements[-7:]]
    }
    val1, val2 = ref_mers[segment]
    return val1, val2


def make_referece_rss(cwd):
    ref_rss = {}
    Path(cwd / "meme").mkdir(parents=True, exist_ok=True)
    rss = cwd / "RSS"
    for rss_file in rss.iterdir():
        stem = rss_file.stem
        out = cwd / "meme" / stem
        meme = out / "meme.txt"
        if not meme.exists():
            run_meme(out, rss_file)
        command = f'cat {meme} | egrep -A2 "regular expression"'
        result = subprocess.run(command, shell=True,
                                capture_output=True, text=True)
        hits = result.stdout.replace(
            "-", "").replace("\t", "").strip().split("\n")
        hits = [hit for hit in hits if hit]
        for i in range(0, len(hits), 2):
            val1, val2 = get_reference_mers(hits[i+1], stem[-1])
            if len(val1) == 7:
                ref_rss.setdefault(stem, {}).setdefault(
                    "heptamer", ''.join(val1))
                ref_rss.setdefault(stem, {}).setdefault(
                    "nonamer", ''.join(val2))
            else:
                ref_rss.setdefault(stem, {}).setdefault(
                    "heptamer", ''.join(val2))
                ref_rss.setdefault(stem, {}).setdefault(
                    "nonamer", ''.join(val1))

    return ref_rss


def check_ref_rss(row, ref_rss):
    region, segment = row[["Region", "Segment"]]
    ref = ref_rss[region+segment]
    for i in ["heptamer", "nonamer"]:
        ref_seq = ref[i]
        query_seq = row[i]
        matches = re.findall(ref_seq, query_seq)
        row[f"ref_{i}"] = ref_seq
        if matches:
            row[f"{i}_matched"] = True
        else:
            row[f"{i}_matched"] = False
    return row


def add_region_segment(row):
    query = row['Old name-like']
    prefix = [i for i in query.split("_") if i.startswith(("TR", "LOC"))][0]
    prefix = re.sub(r"[0-9-]", "", prefix)
    region, segment = prefix[0:3], prefix[3]
    row["Region"], row["Segment"] = region, segment
    return row


def combine_df(original_df, new_df):
    columns_to_merge = [
        'Reference', 'Old name-like',
        'heptamer', 'ref_heptamer', 'heptamer_matched',
        'nonamer', 'ref_nonamer', 'nonamer_matched'
    ]
    combined_df = pd.merge(new_df, original_df[columns_to_merge],
                           on=['Reference', 'Old name-like'],
                           how='left')
    return combined_df


def main():
    cwd = Path.cwd()
    df = pd.read_excel(cwd / 'annotation' / 'RSS_report.xlsx')
    df = df.apply(add_region_segment, axis=1)
    separated_segments = {}
    df = df.apply(lambda row: load_header(row, separated_segments), axis=1)
    original = df.copy()
    df = df.query("Function != 'P'").reset_index(drop=True)
    rss_path = cwd / 'RSS'
    if not rss_path.exists():
        write_fasta_file(separated_segments, rss_path)
    ref_rss = make_referece_rss(cwd)
    original = original.apply(lambda row: check_ref_rss(row, ref_rss), axis=1)
    final_df = combine_df(original, pd.read_excel(
        cwd / 'annotation' / 'annotation_report.xlsx'))
    final_df.to_excel(cwd / 'check.xlsx', index=False)


if __name__ == '__main__':
    main()
