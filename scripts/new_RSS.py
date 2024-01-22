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


def rss_type(start, end, segment, rss_variant, strand):
    variant_12 = 28
    variant_23 = 39
    TR_seg = {
        "V": {23: {"+": [end, str(int(end) + variant_23)],
                   "-": [str(int(start) - variant_23), start]}},
        "D": {12: {"+": [str(int(start) - variant_12), start],
                   "-": [end, str(int(end) + variant_12)]},
              23: {"+": [end, str(int(end) + variant_23)],
                   "-": [str(int(start) - variant_23), start]}},
        "J": {12: {"+": [str(int(start) - variant_12), start],
                   "-": [end, str(int(end) + variant_12)]}}
    }
    return TR_seg[segment][rss_variant][strand]


def fetch_sequence(row, segment, rss_variant):
    query, fasta, start, end, strand, region = row[[
        'Old name-like', 'Path', 'Start coord', 'End coord', 'Strand',
        'Region']]
    with open(fasta, 'r') as fasta_file:
        record = SeqIO.read(fasta_file, 'fasta')
        rss_start, rss_stop = rss_type(start, end, segment, rss_variant, strand)
        rss = record.seq[int(rss_start):int(rss_stop)]
        if strand == '-':
            rss = str(Seq(str(rss)).reverse_complement())
    return query, rss


def add_to_dict(query, dictionary, rss):
    dictionary.setdefault(query, []).append(rss)


def get_mers(segment, rss, rss_variant):
    mers = {
        "V": {23: [rss[0:7], rss[-9:]]},
        "D": {12: [rss[0:9], rss[-7:]], 23: [rss[0:9], rss[-9:]]},
        "J": {12: [rss[0:9], rss[-7:]]}
    }
    return mers[segment].get(rss_variant, ["", ""])


def add_to_row(row, val1, val2, rss_type):
    heptamer, nonamer = (''.join(val1), ''.join(val2)) if len(
        val1) == 7 else (''.join(val2), ''.join(val1))
    row[f"{rss_type}_heptamer"], row[f"{rss_type}_nonamer"] = heptamer, nonamer
    return row


def add_segment(query, segment, region, separated_segments, rss_sequence, function):
    if function != "P":
        add_to_dict(query, separated_segments.setdefault(
            region, {}).setdefault(segment, {}), str(rss_sequence))


def load_header(row, separated_segments):
    region, segment, function = row[["Region", "Segment", "Function"]]
    # Initialize RSS types for D segment
    row["12_heptamer"], row["12_nonamer"], row["23_heptamer"], row["23_nonamer"] = "", "", "", ""

    if segment == "D":
        for rss_variant in [12, 23]:
            query, rss_sequence = fetch_sequence(row, segment, rss_variant)
            val1, val2 = get_mers(segment, rss_sequence, rss_variant)
            add_to_row(row, val1, val2, rss_variant)
            add_segment(query, f"{segment}_{rss_variant}", region,
                        separated_segments, rss_sequence, function)
    else:
        rss_variant = 23 if segment == "V" else 12
        query, rss_sequence = fetch_sequence(row, segment, rss_variant)
        val1, val2 = get_mers(segment, rss_sequence, rss_variant)
        add_to_row(row, val1, val2, rss_variant)
        add_segment(query, f"{segment}_{rss_variant}", region,
                    separated_segments, rss_sequence, function)

    return row


def run_meme(out, rss_file: Path):
    command = f"meme {rss_file} -o {out} -dna -mod zoops -nmotifs 1 -minw 6 -maxw 50"
    subprocess.run(command, shell=True)


def get_reference_mers(pattern, segment, rss_variant):
    rss = re.findall(r'\[[^\]]*\]|.', pattern)
    ref_mers = {
        "V": {23: [rss[0:7], rss[-9:]]},
        "D": {12: [rss[0:9], rss[-7:]],
              23: [rss[0:9], rss[-9:]]},
        "J": {12: [rss[0:9], rss[-7:]]}
    }
    val1, val2 = ref_mers[segment][rss_variant]
    return val1, val2


def make_ref_dict(segment, ref_rss, val1, val2):
    if len(val1) == 7:
        ref_rss.setdefault(segment, {}).setdefault(
            "heptamer", ''.join(val1))
        ref_rss.setdefault(segment, {}).setdefault(
            "nonamer", ''.join(val2))
    else:
        ref_rss.setdefault(segment, {}).setdefault(
            "heptamer", ''.join(val2))
        ref_rss.setdefault(segment, {}).setdefault(
            "nonamer", ''.join(val1))
    return ref_rss


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
        split_stem = stem.split("_")
        segment, rss_variant = split_stem[0][-1], int(split_stem[1])
        val1, val2 = get_reference_mers(hits[1], segment, rss_variant)
        make_ref_dict(stem, ref_rss, val1, val2)
    return ref_rss


def check_ref_rss(row, ref_rss, rss_variant):
    region, segment = row[["Region", "Segment"]]
    ref = ref_rss[f"{region}{segment}_{rss_variant}"]
    for i in ["heptamer", "nonamer"]:
        ref_seq = ref[i]
        query_seq = row[f"{rss_variant}_{i}"]
        matches = re.findall(ref_seq, query_seq)
        row[f"{rss_variant}_ref_{i}"] = ref_seq
        if matches:
            row[f"{rss_variant}_{i}_matched"] = True
        else:
            row[f"{rss_variant}_{i}_matched"] = False
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
        '12_heptamer', '12_heptamer_matched', '12_nonamer',
        '12_nonamer_matched', '12_ref_heptamer', '12_ref_nonamer',
        '23_heptamer', '23_heptamer_matched', '23_nonamer',
        '23_nonamer_matched', '23_ref_heptamer', '23_ref_nonamer'
    ]
    combined_df = pd.merge(new_df, original_df[columns_to_merge],
                           on=['Reference', 'Old name-like'],
                           how='left')
    return combined_df


def apply_check_ref_rss(row, ref_rss, rss_variant):
    segment_last_char = row["Segment"][-1]
    if segment_last_char == "D":
        return check_ref_rss(row, ref_rss, rss_variant)
    if segment_last_char == "V" and rss_variant == 23:
        return check_ref_rss(row, ref_rss, rss_variant)
    if segment_last_char == "J" and rss_variant == 12:
        return check_ref_rss(row, ref_rss, rss_variant)
    return row


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
    for rss_variant in [12, 23]:
        original = original.apply(
            apply_check_ref_rss, args=(ref_rss, rss_variant), axis=1)
    final_df = combine_df(original, pd.read_excel(
        cwd / 'annotation' / 'annotation_report.xlsx'))
    final_df.to_excel(cwd / 'annotation' /
                      'annotation_report_plus.xlsx', index=False)


if __name__ == '__main__':
    main()
