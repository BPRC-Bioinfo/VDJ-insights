import re
import json
import subprocess
import tempfile
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq


def rss_type(start, end, segment, strand):
    TR_seg = {
        "V": {"+": [end, str(int(end) + 39)], "-": [str(int(start) - 39), start]},
        "D": {"+": [str(int(start) - 28), str(int(end) + 28)], "-": [str(int(start) - 28), str(int(end) + 28)]},
        "J": {"+": [str(int(start) - 39), start], "-": [end, str(int(end) + 39)]}
    }
    return TR_seg[segment][strand]


def fetch_sequence(row, ):
    query, fasta, start, end, strand = row[[
        'Old name-like', 'Path', 'Start coord', 'End coord', 'Strand']]
    with open(fasta, 'r') as fasta_file:
        record = SeqIO.read(fasta_file, 'fasta')
        header = record.id
        prefix = [i for i in query.split("_") if i.startswith(("TR", "LOC"))][0]
        prefix = re.sub(r"[0-9-]", "", prefix)
        region, segment = prefix[0:3], prefix[3]
        rss_start, rss_stop = rss_type(start, end, segment, strand)
        rss = record.seq[int(rss_start):int(rss_stop)]
        if strand == '-':
            rss = str(Seq(str(rss)).reverse_complement())
    return region, segment, query, rss


def add_to_dict(query, dictionary, rss):
    dictionary.setdefault(query, []).append(rss)


def load_header(row, separated_segments):
    region, segment, query, rss_sequence = fetch_sequence(
        row)
    add_to_dict(query, separated_segments.setdefault(
        region, {}).setdefault(segment, {}), str(rss_sequence))


def main():
    cwd = Path.cwd()
    df = pd.read_excel(cwd / 'annotation' / 'RSS_report.xlsx')

    separated_segments = {}
    df.apply(lambda row: load_header(row, separated_segments), axis=1)
    print(json.dumps(separated_segments, indent=4))


if __name__ == '__main__':
    main()
