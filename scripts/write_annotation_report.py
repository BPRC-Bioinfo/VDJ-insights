from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq


def main_df(df):
    mask = ~df['query seq'].str.contains(
        '-') & ~df['subject seq'].str.contains('-')
    df = df[mask]
    df['% identity'] = df['% identity'].astype(float)
    reference_df = df.query("`% identity` == 100.000")
    comperere_df = df.query("`% identity` != 100.000")
    merged_df = pd.merge(comperere_df, reference_df[[
                         'query', 'subject']], on='query', suffixes=('', '_100'))
    return merged_df


def add_change_values_df(df):
    df['% Mismatches of total alignment'] = (
        df['mismatches'] / df['alignment length']) * 100
    df['query_seq_length'] = df['query seq'].str.len()
    df['subject_seq_length'] = df['subject seq'].str.len()
    path_df = df['query'].str.split(':', expand=True)
    df[['query', 'start', 'stop', 'strand', 'path']] = path_df[[0, 1, 2, 3, 4]]
    output_df = df[[
        'subject_100', 'query',
        'mismatches', '% Mismatches of total alignment',
        'start', 'stop',
        'subject seq', 'query seq',
        'strand', 'path'
    ]]
    output_df.columns = [
        'Reference', 'Old name-like',
        'Mismatches', '% Mismatches of total alignment',
        'Start coord', 'End coord',
        'Reference seq', 'Old name-like seq',
        'Strand', 'Path'
    ]
    output_df = output_df.sort_values(by="Reference")
    output_df['Old name-like'] = output_df['Old name-like'] + '-like'
    output_df['Count'] = output_df.groupby(
        ['Reference', 'Old name-like']).cumcount() + 1
    output_df['Modified'] = output_df['Old name-like'] + \
        '-' + output_df['Count'].astype(str)
    output_df['Old name-like'] = output_df['Modified']
    output_df.drop(['Modified', 'Count'], axis=1, inplace=True)
    return output_df


def annotation(df, annotation_folder):
    df = df[['Reference', 'Old name-like', 'Mismatches',
             '% Mismatches of total alignment', 'Start coord',
             'End coord', 'Function', 'Path']]
    df.to_excel(annotation_folder / 'annotation_report.xlsx', index=False)


def rss(df, annotation_folder):
    df = df[['Old name-like', 'Start coord',
             'End coord', 'Strand', 'Path']]
    df.to_excel(annotation_folder / 'RSS_report.xlsx', index=False)


def add_orf(row):
    sequence, strand = row[['Old name-like seq', 'Strand']]
    sequence = Seq(sequence)
    if strand == "-":
        sequence = sequence.reverse_complement()
    aa = sequence.translate()
    row["Function"] = "F/ORF" if "*" not in aa else "P"
    return row


def write_report_report(annotation_folder):
    df = pd.read_excel(annotation_folder / "blast_results.xlsx")
    df = main_df(df)
    df = add_change_values_df(df)
    df = df.apply(add_orf, axis=1)
    annotation(df, annotation_folder)
    rss(df, annotation_folder)


if __name__ == '__main__':
    cwd = Path.cwd()
    annotation_folder = cwd / "annotation"
    write_report_report(annotation_folder)
