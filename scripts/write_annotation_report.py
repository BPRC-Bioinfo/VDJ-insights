from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def make_plot(path, value_counts, type):
    plt.figure(figsize=(15, 6))
    value_counts.plot(kind=type)
    plt.xlabel('Mismatches')
    plt.ylabel('Count')
    plt.title('Count of Mismatches')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(path)


def main_df(df):
    mask = ~df['query seq'].str.contains(
        '-') & ~df['subject seq'].str.contains('-')
    df = df[mask]
    # Ensure that '% identity' is a float
    df['% identity'] = df['% identity'].astype(float)

    # Create reference and comparison DataFrames
    reference_df = df.query("`% identity` == 100.000")
    comperere_df = df.query("`% identity` != 100.000")

    # Merge the DataFrames
    merged_df = pd.merge(comperere_df, reference_df[[
                         'query', 'subject']], on='query', suffixes=('', '_100'))
    return merged_df


def add_change_values_df(df,):
    # Calculate % Mismatches of Total Alignment
    df['% Mismatches of total alignment'] = (
        df['mismatches'] / df['alignment length']) * 100

    # Calculate lengths of the sequences
    df['query_seq_length'] = df['query seq'].str.len()
    df['subject_seq_length'] = df['subject seq'].str.len()
    path_df = df['query'].str.split(':', expand=True)
    df['query'] = path_df[0]
    df['start'] = path_df[1]
    df['stop'] = path_df[2]
    df['path'] = path_df[3]
    # Select and rename required columns
    output_df = df[[
        'subject_100', 'query',
        'mismatches', '% Mismatches of total alignment',
        'start', 'stop',
        'path'
    ]]
    output_df.columns = [
        'Reference seq', 'Old name-like',
        'Mismatches', '% Mismatches of total alignment',
        'Start coord', 'End coord',
        'Path'
    ]
    output_df = output_df.sort_values(by="Reference seq")
    output_df['Old name-like'] = output_df['Old name-like'] + '-like'
    return output_df


def main():
    cwd = Path.cwd()
    df = pd.read_excel(cwd / "demo" / "blast_results.xlsx")
    df = main_df(df)
    df = add_change_values_df(df)
    df.to_excel(cwd / 'annotation_report.xlsx', index=False)


if __name__ == '__main__':
    main()
