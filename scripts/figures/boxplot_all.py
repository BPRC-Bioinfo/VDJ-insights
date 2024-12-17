import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from tqdm import tqdm
import os


def open_files(data_path: str) -> pd.DataFrame:
    data = pd.read_excel(data_path)
    return data

def make_pivot_table(data: pd.DataFrame) -> pd.DataFrame:
    pivot_df = data.pivot_table(
        index=['Sample', 'Haplotype', 'Region'],
        columns='Segment',
        aggfunc='size',
        fill_value=0
    )
    return pivot_df


def make_plot(data: pd.DataFrame, output: str, status: str) -> None:
    data.reset_index(inplace=True)
    long_df = data.melt(id_vars=['Sample', 'Haplotype', 'Region'], var_name='Segment', value_name='Count')
    long_df = long_df[long_df['Count'] > 0]

    segmenten = long_df['Segment'].unique()
    regions = long_df['Region'].unique()

    with tqdm(total=len(long_df[['Segment', 'Region']].drop_duplicates()), desc="Creating plots", unit="Plot") as pbar:
        for segment in segmenten:
            for region in regions:

                data = long_df[(long_df['Segment'] == segment) & (long_df['Region'] == region)]
                if len(data) > 0:
                    plt.figure(figsize=(12, 6), dpi=300)
                    ax = sns.stripplot(
                        data=data,
                        x='Population',
                        y='Count',
                        color='black',
                        #hue="Function",
                        #legend=False,
                    )
                    sns.pointplot(
                        data=data,
                        x='Population',
                        y='Count',
                        dodge=True,
                        errorbar='sd',
                        err_kws={'linewidth': 1},
                        linestyle='none',
                        color='black',
                        marker="_",
                        markersize=10,
                        markeredgewidth=1,
                        capsize = .1
                    )

                    #ax = sns.catplot(data=data, x='Population', y='Count',hue='Region', color='black')
                    #ax = sns.catplot(data=data, x='Sample', y='Count', hue='Population', col='Haplotype')
                    #ax = sns.violinplot(data=data, x='Population', y='Count', hue="Function", split=True)

                    plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

                    sample_counts = data.groupby('Population')['Sample'].nunique()
                    x_labels = [f"{pop} ({sample_counts[pop]})" for pop in sample_counts.index]
                    plt.xticks(ticks=range(len(x_labels)), labels=x_labels, rotation=90)

                    plt.title(f'Distribution of {segment} segment counts in {region} region across populations (N={sample_counts.sum()})')
                    plt.xlabel('Population')
                    plt.ylabel('Segment count')
                    plt.xticks(rotation=90)

                    plt.tight_layout()

                    plt.savefig(f"{output}/{region}_{segment}_{status}.pdf", dpi=300)
                    plt.savefig(f'{output}/{region}_{segment}_{status}.svg', dpi=300)
                    pbar.update(1)


def main(path) -> None:
    output = f"{path}/figure/boxplot_all"

    os.makedirs( output, exist_ok=True)

    data = open_files(data_path=f"{path}/annotation/annotation_report_known_rss.xlsx")

    pivot_df = make_pivot_table(data=data)
    make_plot(data=pivot_df, output=output)

if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr"
    main(path)
