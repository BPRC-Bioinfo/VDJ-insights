import os

import pandas as pd
import matplotlib.pyplot as plt
import re

import numpy as np


def open_files(data_path: str) -> pd.DataFrame:
    data = pd.read_excel(data_path)
    return data


def make_pivot_table(data: pd.DataFrame) -> pd.DataFrame:
    pivot_df = data.pivot_table(
        index=['Population','Sample', 'Status', 'Region'],
        columns='Short name',
        aggfunc='size',
        fill_value=0,
    ).reset_index()
    return pivot_df


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]


def create_barplot(data: pd.DataFrame, region: str, output: str, status : str) -> None:
    filtered_data = data[(data['Region'] == region) & (data['Status'] == status.capitalize())]
    reference_labels = filtered_data.columns[6:]
    population_sums = filtered_data.groupby('Population')[reference_labels].sum()
    non_zero_columns = population_sums.loc[:, (population_sums != 0).any(axis=0)]

    if non_zero_columns.empty:
        return

    num_references = len(non_zero_columns.columns)
    populations = non_zero_columns.index
    fig, ax = plt.subplots(figsize=(num_references * 0.35, 8))

    bottom = np.zeros(num_references)

    num_colors = len(populations)
    colors = [plt.cm.tab20(i / 20) for i in range(20)]
    if num_colors > 20:
        extra_colors = [plt.cm.viridis(i / (num_colors - 20)) for i in range(num_colors - 20)]
        colors.extend(extra_colors)
    population_colors = dict(zip(populations, colors[::-1]))

    for population in populations:
        values = non_zero_columns.loc[population].values
        ax.bar(range(num_references), values, bottom=bottom, color=population_colors[population], label=population)
        bottom += values

    for i, total in enumerate(bottom):
        ax.text(i, total, int(total), ha='center', va='bottom', fontsize=13)

    ax.set_xticks(range(num_references))
    ax.set_xticklabels(non_zero_columns.columns, rotation=90, fontsize=13, ha='center')
    ax.set_xlim(-0.5, num_references - 0.5)

    ax.set_title(f"Population distribution {status} segments for region: {region}", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_xlabel('Allele', fontsize=14)


    #ax.legend(title="Population", ncol=len(populations), fontsize=10, frameon=False)

    ax.legend(
        title="Population",
        bbox_to_anchor=(0.5, -0.50),
        loc='upper center',
        ncol=min(len(populations), 10),
        fontsize=10,
        frameon=False
    )
    plt.tight_layout(rect=[0, 0.2, 1, 1])
    os.makedirs(f"{output}", exist_ok=True)

    plt.savefig(f"{output}/{region}_{status}_bar.pdf", dpi=300)
    plt.savefig(f"{output}/{region}_{status}_bar.svg", dpi=300)
    plt.close(fig)


def main(path) -> None:
    output = f"{path}/figure/barplot"
    os.makedirs(output, exist_ok=True)

    for status in ['known', 'novel']:
        data = open_files(data_path=f"{path}/annotation_report_{status}.xlsx")

        pivot_df = make_pivot_table(data)

        regions = pivot_df['Region'].unique()

        for region in regions:
            create_barplot(pivot_df, region=region, output=output, status=status)


if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr/annotation"
    main(path)

