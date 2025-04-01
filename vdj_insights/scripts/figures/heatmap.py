import os

import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
import re

from tqdm import tqdm
import json
from matplotlib import gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def open_files(data_path: str) -> pd.DataFrame:
    data = pd.read_excel(data_path)
    return data


def open_fasta(data_path: str) -> list:
    library = []
    with open(data_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                allel = line.strip().split("_")[1]
                library.append(allel)
    return library


def make_pivot_table(data: pd.DataFrame) -> pd.DataFrame:
    pivot_df = data.pivot_table(
        index=['Population', 'Sample', 'Haplotype', 'Function', 'Region'],
        #index=['Population', 'Sample','accession name', 'Haplotype_meta', 'Region'],
        columns='Short name',
        aggfunc='size',
        fill_value=0,
    ).reset_index()
    #pivot_df = pivot_df.sort_index(axis=1, key=lambda x: x.map(lambda col: libary.index(col) if col in libary else float('inf')))
    #pivot_df.to_csv("output/pivot_genes_samples.csv", sep=',', index=False)
    return pivot_df


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]


def create_matrix(data: pd.DataFrame, function: str, region: str, output: str,  rgb_color: tuple = (0.57, 0.5, 0.120)) -> None:
    filtered_data = data[(data['Function'] == function) & (data['Region'] == region)]
    #filtered_data = data[(data['Region'] == region)]

    #filtered_data.loc[:, 'Sample'] = filtered_data['Sample'].astype(str) + "_" + filtered_data['Sample'].astype(str)

    filtered_data = filtered_data.sort_values(
        by=['Population', 'Sample', 'Haplotype'],
        ascending=[True, True, False],
        key=lambda x: x.map(filtered_data['Population'].value_counts()) if x.name == 'Population' else x
    )
    if filtered_data.empty:
        return
    references = filtered_data.columns[6:]
    presence_absence_matrix = (filtered_data[references] > 0).astype(int)

    non_empty_references = presence_absence_matrix.any(axis=0)
    presence_absence_matrix = presence_absence_matrix.loc[:, non_empty_references]

    count_matrix = filtered_data[references].loc[:, non_empty_references].values
    filtered_data = filtered_data.loc[:, filtered_data.columns[:6]].join(presence_absence_matrix)

    cm = presence_absence_matrix.values

    samples = filtered_data['Sample'].values
    reference_labels = presence_absence_matrix.columns.values
    if cm.size == 0 or len(reference_labels) == 0:
        return

    num_samples = len(samples)
    num_references = len(reference_labels)

    fig = plt.figure(figsize=(num_references * 0.35 + 2, num_samples * 0.35 + 2))
    gs = gridspec.GridSpec(2, 2, width_ratios=[5, 0.5], height_ratios=[0.5, 5], hspace=0.05, wspace=0.05)
    unique_populations = filtered_data['Population'].unique()

    num_colors = len(unique_populations)
    colors = [plt.cm.tab20(i / 20) for i in range(20)]
    if num_colors > 20:
        extra_colors = [plt.cm.viridis(i / (num_colors - 20)) for i in range(num_colors - 20)]
        colors.extend(extra_colors)
    population_colors = dict(zip(unique_populations, colors[::-1]))

    #bovenste bar
    ax_bar_x = plt.subplot(gs[0, 0], sharex=plt.subplot(gs[1, 0]))
    population_sums = filtered_data.groupby('Population')[reference_labels].sum()
    populations = population_sums.index
    bottom = np.zeros(len(reference_labels))
    for population in populations:
        values = population_sums.loc[population].values
        bars_x = ax_bar_x.bar(range(num_references), values, bottom=bottom, color=population_colors[population], label=population)
        bottom += values
    for i, total in enumerate(bottom):
        ax_bar_x.text(i, total, int(total), ha='center', va='bottom', fontsize=14)
    ax_bar_x.set_xlim(-0.5, num_references - 0.5)
    ax_bar_x.axis('off')

    handles, labels = ax_bar_x.get_legend_handles_labels()
    ax_bar_x.legend(handles, labels, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",prop = { "size": 10 }, ncol=len(labels)).set_title("Population")
    ax_bar_x.set_title(f'Region: {region} (N={num_samples})', fontsize=20, y=1.3, weight='bold', loc='left')

    #zijkant bar
    #allele_names = list(set([col.split('-')[0] for col in reference_labels]))
    allele_names = list(set([re.split(r'[-*]', col)[0] for col in reference_labels]))

    allele_names = sorted(allele_names, key=natural_sort_key)

    allele_sums = pd.DataFrame(count_matrix, columns=reference_labels)
    bottom = np.zeros(len(samples))
    ax_bar_y = plt.subplot(gs[1, 1], sharey=plt.subplot(gs[1, 0]))
    for allele_name in allele_names:
        #matching_columns = [col for col in reference_labels if col.split('-')[0] == allele_name]
        matching_columns = [col for col in reference_labels if re.split(r'[-*]', col)[0] == allele_name]

        values = allele_sums[matching_columns].sum(axis=1)
        bars_y = ax_bar_y.barh(range(num_samples), values, left=bottom, color=plt.cm.tab20.colors[hash(allele_name) % 12], label=allele_name)
        bottom += values
    for i, total in enumerate(bottom):
        ax_bar_y.text(total, i, int(total), ha='left', va='center', fontsize=14)
    ax_bar_y.set_ylim(-0.5, num_samples - 0.5)
    ax_bar_y.axis('off')

    ax_legend = plt.subplot(gs[0, 1])
    ax_legend.axis('off')
    handles_right = [plt.Line2D([0], [0], color=plt.cm.tab20.colors[hash(allele_name) % 12], lw=4) for allele_name in allele_names]
    ax_legend.legend(handles_right, allele_names, prop = { "size": 10 }, loc="lower left", ncol=2).set_title("Allele group")


    #heatmap
    ax_heatmap = plt.subplot(gs[1, 0])
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', [(1, 1, 1), rgb_color])
    im = ax_heatmap.imshow(cm, interpolation='nearest', cmap=custom_cmap, aspect='auto')
    for i in range(num_samples):
        for j in range(num_references):
            if count_matrix[i, j] > 1:
                ax_heatmap.text(j, i, int(count_matrix[i, j]), ha='center',color='Black', va='center', fontsize=10)
    ax_heatmap.set_xticks(np.arange(num_references))
    ax_heatmap.set_yticks(np.arange(num_samples))
    ax_heatmap.set_xticklabels(reference_labels, rotation=90, ha='center', fontsize=14)
    ax_heatmap.set_yticklabels(samples, fontsize=14)
    ax_heatmap.set_xticks(np.arange(-.5, num_references, 1), minor=True)
    ax_heatmap.set_yticks(np.arange(-.5, num_samples, 1), minor=True)

    ax_heatmap.grid(which='minor', color='w', linestyle='-', linewidth=4)
    population_change_indices = filtered_data['Population'].ne(filtered_data['Population'].shift()).to_numpy().nonzero()[0]
    for idx in population_change_indices[1:]:
        ax_heatmap.axhline(y=idx - 0.5, color='black', linewidth=0.8)


    ax_heatmap.grid(which='minor', color='w', linestyle='-', linewidth=4)
    ax_heatmap.tick_params(which='minor', bottom=False, left=False)
    ax_heatmap.set_ylabel('Samples', fontsize=14)
    ax_heatmap.set_xlabel('Allele name', fontsize=14)
    f = function.replace("/", "_")
    plt.savefig(f"{output}/{f}_{region}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output}/{f}_{region}.svg", dpi=300, bbox_inches='tight')

    plt.close(fig)



def main(path) -> None:
    data_known = open_files(data_path=f"{path}/annotation_report_known.xlsx")
    data_novel = open_files(data_path=f"{path}/annotation_report_novel.xlsx")
    data = pd.concat([data_known, data_novel])

    output = f"{path}/figure/heatmap"
    os.makedirs(output, exist_ok=True)

    pivot_df = make_pivot_table(data)

    functions = pivot_df['Function'].unique()
    regions = pivot_df['Region'].unique()

    for function, region in product(functions, regions):
        create_matrix(pivot_df, function=function, region=region, output=output)


if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr/annotation"
    main(path)

