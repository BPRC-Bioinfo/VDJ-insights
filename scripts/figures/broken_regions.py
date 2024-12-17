import json
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches


def open_json(data_path: str) -> dict:
    with open(data_path, "r") as file:
        data = json.load(file)
    return data


def make_pivot_table(data: pd.DataFrame) -> pd.DataFrame:
    data = data.groupby(['Flanking_genes', 'Sample'], as_index=False).agg({'N_contigs': 'sum'})
    pivot_table = data.pivot(index='Flanking_genes', columns='Sample', values='N_contigs')
    return pivot_table


def make_plot(pivot_table: pd.DataFrame, df: pd.DataFrame, output: str) -> None:
    data_array = pivot_table.values

    samples = df['Sample'].unique()
    flanking_genes = df['Flanking_genes'].unique()
    num_samples = len(samples)
    num_references = len(flanking_genes)
    fig, ax = plt.subplots(figsize=(num_references * 0.5 + 25, num_samples * 0.3 + 10))

    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', [(1, 1, 1), (1, 1, 1)])
    cax = ax.imshow(data_array, cmap=custom_cmap)
    ax.set_xticks(np.arange(len(pivot_table.columns)))
    ax.set_yticks(np.arange(len(pivot_table.index)))
    ax.set_xticklabels(pivot_table.columns, ha='center', fontsize=10, rotation=90)
    ax.set_yticklabels(pivot_table.index, va='center', fontsize=10)

    for i in range(len(pivot_table.index)):
        for j in range(len(pivot_table.columns)):
            point = data_array[i, j]
            if np.isnan(point):
                ax.text(j, i, "×", ha="center", va="center", color="red", label="Not found", fontweight='bold')
            elif point == 2:
                ax.text(j, i, "-", ha="center", va="center", color="orange", label="Not extract", fontweight='bold')
            else:
                ax.text(j, i, "✓", ha="center", va="center", color="green", label="Extract", fontweight='bold')

    not_found_patch = mpatches.Patch(color='red', label='Not found (×)')
    not_extracted_patch = mpatches.Patch(color='orange', label='Not extracted (-)')
    extracted_patch = mpatches.Patch(color='green', label='Extracted (✓)')
    ax.legend(handles=[extracted_patch, not_extracted_patch, not_found_patch], loc='center left',
              bbox_to_anchor=(1.01, 0.5))

    ax.set_title('Schematic representation of extracted regions with flanking genes per sample')
    plt.tight_layout()

    plt.savefig(f"{output}/broken_regions.svg", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output}/broken_regions.pdf", dpi=300, bbox_inches='tight')


def main(path: str) -> None:
    output = f"{path}/figure/broken_regions"
    os.makedirs(output, exist_ok=True)

    data = open_json(data_path=f"{path}/broken_regions.json")
    rows = []
    for sample, genes in data.items():
        for flanking_gene, contigs in genes.items():
            if contigs['5_Contig'] ==  contigs['3_Contig']:
                num_contigs = 1
            else:
                num_contigs = 2
            receptor_data ={
                "TMEM121--": "IGH",
                "RPIA-LSP1P4": "IGK",
                "GNAZ-TOP3B": "IGL",
                "SALL2-DAD1": "TRA",
                "MGAM2-EPHB6": "TRB",
                "EPDR1-VPS41": "TRG",
                "SALL2-DAD1": "TRD"
            }
            rows.append({
                'Sample': sample.split("-")[1].strip(".fa"),
                'Flanking_genes': receptor_data[flanking_gene],
                'N_contigs': num_contigs
            })
    data = pd.DataFrame(rows)
    pivot_table = make_pivot_table(data)
    #print(pivot_table)
    make_plot(pivot_table, data, output)

if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr_v2"
    main(path)