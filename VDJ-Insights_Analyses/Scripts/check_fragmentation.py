from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from Bio import SeqIO
import re
from matplotlib.lines import Line2D


annotation_dir = Path("/home/susan/PyCharm/output/Human/scaffolds/bcr/annotation/annotation_report_all.xlsx")

annotation_data = pd.read_excel(annotation_dir)

for region_name, regions_df in annotation_data.groupby("Region"):
    fig, axes = plt.subplots(10, 10, figsize=(60, 50))
    axes = axes.flatten()
    for index, (path, region_df) in enumerate(regions_df.groupby("Path")):
        record = SeqIO.read(path, "fasta")
        sequence = str(record.seq).upper()

        contigs = re.finditer(f'N{{{100},}}', sequence)
        N_coords = [(match.start(), match.end() - 1) for match in contigs]

        segment_coords = list(region_df[["Short name", "Segment", "Status", "Start coord", "End coord"]].itertuples(index=False, name=None))
        ax = axes[index]

        for name, segment, status, start, end in segment_coords:
            segment_colors = {'V': '#8f9cd2', 'D': '#c989a1', 'J': '#92c989'}
            status_width = {'Known': 0.20, 'Novel': 0.30}
            ax.vlines([start, end], ymin=status_width[status], ymax=1 - status_width[status], color=segment_colors[segment], linewidth=1, transform=ax.get_xaxis_transform())

        for start, end in N_coords:
            ax.vlines([start, end], ymin=0, ymax=1, color='#000000', linewidth=1, transform=ax.get_xaxis_transform())

        contig_bounds = []
        last_end = 0
        for n_start, n_end in N_coords:
            contig_bounds.append((last_end, n_start - 1))
            last_end = n_end + 1
        contig_bounds.append((last_end, len(sequence) - 1))

        contig_segment_counts = []
        for c_start, c_end in contig_bounds:
            count = sum((s_start >= c_start and s_end <= c_end) for _, _, _, s_start, s_end in segment_coords)
            contig_segment_counts.append(count)

        segment_summary = " | ".join(str(n) for n in contig_segment_counts)

        if len(segment_coords) == 0 and len(N_coords) == 0:
            ax.set_visible(False)
            continue

        sample = region_df["Sample"].unique()[0]
        region = region_df["Region"].unique()[0]
        ax.set_title(f"{sample} | {len(segment_coords)} segmenten | {len(N_coords) + 1} contigs \n {segment_summary}", fontsize=10)

        ax.set_ylim(0, 1.2)
        ax.set_yticks([])
        ax.set_xlim(0, len(sequence))
        ax.set_xlabel("Position", fontsize=6)
        ax.tick_params(axis='x', labelsize=6)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x)}'))

    legend_elements = [
        Line2D([0], [0], color='#8f9cd2', lw=2, label='V-segment'),
        Line2D([0], [0], color='#c989a1', lw=2, label='D-segment'),
        Line2D([0], [0], color='#92c989', lw=2, label='J-segment'),
        Line2D([0], [0], color='#000000', lw=2, label='N-region')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=12, frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    
    fig.savefig(f"/mnt/nanopore/Susan/macacca/human/{region_name}_grid.svg", format="svg")
    plt.savefig(f"/mnt/nanopore/Susan/macacca/human/{region_name}_grid.pdf", format="pdf")
    plt.close(fig)