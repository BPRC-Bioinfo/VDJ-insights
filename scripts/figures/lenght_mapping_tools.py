import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib_venn import venn3


def load_data(path: str) -> pd.DataFrame:
    data_known = pd.read_excel(f"{path}/annotation/annotation_report_known_rss.xlsx")
    data_novel = pd.read_excel(f"{path}/annotation/annotation_report_novel_rss.xlsx")
    combined_data = pd.concat([data_known, data_novel])
    combined_data['sequence_length'] = combined_data['End coord'] - combined_data['Start coord']
    return combined_data


def create_violin_plot(data: pd.DataFrame, output: str) -> None:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(14, 8))
    sns.violinplot(
        data=data,
        x='tool',
        y='sequence_length',
        hue='Segment',
        palette="Set2"
    )
    plt.title('Sequence length by mapping tool', fontsize=16)
    plt.xlabel('Mapping tool', fontsize=14)
    plt.ylabel('Sequence length', fontsize=14)
    plt.legend(title='Segment', fontsize=12, loc='upper right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig(f"{output}/violin_plot.svg", dpi=300)
    plt.savefig(f"{output}/violin_plot.pdf", dpi=300)
    plt.close()


def create_venn_diagram(data: pd.DataFrame, output: str) -> None:
    tool_groups = data.groupby('tool')['Short name'].apply(set)
    if len(tool_groups) != 3:
        raise ValueError("The Venn diagram requires exactly three tools for comparison.")

    (tool1, tool2, tool3) = tool_groups
    plt.figure(figsize=(10, 8))
    venn = venn3(
        [tool1, tool2, tool3],
        set_labels=tool_groups.index.tolist()
    )
    plt.title("Comparing found genenames across tools", fontsize=16, fontweight='bold')
    plt.tight_layout()

    plt.savefig(f"{output}/venn_diagram.svg", dpi=300)
    plt.savefig(f"{output}/venn_diagram.pdf", dpi=300)
    plt.close()


def main(path: str) -> None:
    output = f"{path}/figure/length_mapping_tool"
    os.makedirs(output, exist_ok=True)

    plot_data = load_data(path)
    create_violin_plot(data=plot_data, output=output)
    create_venn_diagram(data=plot_data, output=output)


if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr"
    main(path)

