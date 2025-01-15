import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

plt.style.use('ggplot')


def open_files(data_path: str) -> pd.DataFrame:
    data = pd.read_excel(data_path)
    return data


def main(path) -> None:
    data = open_files(data_path=f"{path}/annotation_report_known_rss.xlsx")

    data['Subfamily'] = data['Short name'].apply(lambda x: re.split(r'[-*]', x)[0])

    regions = data['Region'].unique()
    for region in regions:
        region_data = data[data['Region'] == region]

        df_grouped = region_data.groupby(["Segment", "Subfamily"]).size().reset_index(name='Counts')
        df_grouped['Percentage'] = df_grouped.groupby('Segment')['Counts'].transform(lambda x: x / x.sum() * 100)
        pivot_df = df_grouped.pivot(index="Segment", columns="Subfamily", values="Percentage").fillna(0)

        segment_counts = df_grouped.groupby('Segment')['Counts'].sum()

        ax = pivot_df.plot(
            kind='barh',
            stacked=True,
            figsize=(12, 6),
            width=0.8,
        )

        plt.xlim(0, 100)

        y_labels = [f"{seg} (N={segment_counts[seg]})" for seg in pivot_df.index]
        ax.set_yticks(range(len(pivot_df.index)))
        ax.set_yticklabels(y_labels)

        samples = region_data['Sample'].nunique()
        plt.title(f"Distribution found known V(D)J subfamily of region: {region} (N={samples})", fontsize=16)
        plt.xlabel("Percentage", fontsize=14)
        plt.ylabel("Immune segment", fontsize=14)

        plt.legend(
            title="Subfamily",
            bbox_to_anchor=(0.5, -0.2),
            loc="upper center",
            ncol=9,
            fontsize=10,
            title_fontsize=12,
            frameon=False
        )

        plt.tight_layout()

        os.makedirs(f"{path}/figure/subfamily_barplot", exist_ok=True)
        plt.savefig(f"{path}/figure/subfamily_barplot/{region}_subfamilies.png", bbox_inches="tight", dpi=300)
        plt.close()


if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr"
    main(path)
