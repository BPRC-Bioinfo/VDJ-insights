import json
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def open_json(data_path: str) -> dict:
    with open(data_path, "r") as file:
        data = json.load(file)
    return data


def make_dataframe(data: pd.DataFrame) -> pd.DataFrame:
    all_regions = set()
    for regions_data in data.values():
        all_regions.update(regions_data.keys())
    rows = []
    for assembly, regions_data in data.items():
        found_regions = set()
        for region, info in regions_data.items():
            print(info)
            rows.append({
                "Region": region,
                "Assembly": assembly.split("-")[0],
                "Assembly_Type": info["Extraction status"]
            })
            found_regions.add(region)

        for missing_region in all_regions - found_regions:
            rows.append({
                "Region": missing_region,
                "Assembly": assembly.split("-")[0],
                "Assembly_Type": "Not found"
            })

    df = pd.DataFrame(rows)
    return df

def make_plot(df: pd.DataFrame, output: str) -> None:
    pivot_df = df.groupby(["Region", "Assembly_Type"]).size().unstack(fill_value=0)
    saturated_colors = ["#FC8D59", "#FEE08B", "#4575B4", "#FC8D59", "#FEE08B", "#4575B4"]
    cmap = ListedColormap(saturated_colors[:len(pivot_df.columns)])

    fig, ax = plt.subplots(figsize=(10, 5))
    pivot_df.plot(kind="bar", stacked=True, ax=ax, colormap=cmap)

    for container in ax.containers:
        ax.bar_label(container, label_type="center", fmt="%d",
                     labels=[int(v) if v > 0 else '' for v in container.datavalues])

    plt.title("Found flanking genes in samples for extracting the region of intrest")
    plt.ylabel("Count")
    plt.xlabel("Region")
    plt.xticks(rotation=0)
    legend = plt.legend(title="Flanking genes", bbox_to_anchor=(0.5, -0.10), loc="upper center", ncol=3, frameon=False)
    plt.tight_layout()

    plt.savefig(f"{output}/broken_regions.svg", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output}/broken_regions.pdf", dpi=300, bbox_inches='tight')


def main(path: str) -> None:
    output = f"{path}/annotation/figure/broken_regions"
    os.makedirs(output, exist_ok=True)

    data = open_json(data_path=f"{path}/broken_regions.json")
    df = make_dataframe(data)
    make_plot(df, output)

if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr_v2"
    main(path)