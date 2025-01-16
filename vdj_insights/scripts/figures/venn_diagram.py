import pandas as pd
import json
import os
import matplotlib.pyplot as plt

from venny4py.venny4py import *


def open_files(data_path: str) -> pd.DataFrame:
    data = pd.read_excel(data_path)
    return data


def open_json(data_path: str) -> dict:
    with open(data_path, "r") as file:
        data = json.load(file)
    return data

def make_pivot_table(data: pd.DataFrame) -> pd.DataFrame:
    pivot_table = data.pivot_table(index=['Population', 'Region', 'Segment'],  columns='Short name', aggfunc='size',fill_value=0).reset_index()
    return pivot_table


def make_venn_diagram(data: pd.DataFrame, region: str, status: str, output: str, immuno_region: str) -> None:
    populations = data['Population'].unique()

    sets = {}
    for pop in populations:
        segments = set(data[data['Population'] == pop][data.columns[3:]].columns[
                           data[data['Population'] == pop].iloc[:, 3:].sum() > 0])
        sets[pop] = segments

    plt.figure(figsize=(20, 20))
    venny4py(sets=sets, dpi=600, ext='svg', legend_cols=4, out=f"{output}/{region}")
    plt.title(f"{status.capitalize()} segments of populations across {region}", size=9)
    plt.savefig(f"{output}/{region}/{status}_{region}.svg")
    plt.savefig(f"{output}/{region}/{status}_{region}.png")



def make_overall_venn(data: pd.DataFrame, status: str, output: str, immuno_region: str) -> None:
    populations = data['Population'].unique()

    sets = {}
    for pop in populations:
        segments = set()
        for region in data['Region'].unique():
            region_segments = set(data[data['Population'] == pop][data.columns[3:]].columns[
                                      data[data['Population'] == pop].iloc[:, 3:].sum() > 0])
            segments.update(region_segments)
        sets[pop] = segments

    plt.figure(figsize=(30, 30))

    venny4py(sets=sets, dpi=600, ext='svg', legend_cols=4, out=f"{output}/all_{status}")
    plt.title(f"{status.capitalize()} {immuno_region.upper()} segment in populations", size=9)
    plt.savefig(f"{output}/all_{status}/{status}_all_populations.svg")
    plt.savefig(f"{output}/all_{status}/{status}_all_populations.png")
    plt.close()


def main(path, immuno_region) -> None:
    output = f'{path}/figure/venn_diagram'
    os.makedirs( output, exist_ok=True)
    for status in ['known', 'novel']:
        data = open_files(data_path=f"{path}/annotation_report_{status}_rss.xlsx")
        if data['Population'].nunique() in [2, 3, 4]:
            pivot_df = make_pivot_table(data=data)

            regions = pivot_df['Region'].unique()
            for region in regions:
                region_df = pivot_df[pivot_df['Region'] == region]
                make_venn_diagram(region_df, region, status, output, immuno_region)
            make_overall_venn(pivot_df, status, output, immuno_region)


if __name__ == '__main__':
    path = "/home/jaimy/output/human/bcr"
    main(path, "bcr")


