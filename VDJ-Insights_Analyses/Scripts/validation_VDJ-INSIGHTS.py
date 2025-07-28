import warnings
from unittest.mock import right

warnings.filterwarnings("ignore")

import requests
from bs4 import BeautifulSoup
import os
import pandas as pd
import re
from Bio import SeqIO
import subprocess
from multiprocessing import Pool
import matplotlib.pyplot as plt
from glob import glob
from matplotlib import pyplot as plt
from matplotlib_venn import venn2_unweighted

pd.set_option('display.max_columns', None)

def download_fasta(s_accession_number, receptor, output, region):
    if not os.path.isdir(output):
        os.makedirs(output)

    for extension, fmt in zip(["txt", "fasta"], ["IMGT", "FASTA"]):
        url = f"https://www.imgt.org/ligmdb/view.action?format={fmt}&id={s_accession_number}"
        output_file = f"{output}/{fmt}"

        if not os.path.isdir(output_file):
            os.makedirs(output_file)
        if not os.path.exists(f"{output_file}/{s_accession_number}_{region}.{extension}"):
            response = requests.get(url)
            if response.status_code != 200:
                return False
            soup = BeautifulSoup(response.text, "html.parser")
            pre_tag = soup.find("pre")
            if not pre_tag:
                return False

            data = pre_tag.text.strip()
            with open(f"{output_file}/{s_accession_number}_{region}.{extension}", "w") as file:
                file.write(data + "\n")

    segments = []
    with open(f"{output}/IMGT/{s_accession_number}_{region}.txt", "r") as file:
        for record in SeqIO.parse(file, "embl"):
            for feature in record.features:
                type = feature.type
                if ("REGION" in type or "segment" in type) and (feature.type.startswith("V") or feature.type.startswith("D") or feature.type.startswith("J")):
                    location = feature.location
                    start = location.start
                    end = location.end
                    strand = "+" if location.strand == 1 else "-"
                    allele = feature.qualifiers.get("IMGT_allele", ["-"])[0]
                    if allele == "-":
                        allele = feature.qualifiers.get("allele", ["-"])[0]

                    segment = str(re.search(r'([VDJ])', feature.type).group(1))
                    functional = feature.qualifiers.get("functional", ["-"])[0]
                    pseudo = feature.qualifiers.get("pseudo", ["-"])[0]
                    ORF = feature.qualifiers.get("ORF", ["-"])[0]

                    if functional != "-":
                        function = "functional"
                    elif pseudo != "-":
                        function = "pseudo"
                    elif ORF != "-":
                        function = "ORF"

                    messenger = feature.qualifiers.get(function, ["-"])[0]
                    if allele != "-":
                        segments.append((start, end, allele, segment, strand, function, messenger))

    imgt_data = pd.DataFrame(segments, columns=["Start coord", "End coord", "Short name", "Segment", "Strand", "Function", "Messenger"])
    imgt_data.to_excel(f"{output}/IMGT/{s_accession_number}_{region}.xlsx", index=False)


def process_task(task):
    species, s_accession_number, region = task
    #print(species, s_accession_number, region)

    def find_segment_key(segment_name):
        segments = {
            "TR": ["TRB", "TRA", "TRD", "TRG"],
            "IG": ["IGH", "IGK", "IGL"]
        }
        for key, values in segments.items():
            if segment_name in values:
                return key
        return None

    receptor = find_segment_key(region)
    s_species = species.replace(" ", "_")
    output = os.path.join("validation_v2", s_species, f"{s_accession_number}_{region}")

    download_fasta(s_accession_number, receptor, output, region)


    a = f"{output}/FASTA/"
    o = f"{output}/VDJ-INSIGHTS"
    try:
        cmd = f'python -m tool.vdj_insights.scripts.vdj_insights annotation -i {a} -s "{species}" -r {receptor} --default -o {o} -t 4'
        if not os.path.exists(o + "/annotation/annotation_report_known.xlsx"):
            subprocess.run(cmd , shell=True)
        #subprocess.run(cmd , shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        imgt_data = pd.read_excel(f"{output}/IMGT/{s_accession_number}_{region}.xlsx")
        vdj_insights_data = pd.read_excel(f"{o}/annotation/annotation_report_all.xlsx", usecols=[2, 3, 4, 5, 8, 9, 30, 31])
        compare = imgt_data.merge(vdj_insights_data, on=["Start coord", "End coord", "Strand"], how="outer", suffixes=["_IMGT", "_VDJ"], indicator=True)
        compare["Short name match"] = compare["Short name_IMGT"].fillna("") == compare["Short name_VDJ"].fillna("")
        compare["Function match"] = compare["Function_IMGT"].fillna("") == compare["Function_VDJ"].fillna("")
        for index, row in compare.iterrows():
            if str(row['Short name_IMGT']) in str(row["Similar references"]):
                compare.loc[index, "Short name match"] = True

        def merge_on_coordinate(compare: pd.DataFrame, coordinate: str) -> pd.DataFrame:
            merged_rows = []
            for val, group in compare.groupby(coordinate, sort=False):
                if len(group) > 1:
                    merged = group.iloc[0].copy()
                    if coordinate == 'Start coord':
                        merged['End coord'] = group['End coord'].max()
                    else:
                        merged['Start coord'] = group['Start coord'].min()
                    merged['_merge'] = 'both'
                    merged['Short name match'] = True
                    merged['Function match'] = group['Function match'].all()
                    merged_rows.append(merged)
                else:
                    merged_rows.append(group.iloc[0].copy())
            return pd.DataFrame(merged_rows).reset_index(drop=True)

        compare = merge_on_coordinate(compare, 'Start coord')
        compare = merge_on_coordinate(compare, 'End coord')

        imgt_data_pivot = imgt_data.pivot_table(index=['Segment'], columns=['Segment'], aggfunc='size', fill_value=0)
        vdj_insights_pivot = vdj_insights_data.pivot_table(index=['Segment'], columns=['Segment'], aggfunc='size', fill_value=0)

        function_order = ["functional", "ORF", "pseudo", "Total"]
        segment_order = ["V", "D", "J", "Total"]

        imgt_function_pivot = imgt_data.pivot_table(index=['Function'], columns=['Segment'], aggfunc='size', fill_value=0)
        imgt_function_pivot = imgt_function_pivot.reindex(columns=segment_order[:-1], fill_value=0)
        imgt_function_pivot["Total"] = imgt_function_pivot.sum(axis=1)
        imgt_function_pivot.loc["Total"] = imgt_function_pivot.sum(axis=0)
        imgt_function_pivot = imgt_function_pivot.reset_index()
        imgt_function_pivot["Function"] = pd.Categorical(imgt_function_pivot["Function"], categories=function_order, ordered=True)
        imgt_function_pivot = imgt_function_pivot.sort_values("Function")
        imgt_function_pivot = imgt_function_pivot[["Function"] + segment_order]

        vdj_insights_function_pivot = vdj_insights_data.pivot_table(index=['Function'], columns=['Segment'], aggfunc='size',fill_value=0)
        vdj_insights_function_pivot = vdj_insights_function_pivot.reindex(columns=segment_order[:-1], fill_value=0)
        vdj_insights_function_pivot["Total"] = vdj_insights_function_pivot.sum(axis=1)
        vdj_insights_function_pivot.loc["Total"] = vdj_insights_function_pivot.sum(axis=0)
        vdj_insights_function_pivot = vdj_insights_function_pivot.reset_index()
        vdj_insights_function_pivot["Function"] = pd.Categorical(vdj_insights_function_pivot["Function"], categories=function_order, ordered=True)
        vdj_insights_function_pivot = vdj_insights_function_pivot.sort_values("Function")
        vdj_insights_function_pivot = vdj_insights_function_pivot[["Function"] + segment_order]


        combined_false = compare[(compare['Function match'] == False) & (compare['_merge'] == "both")]
        combined_false_pivot = combined_false.pivot_table(index=['Function_IMGT'], columns=['Segment_IMGT'], aggfunc='size', fill_value=0)
        combined_false_pivot = combined_false_pivot.reindex(columns=segment_order[:-1], fill_value=0)
        combined_false_pivot["Total"] = combined_false_pivot.sum(axis=1)
        combined_false_pivot.loc["Total"] = combined_false_pivot.sum(axis=0)
        combined_false_pivot = combined_false_pivot.reset_index()
        combined_false_pivot["Function_IMGT"] = pd.Categorical(combined_false_pivot["Function_IMGT"],categories=function_order, ordered=True)
        combined_false_pivot = combined_false_pivot.sort_values("Function_IMGT")
        combined_false_pivot = combined_false_pivot[["Function_IMGT"] + segment_order]


        combined_true = compare[(compare['Function match'] == True) & (compare['_merge'] == "both")]
        combined_true_pivot = combined_true.pivot_table(index=['Function_IMGT'], columns=['Segment_IMGT'], aggfunc='size',fill_value=0)
        combined_true_pivot = combined_true_pivot.reindex(columns=segment_order[:-1], fill_value=0)
        combined_true_pivot["Total"] = combined_true_pivot.sum(axis=1)
        combined_true_pivot.loc["Total"] = combined_true_pivot.sum(axis=0)
        combined_true_pivot = combined_true_pivot.reset_index()
        combined_true_pivot["Function_IMGT"] = pd.Categorical(combined_true_pivot["Function_IMGT"],categories=function_order, ordered=True)
        combined_true_pivot = combined_true_pivot.sort_values("Function_IMGT")
        combined_true_pivot = combined_true_pivot[["Function_IMGT"] + segment_order]

        plt.figure(figsize=(10, 6))
        bars_total = plt.bar(combined_true_pivot["Function_IMGT"], combined_true_pivot["Total"],color="lightgray", label="Totaal IMGT", alpha=0.7)
        bars_false = plt.bar(combined_false_pivot["Function_IMGT"], combined_false_pivot["Total"],color="red", label="Mismatch VDJ-Insights", alpha=0.8)
        for bar in bars_total:
            plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{int(bar.get_height())}",
                     ha='center', va='bottom', fontsize=10, color='black')

        for bar in bars_false:
            plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{int(bar.get_height())}",
                     ha='center', va='bottom', fontsize=10, color='black')

        plt.xlabel("Function IMGT")
        plt.ylabel("Aantal")
        plt.title("Vergelijking VDJ-Insights vs IMGT (Totaal vs Mismatches)")
        plt.legend()
        plt.xticks(rotation=45)
        plt.savefig(output, dpi=300, bbox_inches="tight")


        count = compare[(compare['_merge'] == "both") & (compare['Segment_IMGT'] == "V")]['Function match'].value_counts().to_dict()

        counts = compare['_merge'].value_counts()
        both = counts.get('both', 0)
        left_only = counts.get('left_only', 0)
        right_only = counts.get('right_only', 0)

        total = both + left_only
        pct = (both / total * 100) if total > 0 else 0.0

        print(f'{species:<25} {s_accession_number:<15} {region:<5} both: {both:<5} {pct:<10.2f}% {right_only:<2} extra indentified by VDJ-Insights')
        #print(f"{species} {s_accession_number} {region} only V: {count} {(count[True] / (count[True] + count[False])) * 100:.2f}% {found}")


        with pd.ExcelWriter(f"{output}/compair.xlsx") as writer:
            compare.to_excel(writer, sheet_name="Segment", index=False)
            imgt_data.to_excel(writer, sheet_name="Data IMGT", index=False)
            vdj_insights_data.to_excel(writer, sheet_name="Data VDJ-Insights", index=False)
            imgt_data_pivot.to_excel(writer, sheet_name="Count IMGT", index=False)
            vdj_insights_pivot.to_excel(writer, sheet_name="Count VDJ-Insights", index=False)
            imgt_function_pivot.to_excel(writer, sheet_name="IMGT functions", index=True)
            vdj_insights_function_pivot.to_excel(writer, sheet_name="VDJ-Insights functions", index=True)
            combined_false.to_excel(writer, sheet_name="False functions", index=False)
            combined_false_pivot.to_excel(writer, sheet_name="False predict", index=True)
            combined_true_pivot.to_excel(writer, sheet_name="True predict", index=True)

    except Exception as e:
        print(f"error: {e}")




def main():
    tasks = [
        ('Mus musculus', 'BK063712', 'IGH'),
        ('Mus musculus', 'BK063713', 'IGH'),
        ('Mus musculus', 'BK063714', 'IGH'),
        ('Mus musculus', 'IMGT000188', 'IGH'),

        ('Mus musculus', 'IMGT000117', 'IGL'),
        ('Mus musculus', 'IMGT000126', 'IGL'),
        ('Mus musculus', 'IMGT000134', 'IGL'),
        ('Mus musculus', 'IMGT000211', 'IGL')
    ]
    tasks = [('Gorilla gorilla gorilla', 'BK063604', 'IGH'),
             ('Gorilla gorilla gorilla', 'BK063605', 'IGH'),
             ('Gorilla gorilla gorilla', 'BK063606', 'IGH'),
             ('Gorilla gorilla gorilla', 'BK063607', 'IGH'),
             ('Gorilla gorilla gorilla', 'BK063600', 'IGK'),
             ('Gorilla gorilla gorilla', 'BK063603', 'IGK'),
             ('Gorilla gorilla gorilla', 'BK063601', 'IGK'),
             ('Gorilla gorilla gorilla', 'BK063602', 'IGK'),
             ('Gorilla gorilla gorilla', 'BK063608', 'IGL'),
             ('Gorilla gorilla gorilla', 'BK063609', 'IGL'),
             ('Gorilla gorilla gorilla', 'BK063610', 'IGL'),
             ('Gorilla gorilla gorilla', 'BK063611', 'IGL'),

             ('Homo sapiens', 'BK063799', 'IGH'),
             ('Homo sapiens', 'IMGT000024', 'TRA'),
             ('Homo sapiens', 'IMGT000021', 'TRB'),
             ('Homo sapiens', 'IMGT000011', 'TRG'),

             ('Macaca mulatta', 'BK063715', 'IGH'),
             ('Macaca mulatta', 'BK063717', 'IGL'),
             ('Macaca mulatta', 'BK063716', 'IGK'),

             #('Macaca mulatta', 'IMGT000076', 'TRA'),
             #('Macaca mulatta', 'IMGT000073', 'TRB'),
             #('Macaca mulatta', 'IMGT000059', 'TRG'),

             ('Mus musculus', 'BK063712', 'IGH'),
             ('Mus musculus', 'BK063713', 'IGH'),
             ('Mus musculus', 'BK063714', 'IGH'),
             ('Mus musculus', 'IMGT000188', 'IGH'),

             ('Mus musculus', 'IMGT000117', 'IGL'),
             ('Mus musculus', 'IMGT000126', 'IGL'),
             ('Mus musculus', 'IMGT000134', 'IGL'),
             ('Mus musculus', 'IMGT000211', 'IGL')
             ]

    with Pool(processes=6) as pool:
        pool.map(process_task, tasks)


    show_data = True
    if show_data:
        species_totals = {}

        right_only = pd.DataFrame()
        left_only = pd.DataFrame()

        for species in glob("validation_v2/*"):
            species_name = os.path.basename(species)
            for region in glob(f"{species}/*/compair.xlsx"):
                id = region.split("/")[2]
                if species_name == "Macaca_mulatta" and "TR" in region: #niet mee voor validatie
                    continue
                segments = pd.read_excel(region, sheet_name="Segment")
                both_df = segments[(segments['_merge'] == "both")]
                for group_name, s_df in both_df.groupby("Function match"):
                    values = s_df["Function_IMGT"].value_counts().to_dict()
                    if species_name not in species_totals:
                        species_totals[species_name] = {"Functional_True": 0,
                                                        "Functional_False": 0,
                                                        "ORF_True": 0,
                                                        "ORF_False": 0,
                                                        "Pseudo_True": 0,
                                                        "Pseudo_False": 0,
                                                        "Totaal": 0,
                                                        "True": 0,
                                                        "False": 0,
                                                        "only_imgt": 0,
                                                        "only_vdj": 0,
                                                        "Totaal_rows": 0,
                                                        "Totaal_both": 0,
                                                        "Totaal_vdj": 0
                                                        }


                    species_totals[species_name][f"Functional_{group_name}"] += values.get("functional", 0)
                    species_totals[species_name][f"ORF_{group_name}"] += values.get("ORF", 0)
                    species_totals[species_name][f"Pseudo_{group_name}"] += values.get("pseudo", 0)

                percentage_gevonden_segmenten = segments['_merge'].value_counts().to_dict()
                count = segments[segments['_merge'] == "both"]['Function match'].value_counts().to_dict()
                species_totals[species_name]["True"] += count.get(True, 0)
                species_totals[species_name]["False"] += count.get(False, 0)
                species_totals[species_name]["Totaal"] += (count.get(True, 0) + count.get(False, 0))

                species_totals[species_name]["only_imgt"] += segments[segments['_merge'] == "left_only"].shape[0]
                species_totals[species_name]["only_vdj"] += segments[segments['_merge'] == "right_only"].shape[0]

                species_totals[species_name]["Totaal_rows"] += segments[segments['_merge'] != "right_only"].shape[0]
                species_totals[species_name]["Totaal_both"] += segments[segments['_merge'] == "both"].shape[0]

                total_vdj_df = segments[(segments['_merge'] == "both") | (segments['_merge'] == "right_only")]
                species_totals[species_name]["Totaal_vdj"] += total_vdj_df.shape[0]

                left_only_df = segments[segments['_merge'] == "left_only"]
                right_only_df = segments[segments['_merge'] == "right_only"]
                left_only_df['ID'] = id
                right_only_df['ID'] = id
                left_only = pd.concat([left_only, left_only_df])
                right_only = pd.concat([right_only, right_only_df])

        for key, value in species_totals.items():
            venn2_unweighted(subsets=(value['only_imgt'], value['only_vdj'], value['Totaal_both']), set_labels=('IMGT', 'VDJ-Insights'))
            title = key.replace("_", " ")
            plt.title(f"{title}")
            plt.savefig(f"{key}.svg", dpi=300, bbox_inches='tight')
            plt.close()

            print("-" * 100)
            print(key.replace("_", " "))
            print(f"total rows: {value['Totaal_rows']}  total VDJ: {value['Totaal_vdj']} total found both: {value['Totaal_both']} acc: {(value['Totaal_both'] / value['Totaal_rows']) * 100:.2f}%")
            print(f"VDJ-Insights correctly annotated {value['Functional_True']} functional, {value['ORF_True']} ORF, and {value['Pseudo_True']} "
                  f"pseudo segments, with an overall accuracy of {(value['True'] / (value['True'] + value['False'])) * 100:.2f}% ({value['True']}/{value['True'] + value['False']}). Out of"
                  f" {value['Totaal']} total segments in {key}, {value['False']} segments were assigned with an alternative functionality.")



        #overall
        tr = 0
        fl = 0
        both = 0
        total = 0
        total_vdj = 0
        for key, value in species_totals.items():
            tr += value['True']
            fl += value['False']
            both += value['Totaal_both']
            total += value['Totaal_rows']
            total_vdj += value['Totaal_vdj']


        print("-" * 100)
        print(
            f"With a total of {total} segment annotated by IMGT, {both} were identified by VDJ-Insights ({(both / total) * 100:.2f}%). "
            f"{left_only.shape[0]} segments annotated by IMGT were not detected by VDJ-Insights."
            f"VDJ-Insights detected {total_vdj} segments, which is {right_only.shape[0]} more than IMGT."
            f"These primarily consist of small D segments that are randomly mapped to the genome, potentially leading to false positives. "
            f"VDJ-Insights correctly annotated the functionality of {tr} out of {tr + fl} segments, achieving an overall functional annotation accuracy of ({(tr / (tr + fl)) * 100:.2f}%)."
        )
        #print(left_only[["Start coord", "End coord", "Short name_IMGT", "ID"]])
        print(f"The overall functionality: {(tr / (tr + fl)) * 100:.2f}% (True: {tr} / Total: {tr + fl})")
        print(f"The overall found segments: {(both / total) * 100:.2f}% (Both found: {both} / Total: {total})")
        #print(f"The overall found segments: {(both / total_vdj) * 100:.2f}% (Both found: {both} / Total VDJ: {total_vdj})")





if __name__ == "__main__":
    main()