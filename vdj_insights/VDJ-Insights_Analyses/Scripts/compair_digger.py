# Annotation comparison script for IMGT, Digger, and VDJ-Insights
# ---------------------------------------------------------------
# This script compares annotated gene segments from three tools:
# IMGT, Digger, and VDJ-Insights. It performs the following:
#
# - Reads configuration and annotations
# - Merges and compares coordinates and annotations
# - Identifies overlaps, differences, and unique annotations
# - Generates comparison tables and Venn diagrams
# - Aggregates results into summary Excel outputs
#
# Author: [Your Name]
# Date: [Date of publication]
# ----------------------------------------------------------------

import pandas as pd
import numpy as np
import glob
import warnings
from itertools import product
import os
import yaml
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted

warnings.simplefilter(action='ignore', category=FutureWarning)

# --------------------------
# Load configuration

# Initialize sets for IG and TR tracking across samples
imgt_IG_set = set()
imgt_TR_set = set()
imgt_total_set = set()
vdj_IG_set = set()
vdj_TR_set = set()
vdj_total_set = set()
digger_IG_set = set()
digger_TR_set = set()
digger_total_set = set()
# --------------------------
# Initialize result storage lists
compair_list = []
only_digger_list = []
only_vdj_insights_list = []
imgt_vdj_list = []
imgt_digger_list = []
vdj_digger_list = []
imgt_vdj_f_list = []
imgt_digger_f_list = []
vdj_digger_f_list = []

with open("config_human.yaml", "r") as file:
    config = yaml.safe_load(file)

VDJ_IMGT_path = config["Result_path"]
Digger_path = config["Digger_path"]
samples = config["samples"]
Region = config["region"]
species = config["species"].capitalize()

os.makedirs(os.path.join(VDJ_IMGT_path, "compair"), exist_ok=True)

# Determine regions based on IG or TR input
if Region == "IG":
    regions = ["IGH", "IGK", "IGL"]
elif Region == "TR":
    regions = ["TRA", "TRB", "TRA&TRD", "TRG"]
else:
    regions = ["IGH", "IGK", "IGL", "TRA&TRD", "TRB", "TRG"]

# Helper functions for plotting with structured inputs
def plot_venn(imgt_s, digger_s, vdj_s, sample, region, species, filename_prefix):
    region_title = region.upper()
    species_title = species.replace("_", " ").capitalize()

    if sample:
        title = f"Tool comparison of {region_title} in {sample}"
    elif region.lower() == "all":
        title = species_title
    else:
        title = f"{region_title}-region in the {species_title}"

    plt.figure(figsize=(6, 6))
    venn3_unweighted([imgt_s, digger_s, vdj_s],
                     ('IMGT', 'Digger', 'VDJ-Insights'),
                     set_colors=('#1B9E77', '#D95F02', '#7570B3'))
    if region == "all" and sample == None:
        title = species.replace("_", " ")
        title = title.capitalize()
    plt.title(title, fontsize=20)
    plt.savefig(f"{filename_prefix}.svg", format="svg", bbox_inches="tight")
    plt.savefig(f"{filename_prefix}.pdf", format="pdf", bbox_inches="tight")
    plt.close()

def plot_function_ven(df, sample=None, region=None, species=None, filename=None):
    title = None
    if region and species:
        region_title = region.upper()
        species_title = species.replace("_", " ").capitalize()
        
        if sample:
            title = f"Functionality comparison of {region_title} in {sample}"
        else:
            if species != "Human": 
                if region.lower() == "all":
                    title = species.replace("_"," ").capitalize()
                else:
                    title = f"Functionality comparison of {region_title} in {species_title}"
            elif species == "Human" and region.lower() == "all":
                if region.lower() == "all":
                    title = None
                else:
                    title = f"Functionality comparison of {region_title} in {species_title}"
            else:
                title = f"Functionality comparison of {region_title} in {species_title}"
    else:
        title = "Functionality comparison"

    plt.figure(figsize=(20, 6), facecolor='white')
    if "Sample" not in df.columns:
        df["Sample"] = "dummy"
    labels = ['A.', 'B.', 'C.']
    Functions = ["Functional", "ORF", "Pseudo"]
    for idx, Function in enumerate(Functions):
        Function_imgt_set = set(zip(
            df.loc[df["Function_imgt"] == Function, "Start coord"],
            df.loc[df["Function_imgt"] == Function, "End coord"],
            df.loc[df["Function_imgt"] == Function, "Function_imgt"],
            df.loc[df["Function_imgt"] == Function, "Sample"]))
        Function_vdj_set = set(zip(
            df.loc[df["Function_vdj-insights"] == Function, "Start coord"],
            df.loc[df["Function_vdj-insights"] == Function, "End coord"],
            df.loc[df["Function_vdj-insights"] == Function, "Function_vdj-insights"],
            df.loc[df["Function_vdj-insights"] == Function, "Sample"]))
        Function_digger_set = set(zip(
            df.loc[df["Function_digger"] == Function, "Start coord"],
            df.loc[df["Function_digger"] == Function, "End coord"],
            df.loc[df["Function_digger"] == Function, "Function_digger"],
            df.loc[df["Function_digger"] == Function, "Sample"]))

        plt.subplot(1, 3, idx + 1, facecolor='white')

        if species == "Human" and region == "all":
            plt.text(-0.15, 1.05, labels[idx], transform=plt.gca().transAxes, fontsize=28, fontweight='bold', va='top', ha='left')
            
        venn3_unweighted([Function_imgt_set, Function_digger_set, Function_vdj_set], ('IMGT', 'Digger', 'VDJ-Insights'), set_colors=('#1B9E77', '#D95F02', '#7570B3'))
        if Function == "Pseudo":
            #plt.title("Pseudogene segment comparison")
            plt.title("Pseudogene segment", fontsize=20)
        else:
            #plt.title(f"{Function} gene segment comparison")
            plt.title(f"{Function} gene segment", fontsize=20)

    if title is not None:
        plt.suptitle(title, fontsize=25)

    plt.tight_layout()
    plt.savefig(f"{filename}.svg", format='svg')
    plt.savefig(f"{filename}.pdf", format='pdf')
    plt.close()

for sample, region in product(samples, regions):

    #If samples IMGT and VDJ-Insight match with Digger
    #inladen
    if region == "TRA&TRD":
        imgt_path = glob.glob(os.path.join(VDJ_IMGT_path, f"{sample}_TRA", "IMGT", f"*_TRA.xlsx"))
        path_vdj = glob.glob(os.path.join(VDJ_IMGT_path, f"{sample}_TRA", "VDJ-INSIGHTS", "annotation", "annotation_report_all.xlsx"))
    else:
        imgt_path = glob.glob(os.path.join(VDJ_IMGT_path, f"{sample}_{region}", "IMGT", f"*_{region}.xlsx"))
        path_vdj = glob.glob(os.path.join(VDJ_IMGT_path, f"{sample}_{region}", "VDJ-INSIGHTS", "annotation", "annotation_report_all.xlsx"))

    '''
    #If samples IMGT and VDJ-Insights do not match with Digger
    if region == "TRA&TRD":
        imgt_path = glob.glob(f"{VDJ_IMGT_path}*_TRA/IMGT/*_TRA.xlsx")
        path_vdj = glob.glob(f"{VDJ_IMGT_path}*_TRA/VDJ-INSIGHTS/annotation/annotation_report_all.xlsx")
    else:
        imgt_path = glob.glob(f"{VDJ_IMGT_path}*_{region}/IMGT/*_{region}.xlsx")
        path_vdj = glob.glob(f"{VDJ_IMGT_path}*_{region}/VDJ-INSIGHTS/annotation/annotation_report_all.xlsx")
    '''
    imgt_file = pd.read_excel(imgt_path[0]) if imgt_path else None


    if path_vdj:
        vdj_file = pd.read_excel(path_vdj[0])
        if region == "TRA&TRD":
            digger_file_TRA = pd.read_csv(glob.glob(f"{Digger_path}{sample}_TRA_digger_combined.csv")[0])
            digger_file_TRD = pd.read_csv(glob.glob(f"{Digger_path}{sample}_TRD_digger_combined.csv")[0])
            digger_file = pd.concat([digger_file_TRA, digger_file_TRD], ignore_index=True)
        else:
            digger_file = pd.read_csv(glob.glob(f"{Digger_path}{sample}_{region}_digger_combined.csv")[0])

        digger_file = digger_file.drop_duplicates(subset=["start", "end"], keep="first")
        digger_file.rename(columns={"blast_match": "Segment_digger", "start": "Start coord", "end": "End coord","functional":"Function_digger"}, inplace=True)
        digger_file.loc[digger_file["Function_digger"] == "Functional", "Function_digger"] = "functional"
        digger_file["Function_digger"] = digger_file["Function_digger"].str.capitalize()
        digger_file["Function_digger"] = digger_file["Function_digger"].replace({"Orf": "ORF"})
        digger_file["Start coord"] = digger_file["Start coord"] - 1
        if region[0:2] == "IG":
            digger_file["Segment_digger"] = digger_file['Segment_digger'].str.extract(r'(IG[^_]+)')
        else:
            digger_file["Segment_digger"] = digger_file['Segment_digger'].str.extract(r'(TR[^_]+)')

        imgt_file.rename(columns={"Short name": "Segment_imgt", "Function":"Function_imgt"}, inplace=True)
        vdj_file.rename(columns={"Short name": "Segment_vdj-insights", "Function":"Function_vdj-insights"}, inplace=True)

        imgt = imgt_file[["Segment_imgt", "Start coord", "End coord","Function_imgt"]]
        vdj = vdj_file[["Segment_vdj-insights", "Start coord", "End coord", "Function_vdj-insights"]]
        digger = digger_file[["Segment_digger", "Start coord", "End coord", "Function_digger"]]

        def check(row, tool1, tool2):
            if pd.notna(row[f"Segment_{tool1}"]) and pd.isna(row[f"Segment_{tool2}"]):
                return f"Found in {tool1} only"
            elif pd.isna(row[f"Segment_{tool1}"]) and pd.notna(row[f"Segment_{tool2}"]):
                return f"Found in {tool2} only"
            else:
                if row[f"Segment_{tool1}"] == row[f"Segment_{tool2}"]:
                    return "Found in both tools"
                else:
                    return "Found in both tools, different segment name"

        def check_function(row ,tool1, tool2):
            if pd.notna(row[f"Function_{tool1}"]) and pd.isna(row[f"Function_{tool2}"]):
                return f"Only found in {tool1}"
            elif pd.isna(row[f"Function_{tool1}"]) and pd.notna(row[f"Function_{tool2}"]):
                return f"Only found in {tool2}"
            else:
                if pd.notna(row[f"Function_{tool1}"]) and row[f"Function_{tool1}"] == row[f"Function_{tool2}"]:
                    return "Functionality the same"
                elif pd.notna(row[f"Function_{tool1}"]) and row[f"Function_{tool1}"] != row[f"Function_{tool2}"]:
                    return "Functionality different"
                else:
                    return
        def extract_vdj(segment):
            if pd.notna(segment) and segment != "Missing segment":
                return segment[3]
            return "No Segment"

        def pivot_adds(ful_file, pivot_file, tool1, tool2):
            ful_file["VDJ"] = ful_file.apply(lambda row: extract_vdj(row[f"Segment_{tool1}"])
            if pd.notna(row[f"Segment_{tool1}"]) else extract_vdj(row[f"Segment_{tool2}"]), axis=1)

            pivot_counts = ful_file.groupby(['Short name match', 'VDJ']).size().unstack(fill_value=0)
            pivot_counts.rename(columns={"V": "V_count", "D": "D_count", "J": "J_count"}, inplace=True)
            pivot_file = pivot_file.merge(pivot_counts, on="Short name match", how="left").fillna(0)
            return pivot_file

        def pivot_function_adds(ful_file, pivot_file, tool1, tool2):
            ful_file["Function"] = ful_file.apply(lambda row: (row[f"Function_{tool1}"])
            if pd.notna(row[f"Function_{tool1}"]) else row[f"Function_{tool2}"], axis=1)

            pivot_counts = ful_file.groupby(['Function match', 'Function']).size().unstack(fill_value=0)
            pivot_file = pivot_file.merge(pivot_counts, on="Function match", how="left").fillna(0)
            pivot_file[f"<- these functionalities are based on the functionalities of {tool1}"] = ''
            return pivot_file

        def Alternative_coord(compair_file, tool):
            compair_file[f"{tool}_start"] = np.nan
            compair_file[f"{tool}_end"] = np.nan
            compair_file = compair_file.sort_values(by="End coord")
            compair_file = compair_file.reset_index(drop=True)
            rows_to_remove = []
            for i in range(len(compair_file) - 1):
                current_row = compair_file.iloc[i]
                next_row = compair_file.iloc[i + 1]
                if current_row["End coord"] == next_row["End coord"] and current_row["Start coord"] != next_row["Start coord"]:
                    # If Segment_digger is empty, store alternative coords
                    if pd.notna(current_row["Segment_imgt"]):
                        if pd.isna(current_row[f"Segment_{tool}"]):  # Proper check for NaN
                            compair_file.at[i, f"Segment_{tool}"] = next_row[f"Segment_{tool}"]
                            compair_file.at[i, f"Function_{tool}"] = next_row[f"Function_{tool}"]
                            compair_file.at[i, f"{tool}_start"] = next_row["Start coord"]
                            rows_to_remove.append(i + 1)
                    if pd.notna(next_row["Segment_imgt"]):
                        if pd.isna(next_row[f"Segment_{tool}"]):  # Proper check for NaN
                            compair_file.at[i + 1, f"Segment_{tool}"] = current_row[f"Segment_{tool}"]
                            compair_file.at[i + 1, f"Function_{tool}"] = current_row[f"Function_{tool}"]
                            compair_file.at[i + 1, f"{tool}_start"] = current_row["Start coord"]
                            rows_to_remove.append(i)
            compair_file = compair_file.drop(rows_to_remove).reset_index(drop=True)
            compair_file = compair_file.sort_values(by="Start coord")
            compair_file = compair_file.reset_index(drop=True)
            rows_to_remove = []

            for i in range(len(compair_file) - 1):
                current_row = compair_file.iloc[i]
                next_row = compair_file.iloc[i + 1]
                if current_row["Start coord"] == next_row["Start coord"] and current_row["End coord"] != next_row["End coord"]:
                    if pd.notna(current_row["Segment_imgt"]):
                        if pd.isna(current_row[f"Segment_{tool}"]):  # Proper check for NaN
                            compair_file.at[i, f"Segment_{tool}"] = next_row[f"Segment_{tool}"]
                            compair_file.at[i, f"Function_{tool}"] = next_row[f"Function_{tool}"]
                            compair_file.at[i, f"{tool}_end"] = next_row["End coord"]
                            rows_to_remove.append(i + 1)
                    if pd.notna(next_row["Segment_imgt"]):
                        if pd.isna(next_row[f"Segment_{tool}"]):  # Proper check for NaN
                            compair_file.at[i + 1, f"Segment_{tool}"] = current_row[f"Segment_{tool}"]
                            compair_file.at[i + 1, f"Function_{tool}"] = current_row[f"Function_{tool}"]
                            compair_file.at[i + 1, f"{tool}_end"] = current_row["End coord"]
                            rows_to_remove.append(i)
            compair_file = compair_file.drop(rows_to_remove).reset_index(drop=True)

            for i in range(len(compair_file)):
                current_row = compair_file.iloc[i]
                if pd.notna(current_row[f"Function_{tool}"]) and pd.isna(current_row[f"Segment_{tool}"]):
                    compair_file.at[i, f"Segment_{tool}"] = "Missing segment"
            return compair_file

        compair_vdj_imgt = pd.merge(imgt, vdj, on=["Start coord", "End coord"], how="outer")
        compair_vdj_imgt = Alternative_coord(compair_vdj_imgt,"vdj-insights")
        vdj_2 = compair_vdj_imgt.drop(columns=['Segment_imgt','Function_imgt'])

        compair_vdj_imgt["Short name match"] = compair_vdj_imgt["Segment_imgt"].fillna("") == compair_vdj_imgt["Segment_vdj-insights"].fillna("")
        compair_vdj_imgt["Short name match"] = compair_vdj_imgt.apply(lambda row: check(row, "imgt", "vdj-insights"), axis=1)
        compair_vdj_imgt["Function match"] = compair_vdj_imgt["Function_imgt"].fillna("") == compair_vdj_imgt["Function_vdj-insights"].fillna("")
        compair_vdj_imgt["Function match"] = compair_vdj_imgt.apply(lambda row: check_function(row, 'imgt', 'vdj-insights'), axis=1)
        compair_vdj_imgt.dropna(axis=1, how='all', inplace=True)

        imgt_vdj_data_pivot = compair_vdj_imgt.groupby("Short name match").size().reset_index(name="total_count")
        imgt_vdj_data_pivot = pivot_adds(compair_vdj_imgt, imgt_vdj_data_pivot, "imgt", "vdj-insights")
        imgt_vdj_function_pivot = compair_vdj_imgt.groupby("Function match").size().reset_index(name="total_count")
        imgt_vdj_function_pivot = pivot_function_adds(compair_vdj_imgt, imgt_vdj_function_pivot, "imgt", "vdj-insights")

        compair_digger_imgt = pd.merge(imgt, digger, on=["Start coord", "End coord"], how="outer")
        compair_digger_imgt = Alternative_coord(compair_digger_imgt, "digger")

        digger_2 = compair_digger_imgt.drop(columns=['Segment_imgt','Function_imgt'])
        compair_digger_imgt["Short name match"] = compair_digger_imgt["Segment_imgt"].fillna("") == compair_digger_imgt["Segment_digger"].fillna("")
        compair_digger_imgt["Short name match"] = compair_digger_imgt.apply(lambda row: check(row, "imgt", "digger"), axis=1)
        compair_digger_imgt["Function match"] = compair_digger_imgt["Function_imgt"].fillna("") == compair_digger_imgt["Function_digger"].fillna("")
        compair_digger_imgt["Function match"] = compair_digger_imgt.apply(lambda row: check_function(row, 'imgt', 'digger'), axis=1)
        compair_digger_imgt.dropna(axis=1, how='all', inplace=True)
        imgt_digger_data_pivot = compair_digger_imgt.groupby("Short name match").size().reset_index(name="total_count")
        imgt_digger_data_pivot = pivot_adds(compair_digger_imgt,imgt_digger_data_pivot,"imgt", "digger")
        imgt_digger_function_pivot = compair_digger_imgt.groupby("Function match").size().reset_index(name="total_count")
        imgt_digger_function_pivot = pivot_function_adds(compair_digger_imgt, imgt_digger_function_pivot, "imgt", "digger")

        compair_vdj_digger = pd.merge(vdj_2, digger_2, on=["Start coord", "End coord"], how="outer")
        compair_all = pd.merge(imgt, compair_vdj_digger, on=["Start coord", "End coord"], how="outer")
        compair_vdj_digger["Short name match"] = compair_vdj_digger["Segment_vdj-insights"].fillna("") == compair_vdj_digger["Segment_digger"].fillna("")
        compair_vdj_digger["Short name match"] = compair_vdj_digger.apply(lambda row: check(row, "vdj-insights", "digger"), axis=1)
        compair_vdj_digger["Function match"] = compair_vdj_digger["Function_vdj-insights"].fillna("") == compair_vdj_digger["Function_digger"].fillna("")
        compair_vdj_digger["Function match"] = compair_vdj_digger.apply(lambda row: check_function(row, 'vdj-insights', 'digger'), axis=1)
        compair_vdj_digger.dropna(axis=1, how='all', inplace=True)
        vdj_digger_data_pivot = compair_vdj_digger.groupby("Short name match").size().reset_index(name="total_count")
        vdj_digger_data_pivot = pivot_adds(compair_vdj_digger, vdj_digger_data_pivot, "vdj-insights", "digger")
        vdj_digger_function_pivot = compair_vdj_digger.groupby("Function match").size().reset_index(name="total_count")
        vdj_digger_function_pivot = pivot_function_adds(compair_vdj_digger, vdj_digger_function_pivot, "vdj-insights", "digger")

        merged_df = pd.merge(imgt, vdj_2, on=["Start coord", "End coord"], how="outer")
        compair = pd.merge(merged_df, digger_2, on=["Start coord", "End coord"], how="outer")
        compair = compair[['Start coord', 'End coord', 'Segment_imgt', 'Segment_vdj-insights', 'Segment_digger','vdj-insights_start', 'vdj-insights_end', 'digger_start','digger_end', 'Function_imgt', 'Function_vdj-insights','Function_digger']]
        compair.dropna(axis=1, how='all', inplace=True)

        compair["imgt_vdj-insights"] = ((compair["Segment_imgt"].isna() & compair["Segment_vdj-insights"].isna()) | (compair["Segment_imgt"].notna() & compair["Segment_vdj-insights"].notna()))
        compair["imgt_digger"] = ((compair["Segment_imgt"].isna() & compair["Segment_digger"].isna()) | (compair["Segment_imgt"].notna() & compair["Segment_digger"].notna()))
        compair["vdj-insights_digger"] = ((compair["Segment_vdj-insights"].isna() & compair["Segment_digger"].isna()) | (compair["Segment_vdj-insights"].notna() & compair["Segment_digger"].notna()))
        
        compair["all_true"] = ((compair["Segment_imgt"].notna() & compair["Segment_vdj-insights"].notna() & compair["Segment_digger"].notna()))
        compair["Function_imgt_vdj-insights"] = compair["Function_imgt"].fillna("MISSING") == compair["Function_vdj-insights"].fillna("MISSING")
        compair["Function_imgt_digger"] = compair["Function_imgt"].fillna("MISSING") == compair["Function_digger"].fillna("MISSING")
        compair["Function_vdj-insights_digger"] = compair["Function_vdj-insights"].fillna("MISSING") == compair["Function_digger"].fillna("MISSING")
        compair["Function_all_true"] = (
                (compair["Function_imgt"].fillna("MISSING") == compair["Function_vdj-insights"].fillna("MISSING")) &
                (compair["Function_vdj-insights"].fillna("MISSING") == compair["Function_digger"].fillna("MISSING")))

        digger_file["Segment"] = digger_file["gene_type"].str[-1]
        imgt_pivot = imgt_file.pivot_table(index=['Segment'], columns=['Segment'], aggfunc='size', fill_value=0)
        vdj_insights_pivot = vdj_file.pivot_table(index=['Segment'], columns=['Segment'], aggfunc='size', fill_value=0)
        digger_pivot = digger_file.pivot_table(index=['Segment'], columns=['Segment'], aggfunc='size', fill_value=0)

        digger_table_only = compair[compair["imgt_digger"] == False]
        vdj_table_only = compair[compair["imgt_vdj-insights"] == False]

        os.makedirs(f"{VDJ_IMGT_path}/compair/", exist_ok=True)
        with pd.ExcelWriter(f"{VDJ_IMGT_path}/compair/{sample}_{region}.xlsx") as writer:
            compair.to_excel(writer, sheet_name="compair", index=False)
            imgt_vdj_data_pivot.to_excel(writer, sheet_name="IMGT_VDJ-Insights", index=False)
            imgt_vdj_function_pivot.to_excel(writer, sheet_name="IMGT_VDJ-Insights_function", index=False)
            imgt_digger_data_pivot.to_excel(writer, sheet_name="IMGT_Digger", index=False)
            imgt_digger_function_pivot.to_excel(writer, sheet_name="IMGT_Digger_function", index=False)
            vdj_digger_data_pivot.to_excel(writer, sheet_name="VDJ-Insights_Digger", index=False)
            vdj_digger_function_pivot.to_excel(writer, sheet_name="VDJ-Insights_Digger_function", index=False)
            imgt_pivot.to_excel(writer, sheet_name="Count_IMGT", index=False)
            vdj_insights_pivot.to_excel(writer, sheet_name="Count_VDJ-Insights", index=False)
            digger_pivot.to_excel(writer, sheet_name="Count_Digger", index=False)
            digger_table_only.to_excel(writer, sheet_name="Only_Digger", index=False)
            vdj_table_only.to_excel(writer, sheet_name="Only_VDJ-Insights", index=False)

        #Make ven diagrams

        tools = ["Segment_imgt","Segment_digger","Segment_vdj-insights"]
        compair_all[tools] = compair_all[tools].applymap(lambda x: "segment" if pd.notnull(x) else x)

        imgt_set = set(zip(
            compair_all.loc[compair_all["Segment_imgt"].notna(), "Start coord"],
            compair_all.loc[compair_all["Segment_imgt"].notna(), "End coord"],
            compair_all.loc[compair_all["Segment_imgt"].notna(), "Segment_imgt"],
            [sample] * compair_all["Segment_imgt"].notna().sum()))
        vdj_set = set(zip(
            compair_all.loc[compair_all["Segment_vdj-insights"].notna(), "Start coord"],
            compair_all.loc[compair_all["Segment_vdj-insights"].notna(), "End coord"],
            compair_all.loc[compair_all["Segment_vdj-insights"].notna(), "Segment_vdj-insights"],
            [sample] * compair_all["Segment_vdj-insights"].notna().sum()))
        digger_set = set(zip(
            compair_all.loc[compair_all["Segment_digger"].notna(), "Start coord"],
            compair_all.loc[compair_all["Segment_digger"].notna(), "End coord"],
            compair_all.loc[compair_all["Segment_digger"].notna(), "Segment_digger"],
            [sample] * compair_all["Segment_digger"].notna().sum()))

        plot_venn(imgt_set, digger_set, vdj_set, sample, region, species, f"{VDJ_IMGT_path}/compair/{sample}_{region}_ven-comparison")
        plot_function_ven(compair_all, sample, region, species, f"{VDJ_IMGT_path}/compair/{sample}_{region}_function_comparison")

        #make total Vendiagrams
        if region[0:2] == "IG":
            imgt_IG_set = imgt_IG_set | imgt_set
            vdj_IG_set = vdj_IG_set | vdj_set
            digger_IG_set = digger_IG_set | digger_set
        else:
            imgt_TR_set = imgt_TR_set | imgt_set
            vdj_TR_set = vdj_TR_set | vdj_set
            digger_TR_set = digger_TR_set | digger_set
        imgt_total_set = imgt_total_set | imgt_set
        vdj_total_set = vdj_total_set | vdj_set
        digger_total_set = digger_total_set | digger_set

#plot_venn(imgt_IG_set, digger_IG_set, vdj_IG_set, f"Tool comparison of IG-region in the\n{species}", f"{VDJ_IMGT_path}/compair/{species}_IG_ven-comparison")
#plot_venn(imgt_TR_set, digger_TR_set, vdj_TR_set, f"Tool comparison of TR-region in the\n{species}", f"{VDJ_IMGT_path}/compair/{species}_TR_ven-comparison")
#plot_venn(imgt_total_set, digger_total_set, vdj_total_set, f"Tool comparison in the {species}\n", f"{VDJ_IMGT_path}/compair/{species}_total_ven-comparison")

plot_venn(imgt_IG_set, digger_IG_set, vdj_IG_set, sample="", region="IG", species=species, filename_prefix=f"{VDJ_IMGT_path}/compair/{species}_IG_ven-comparison")
plot_venn(imgt_TR_set, digger_TR_set, vdj_TR_set, sample="", region="TR", species=species, filename_prefix=f"{VDJ_IMGT_path}/compair/{species}_TR_ven-comparison")
plot_venn(imgt_total_set, digger_total_set, vdj_total_set, sample="", region="ALL", species=species, filename_prefix=f"{VDJ_IMGT_path}/compair/{species}_total_ven-comparison")

#combine results in seperated excl file
def read_excel (file,tool1,tool2, sample, region):
    if tool2 == '':
        if tool1 == 'compair':
            df = pd.read_excel(file, sheet_name=tool1)
        else:
            df = pd.read_excel(file, sheet_name=f"Only_{tool1}")
        df_f = None
    else:
        df = pd.read_excel(file, sheet_name=f"{tool1}_{tool2}")
        df_f = pd.read_excel(file,sheet_name=f"{tool1}_{tool2}_function")
        df_f.insert(0, "Sample", sample)
        df_f.insert(0, "Region", region)
    df.insert(0, "Sample", sample)
    df.insert(0, "Region", region)
    return df, df_f

def sum_region_group(df, region_list, total_label):
    """
    Groups the DataFrame by 'Short name match' for the given region_list,
    sums up numeric columns, and assigns a label in the 'Region' column.
    """
    # Filter rows for the region_list
    subset = df[df["Region"].isin(region_list)]
    if subset.empty:
        return pd.DataFrame()  # No rows for these regions

    # Group by 'Short name match' and sum the numeric columns
    if "Short name match" in df:
        grouped = subset.groupby("Short name match", as_index=False)[["total_count", "D_count", "J_count", "V_count"]].sum()
    else:
        grouped = subset.groupby("Function match", as_index=False)[["total_count", "Functional", "ORF", "Pseudo"]].sum()

    # Insert the Region and Sample columns at the front
    grouped.insert(0, "Region", total_label)
    grouped.insert(1, "Sample", "")  # or "All_Samples" if you prefer

    return grouped

compair_dir = os.path.join(VDJ_IMGT_path, "compair")
summary_path = os.path.join(VDJ_IMGT_path, f"compair/{species}_comparison_summary.xlsx")
if os.path.exists(summary_path):
    os.remove(summary_path)

summary_files = [
    f for f in glob.glob(os.path.join(compair_dir, "*.xlsx"))
    if not f.endswith(f"{species}_comparison_summary.xlsx")]



for file in summary_files:
    base = os.path.basename(file)
    # Assume the file is named like "sample_region.xlsx"
    try:
        sample, region_with_ext = base.split('_', 1)
    except ValueError:
        print(f"Skipping file {base}: filename does not match expected pattern 'sample_region.xlsx'")
        continue
    region = region_with_ext.replace('.xlsx', '')
    compair = 'compair'

    compair_df, nope = read_excel(file, compair, '', sample, region)
    compair_list.append(compair_df)
    only_digger_df, nope = read_excel(file, "Digger", '', sample, region)
    only_digger_list.append(only_digger_df)
    only_vdj_df, nope = read_excel(file, "VDJ-Insights", '', sample, region)
    only_vdj_insights_list.append(only_vdj_df)

    imgt_vdj_df, imgt_vdj_f_df = read_excel(file, "IMGT", "VDJ-Insights", sample, region)
    imgt_vdj_list.append(imgt_vdj_df)
    imgt_vdj_f_list.append(imgt_vdj_f_df)
    imgt_digger_df, imgt_digger_f_df = read_excel(file, "IMGT", "Digger", sample, region)
    imgt_digger_list.append(imgt_digger_df)
    imgt_digger_f_list.append(imgt_digger_f_df)
    vdj_digger_df, vdj_digger_f_df = read_excel(file, "VDJ-Insights", "Digger", sample, region)
    vdj_digger_list.append(vdj_digger_df)
    vdj_digger_f_list.append(vdj_digger_f_df)


compair_combined = pd.concat(compair_list, ignore_index=True) if compair_list else None
only_vdj_combined = pd.concat(only_vdj_insights_list, ignore_index=True) if only_vdj_insights_list else None
only_digger_combined = pd.concat(only_digger_list, ignore_index=True) if only_digger_list else None
imgt_vdj_combined = pd.concat(imgt_vdj_list, ignore_index=True) if imgt_vdj_list else None
imgt_digger_combined = pd.concat(imgt_digger_list, ignore_index=True) if imgt_digger_list else None
vdj_digger_combined = pd.concat(vdj_digger_list, ignore_index=True) if vdj_digger_list else None
imgt_vdj_f_combined = pd.concat(imgt_vdj_f_list, ignore_index=True) if imgt_vdj_f_list else None
imgt_digger_f_combined = pd.concat(imgt_digger_f_list, ignore_index=True) if imgt_digger_f_list else None
vdj_digger_f_combined = pd.concat(vdj_digger_f_list, ignore_index=True) if vdj_digger_f_list else None

compair_combined_IG = compair_combined[compair_combined["Region"].str.upper().str.startswith("IG")]
compair_combined_TR = compair_combined[compair_combined["Region"].str.upper().str.startswith("TR")]

#plot_function_ven(compair_combined_IG, f"Functionality comparison IG regions in the {species}", f"{VDJ_IMGT_path}/compair/{species}_IG_ven-function_comparison")
#plot_function_ven(compair_combined_TR, f"Functionality comparison TR regions in the {species}", f"{VDJ_IMGT_path}/compair/{species}_TR_ven-function_comparison")
#plot_function_ven(compair_combined, f"Functionality comparison in the {species}", f"{VDJ_IMGT_path}/compair/{species}_total_ven-function_comparison")

plot_function_ven(compair_combined_IG, sample=None, region="IG", species=species, filename=f"{VDJ_IMGT_path}/compair/{species}_IG_ven-function_comparison")
plot_function_ven(compair_combined_TR, sample=None, region="TR", species=species, filename=f"{VDJ_IMGT_path}/compair/{species}_TR_ven-function_comparison")
plot_function_ven(compair_combined, sample=None, region="all", species=species, filename=f"{VDJ_IMGT_path}/compair/{species}_total_ven-function_comparison")

# For a "total" result, you might want to combine the three main comparisons.
df_dict = {"imgt_vdj_combined": imgt_vdj_combined,
    "imgt_digger_combined": imgt_digger_combined,
    "vdj_digger_combined": vdj_digger_combined}

ig_regions = ["IGH", "IGK", "IGL"]
tr_regions = ["TRA", "TRB", "TRG", "TRA&TRD"]

for name, df in df_dict.items():
    # If the columns might contain NaN, fill with zero first
    df[["total_count", "D_count", "J_count", "V_count"]] = \
        df[["total_count", "D_count", "J_count", "V_count"]].fillna(0)

    # Summation across IG
    df_total_ig = sum_region_group(df, ig_regions, "Total_IG")
    # Summation across TR
    df_total_tr = sum_region_group(df, tr_regions, "Total_TR")
    # Summation across all regions
    df_total = df.groupby("Short name match", as_index=False)[["total_count", "D_count", "J_count", "V_count"]].sum()
    df_total.insert(0, "Region", "Total")
    df_total.insert(1, "Sample", "")

    # Concatenate original + IG total + TR total + overall total
    df_final = pd.concat([df, df_total_ig, df_total_tr, df_total], ignore_index=True)
    df_final.rename(columns={"Short name match": "Comparison"}, inplace=True)

    # 4) Reassign the result back to the dictionary
    df_dict[name] = df_final

imgt_vdj_combined = df_dict["imgt_vdj_combined"]
imgt_digger_combined = df_dict["imgt_digger_combined"]
vdj_digger_combined = df_dict["vdj_digger_combined"]

df_dict = {"imgt_vdj_f_combined": imgt_vdj_f_combined,
    "imgt_digger_f_combined": imgt_digger_f_combined,
    "vdj_digger_f_combined": vdj_digger_f_combined}

ig_regions = ["IGH", "IGK", "IGL"]
tr_regions = ["TRA", "TRB", "TRG", "TRA&TRD"]

for name, df in df_dict.items():
    # If the columns might contain NaN, fill with zero first
    df[["total_count", "Functional", "ORF", "Pseudo"]] = \
        df[["total_count", "Functional", "ORF", "Pseudo"]].fillna(0)

    df_total_ig = sum_region_group(df, ig_regions, "Total_IG")
    df_total_tr = sum_region_group(df, tr_regions, "Total_TR")
    df_total = df.groupby("Function match", as_index=False)[["total_count", "Functional", "ORF", "Pseudo"]].sum()
    df_total.insert(0, "Region", "Total")
    df_total.insert(1, "Sample", "")

    # Concatenate original + IG total + TR total + overall total
    df_final = pd.concat([df, df_total_ig, df_total_tr, df_total], ignore_index=True)

    # 4) Reassign the result back to the dictionary
    df_dict[name] = df_final
# Write all combined DataFrames into a summary Excel file with multiple sheets.
imgt_vdj_f_combined = df_dict["imgt_vdj_f_combined"]
imgt_digger_f_combined = df_dict["imgt_digger_f_combined"]
vdj_digger_f_combined = df_dict["vdj_digger_f_combined"]


summary_filepath = os.path.join(VDJ_IMGT_path, f"compair/{species}_comparison_summary.xlsx")
with pd.ExcelWriter(summary_filepath) as writer:
    #total_combined.to_excel(writer, sheet_name="Total", index=False)
    compair_combined.to_excel(writer, sheet_name="compair_all", index=False)
    imgt_vdj_combined.to_excel(writer, sheet_name="IMGT_VDJ-Insights", index=False)
    imgt_vdj_f_combined.to_excel(writer, sheet_name="IMGT_VDJ-Insights_Functionality", index=False)
    imgt_digger_combined.to_excel(writer, sheet_name="IMGT_Digger", index=False)
    imgt_digger_f_combined.to_excel(writer, sheet_name="IMGT_Digger_Functionality", index=False)
    vdj_digger_combined.to_excel(writer, sheet_name="VDJ-Insights_Digger", index=False)
    vdj_digger_f_combined.to_excel(writer, sheet_name="VDJ-Insights_Digger_Functionality", index=False)
    only_vdj_combined.to_excel(writer, sheet_name="Only_VDJ-Insights", index=False)
    only_digger_combined.to_excel(writer, sheet_name="Only_Digger", index=False)


print(f"Summary file is writen!")
