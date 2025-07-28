import pandas as pd
import numpy as np
import seaborn as sns
import os
import yaml
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted
import re

print ("The code starts!")
with open("config_runs.yaml", "r") as file:
    config = yaml.safe_load(file)

run_paths = [config["run_1"], config["run_2"], config["run_3"]]
run_completest = config["name_run3"]
run_names = [config["name_run1"], config["name_run2"], config["name_run3"]]

output_path = config["output_dir"]  # e.g. "/path/to/output.xlsx"
os.makedirs(output_path, exist_ok=True)

# Start with an empty DataFrame that has the "Start coord" and "End coord" columns
compair_df = pd.DataFrame(columns=["Sample", "Region", "Start coord", "End coord"])

def Alternative_coord(compair_file, tool, position):
    # check-counts
    if position == "Start":
        opposite = "End"
    else:
        opposite = "Start"
    c_end = 0
    c_beide = 0
    compair_file = compair_file.sort_values(by=["Sample", "Region", f"{position} coord", f"{opposite} coord"]).reset_index(drop=True)
    rows_to_remove = []
    for i in range(len(compair_file) - 1):
        current_row = compair_file.iloc[i]
        next_row = compair_file.iloc[i + 1]
        if current_row[f"{position} coord"] == next_row[f"{position} coord"] and current_row[f"{opposite} coord"] != next_row[f"{opposite} coord"]:
            c_end = c_end + 1
            if pd.notna(current_row[f"Segment_{run_completest}"]):
                #print(current_row[f"Segment_{run_completest}"])
                #print(current_row[f"Segment_{tool}"])
                if pd.isna(current_row[f"Segment_{tool}"]):  # Proper check for NaN
                    compair_file.at[i, f"Segment_{tool}"] = next_row[f"Segment_{tool}"]
                    compair_file.at[i, f"Function_{tool}"] = next_row[f"Function_{tool}"]
                    compair_file.at[i, f"{tool}_{opposite}"] = next_row[f"{opposite} coord"]
                    rows_to_remove.append(i + 1)
            #print(next_row[f"Segment_{run_completest}"])
            #print(next_row[f"Segment_{tool}"])
            if pd.notna(next_row[f"Segment_{run_completest}"]):
                if pd.isna(next_row[f"Segment_{tool}"]):  # Proper check for NaN
                    compair_file.at[i + 1, f"Segment_{tool}"] = current_row[f"Segment_{tool}"]
                    compair_file.at[i + 1, f"Function_{tool}"] = current_row[f"Function_{tool}"]
                    compair_file.at[i + 1, f"{tool}_{opposite}"] = current_row[f"{opposite} coord"]
                    rows_to_remove.append(i)
        elif ((next_row[f"{position} coord"] - current_row[f"{position} coord"]) >0 and (next_row[f"{position} coord"] - current_row[f"{position} coord"]) <=5) \
                or (current_row["Region"] == next_row["Region"] and current_row["End coord"] >= next_row["Start coord"]):
            if current_row[f"{opposite} coord"] != next_row[f"{opposite} coord"]:
                c_beide = c_beide + 1
                if pd.notna(current_row[f"Segment_{run_completest}"]):
                    if pd.isna(current_row[f"Segment_{tool}"]):  # Proper check for NaN
                        compair_file.at[i, f"Segment_{tool}"] = next_row[f"Segment_{tool}"]
                        compair_file.at[i, f"Function_{tool}"] = next_row[f"Function_{tool}"]
                        compair_file.at[i, f"{tool}_{position}"] = next_row[f"{position} coord"]
                        compair_file.at[i, f"{tool}_{opposite}"] = next_row[f"{opposite} coord"]
                        rows_to_remove.append(i + 1)
                if pd.notna(next_row[f"Segment_{run_completest}"]):
                    if pd.isna(next_row[f"Segment_{tool}"]):
                        compair_file.at[i + 1, f"Segment_{tool}"] = current_row[f"Segment_{tool}"]
                        compair_file.at[i + 1, f"Function_{tool}"] = current_row[f"Function_{tool}"]
                        compair_file.at[i + 1, f"{tool}_{position}"] = current_row[f"{position} coord"]
                        compair_file.at[i + 1, f"{tool}_{opposite}"] = current_row[f"{opposite} coord"]
                        rows_to_remove.append(i)

    #print(f"count_run = {c_end}")
    #print(f"count_beide = {c_beide}")
    compair_file = compair_file.drop(rows_to_remove).reset_index(drop=True)

    for i in range(len(compair_file)):
        current_row = compair_file.iloc[i]
        if pd.notna(current_row[f"Function_{tool}"]) and pd.isna(current_row[f"Segment_{tool}"]):
            compair_file.at[i, f"Segment_{tool}"] = "Missing segment"
    return compair_file


def pivot_adds(row):
    combi = pd.notna(row["Segment_Combi-lib"])
    imgt = pd.notna(row["Segment_IMGT-lib"])
    vdj = pd.notna(row["Segment_VDJbase-lib"])
    if combi and imgt and vdj:
        return "Annotated using all libraries"
    elif combi and imgt and not vdj:
        return "Annotated using the IMGT library"
    elif combi and vdj and not imgt:
        return "Annotated using the VDJ library"
    else:
        return ''
    '''elif combi and imgt and not vdj:
        return "Annotated using the Combi and IMGT library"
    elif combi and vdj and not imgt:
        return "Annotated using the Combi and VDJbase library"
    elif vdj and imgt and not combi:
        return "Annotated using the IMGT and VDJbase library"
    elif combi and not imgt and not vdj:
        return "Annotated using the Combi library"
    elif vdj and not imgt and not combi:
        return "Annotated using the VDJbase library"
    elif imgt and not vdj and not combi:
        return "Annotated using IMGT library"
    else:
        return "Missing in all libraries"'''

    ful_file["VDJ"] = ful_file.apply(lambda row: extract_vdj(row[f"Segment_{tool1}"])
    if pd.notna(row[f"Segment_{tool1}"]) else extract_vdj(row[f"Segment_{tool2}"]), axis=1)

    pivot_counts = ful_file.groupby(['Short name match', 'VDJ']).size().unstack(fill_value=0)
    pivot_counts.rename(columns={"V": "V_count", "D": "D_count", "J": "J_count"}, inplace=True)
    pivot_file = pivot_file.merge(pivot_counts, on="Short name match", how="left").fillna(0)
    return pivot_file


# Iterate over the runs and names in parallel
for path, name in zip(run_paths, run_names):
    # Build the Excel file path
    excel_path = os.path.join(path, "annotation/annotation_report_all.xlsx")

    # Read the Excel
    run_df = pd.read_excel(excel_path)

    run_df['Region'] = run_df['Path'].str.split("_").str[-1].str.replace(".fasta", "")

    # Subset the columns we need, using double brackets
    run_df.rename(columns={"Target name": f"Segment_{name}", "Function": f"Function_{name}"}, inplace=True)
    subset_df = run_df[["Sample", "Region", "Start coord", "End coord", f"Segment_{name}", f"Function_{name}"]]

    # Merge with our compair_df
    compair_df = pd.merge(compair_df, subset_df, on=["Sample", "Region", "Start coord", "End coord"], how="outer")

#change order columns
cols = compair_df.columns.tolist()
new_cols = ["Sample", "Region", "Start coord", "End coord"] + [c for c in cols if c not in ["Sample", "Region", "Start coord", "End coord"]]
compair_df = compair_df[new_cols]
compair_df['Sample'] = compair_df['Sample'].str.extract(r'^(GCA_\d+\.\d+)')

for name in run_names:
    if name != run_completest:
        compair_df[f"{name}_Start"] = np.nan
        compair_df[f"{name}_End"] = np.nan
        compair_df = Alternative_coord(compair_df, name, "Start")
        compair_df = Alternative_coord(compair_df, name, "End")

compair_df["Function_all_true"] = compair_df[[f"Function_{rn}" for rn in run_names]].fillna("MISSING").nunique(axis=1) == 1

function_cols = [f"Function_{rn}" for rn in run_names]
mask_missing = compair_df[function_cols].isna().any(axis=1)  # True if any column is NaN
compair_df.loc[mask_missing, "Function_all_true"] = np.nan
compair_df.dropna(axis=1, how='all', inplace=True)

#column_order = ["Annotated using all libraries", "Annotated using the Combi and IMGT library", "Annotated using the Combi and VDJbase library", "Annotated using the IMGT and VDJbase library",
 #               "Annotated using the Combi library", "Annotated using IMGT library", "Annotated using the VDJbase library", "Missing in all libraries"]
column_order = ["Annotated using all libraries", "Annotated using the IMGT library", "Annotated using the VDJ library"]

# Generate the library presence column
compair_df["Library_presence"] = compair_df.apply(pivot_adds, axis=1)

# Create the pivot table from the 'Library_presence' column
pivot_table = compair_df.groupby(["Sample", "Library_presence"]).size().reset_index(name="counts")
pivot_table = pivot_table.pivot(index="Sample", columns="Library_presence", values="counts").fillna(0)

# Reorder columns based on the predefined column order
pivot_table = pivot_table[[col for col in column_order if col in pivot_table.columns]]
pivot_table = pivot_table.reset_index()

# Define final column ordering and reorder 'pivot_table'
columns = ["Sample"] + [col for col in column_order if col in pivot_table.columns]
pivot_table = pivot_table[columns]

# Sort 'pivot_table' by the total numeric sum in descending order
pivot_table["Total"] = pivot_table.select_dtypes(include="number").sum(axis=1)
pivot_table = pivot_table.sort_values(by="Total", ascending=True).drop(columns=["Total"])

# Ensure numeric columns are properly typed as float
numeric_columns = pivot_table.select_dtypes(include=[np.number]).columns
pivot_table[numeric_columns] = pivot_table[numeric_columns].astype(float)

# Prepare for plotting
categories = numeric_columns  # Numeric columns are the stacked bar chart's categories
bar_positions = np.arange(len(pivot_table))  # Positions for bars (Samples)
colors = sns.color_palette("Set2", n_colors=len(categories))  # Seaborn palette for colors

# Start plotting a stacked bar chart
plt.figure(figsize=(10, max(2, len(pivot_table) * 0.25)))
left = np.zeros(len(pivot_table))# Dynamic width adjustment
bottom = np.zeros(len(pivot_table))  # Initial bottom for stacking bars

# Plot each category as a stacked bar
for cat, color in zip(categories, colors):
    plt.barh(
        pivot_table["Sample"],           # Y-axis: sample names
        pivot_table[cat],                # Widths of the bars
        label=cat,
        color=color,
        left=left,                       # Stack from the left
        edgecolor="white",
        height=0.8
    )

    left += pivot_table[cat]
    bottom += pivot_table[cat]  # Update bottom for the next category

# Configure the y-axis ticks
plt.yticks(
    ticks=np.arange(len(pivot_table)),
    labels=pivot_table["Sample"],
    rotation=0
)

plt.ylim(-0.5, len(pivot_table) - 0.5)

# Align the X-axis limits perfectly with the first and last bar
max_total = pivot_table[categories].sum(axis=1).max()
plt.xlim(0, max_total * 1.1)  # 10% extra space on the right

# Axis labels
plt.ylabel("Samples")
plt.xlabel("Number of gene segments")
#plt.title("Number of gene segments annotated per library")
plt.legend(title="Library Categories", bbox_to_anchor=(1.01, 1), loc="upper left")

# Fine-tune layout and save the chart
plt.tight_layout()
# plt.show()  # Uncomment if you want to display the chart in the script's output
plt.savefig(f"{output_path}VDJ-Insights_run-comparison_barchart.svg", format="svg", bbox_inches="tight")
plt.savefig(f"{output_path}VDJ-Insights_run-comparison_barchart.pdf", format="pdf", bbox_inches="tight")

compair_df = compair_df.drop(columns=["Library_presence"])

# Write the cleaned data and pivot table to an Excel file
with pd.ExcelWriter(f"{output_path}VDJ-Insights_run-comparison.xlsx") as writer:
    compair_df.to_excel(writer, sheet_name="compair", index=False)  # Save compair_df
    pivot_table.to_excel(writer, sheet_name="library_table", index=False)  # Save pivot_table


print(f"File is written to {output_path}!")