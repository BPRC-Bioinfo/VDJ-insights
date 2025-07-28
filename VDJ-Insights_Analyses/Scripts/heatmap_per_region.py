import os
import pandas as pd
import matplotlib.pyplot as plt
import re
from Check_dup_haps import find_duplicates_by_region

import itertools
from matplotlib import gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
novel_dict = {}
subgroup_dict = {}


def make_pivot_table(data: pd.DataFrame,output):
    pivot_df = data.pivot_table(
        index=['Sample', 'Region', 'Function'],
        columns='Subgroup',
        aggfunc=lambda x: ', '.join(x),
        values='Allele_status',
        fill_value=''
    ).reset_index()

    # Sort the columns of the pivot table
    sorted_columns = sort_segments(pivot_df.columns)
    pivot_df = pivot_df[sorted_columns]  # Reorder columns in sorted order

    #Reset the index for compatibility
    #pivot_df = pivot_df.reset_index()
    #pivot_df = pivot_df.drop(columns=['index'])

    #pivot_df.to_csv(os.path.join(output, "pivot_genes_samples_collapsed.csv"), sep=',', index=False)
    return pivot_df

def remove_missing_data(pivot_df: pd.DataFrame, missing_threshold: pd.DataFrame):
    incomplete_combinations = missing_threshold[missing_threshold["Complete"] == "Incomplete"]
    incomplete_set = set(zip(incomplete_combinations["Sample"], incomplete_combinations["Region"]))
    cleaned_pivot_df = pivot_df[~pivot_df[["Sample", "Region"]].apply(tuple, axis=1).isin(incomplete_set)]
    return cleaned_pivot_df


def merge_duplicates(pivot_df: pd.DataFrame, duplicates_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each Short_name in duplicates_df, find all (Sample, Region) pairs that share it.
    In pivot_df there should be exactly one row per (Sample, Region, Function).
    Only if *all* Function-rows across those duplicate samples are identical (aside from Sample)
    do we merge them into one row per Function (with Sample=comma-joined list).
    """
    merged_rows = []
    rows_to_drop = set()

    # work on a copy to avoid mutating the input
    df = pivot_df.copy()

    # get the unique Short_name→(Sample,Region) mapping
    dup = duplicates_df[['Short name', 'Sample', 'Region']].drop_duplicates()

    # iterate each clone/short_name
    for shortname, group in dup.groupby('Short name'):
        samples = group['Sample'].unique().tolist()
        regions = group['Region'].unique().tolist()

        # need at least 2 samples and exactly one region
        if len(samples) < 2 or len(regions) != 1:
            continue
        region = regions[0]

        # pull all pivot rows for those samples in that region
        mask = (
            df['Sample'].isin(samples) &
            (df['Region'] == region)
        )
        subset = df[mask]

        # require exactly 3 functions, each appearing once per sample
        func_counts = subset['Function'].value_counts()
        if not (len(func_counts) == 3 and all(cnt == len(samples) for cnt in func_counts)):
            continue

        # check every function-row is identical across samples (aside from Sample/Region)
        cols_to_check = [c for c in df.columns if c not in ('Sample', 'Region')]
        all_match = True
        for func in func_counts.index:
            subf = subset[subset['Function'] == func]
            ref = subf.iloc[0][cols_to_check]
            # .all().all() flattens row & column boolean matrix
            if not subf[cols_to_check].eq(ref).all().all():
                all_match = False
                break

        if not all_match:
            continue

        # merge: one new row per Function, with combined Sample list
        for func in func_counts.index:
            subf = subset[subset['Function'] == func]
            ref = subf.iloc[0][cols_to_check]
            new_row = ref.to_dict()
            new_row['Sample'] = ', '.join(sorted(samples))
            new_row['Region'] = region
            merged_rows.append(new_row)

        # mark originals for removal
        rows_to_drop.update(subset.index)

    # drop originals and append merged rows
    df = df.drop(index=rows_to_drop)
    if merged_rows:
        df = pd.concat([df, pd.DataFrame(merged_rows)], ignore_index=True, sort=False)

    # restore sort order
    df = df.sort_values(['Sample', 'Region', 'Function']).reset_index(drop=True)
    return df

def sort_segments(segment_labels):
    """
    Sorts segments by:
      1. Type priority J(1) > D(2) > V(3)
      2. For V-segments only: number after first ‘-’ (ascending),
         then number after second ‘-’ (variant, ascending).
      3. For J/D: family number (digits after letter) ascending,
         then number after first ‘-’ (variant, ascending).
    """
    def parse_key(label):
        # 1) segment‐type priority
        m = re.search(r"([JDV])", label)
        seg = m.group(1) if m else ''
        prio = {"J": 1, "D": 2, "V": 3}.get(seg, 99)

        # pull out all hyphen‐numbers
        hyphens = re.findall(r"-(\d+)", label)
        if seg == "V":
            # for V: primary = first hyphen number, variant = second (if any)
            primary = int(hyphens[0]) if hyphens else float("inf")
            variant = int(hyphens[1]) if len(hyphens) > 1 else 0
            return (prio, primary, variant)
        else:
            # for J/D: primary = family (digits after J/D), variant = first hyphen
            fam_m = re.search(rf"{seg}(\d+)", label)
            fam = int(fam_m.group(1)) if fam_m else float("inf")
            variant = int(hyphens[0]) if hyphens else 0
            return (prio, fam, variant)

    return sorted(segment_labels, key=parse_key)

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

def reorder_matrix(matrix):
    row_linkage = linkage(matrix, method='average', metric='euclidean')
    row_order = leaves_list(row_linkage)
    reordered_matrix = matrix.iloc[row_order, :]
    reordered_rows = matrix.index[row_order].tolist()
    return reordered_matrix, reordered_rows

def clean_segment_name(name):
    match = re.match(r"([A-Z]+(?:\(\w+\))?)(\d+)?(?:-.+?)?(?:S\d+)?$", name)
    if match:
        gene = match.group(1)
        number = match.group(2)
        if number:
            return f"{gene}{number}"
        return gene
    return name

'''
Other color options
color_viridis = rgb_known: tuple = (0.45, 0.30, 0.55), rgb_intermediate: tuple = (0.60, 0.75, 0.50), rgb_novel: tuple = (0.993248, 0.906157, 0.143936)
color_inferno = rgb_known: tuple = (0.35, 0.30, 0.25), rgb_intermediate: tuple = (0.70, 0.45, 0.30), rgb_novel: tuple = (0.99, 0.98, 0.64)
color_inferno2 = rgb_known: tuple = (0.50, 0.45, 0.40), rgb_intermediate: tuple = (0.75, 0.60, 0.50), rgb_novel: tuple = (0.92, 0.88, 0.55)
color_cool-warm = rgb_known: tuple = (0.55, 0.75, 0.85), rgb_intermediate: tuple = (0.98, 0.90, 0.55), rgb_novel: tuple = (0.95, 0.55, 0.55)
color_blue_scale = rgb_known: tuple = (0.85, 0.92, 1.0), rgb_intermediate: tuple = (0.55, 0.70, 0.90), rgb_novel: tuple = (0.30, 0.50, 0.75)
color_green_red = rgb_known: tuple = (0.6, 0.2, 0.2), rgb_intermediate: tuple = (0.9, 0.9, 0.2), rgb_novel: tuple = (0.2, 0.6, 0.2)
original_colors = rgb_known: tuple = (0.57, 0.5, 0.120), rgb_intermediate: tuple = (0.55, 0.25, 0.65), rgb_novel: tuple = (0.0, 0.5, 0.5)
'''
#def create_matrix(data: pd.DataFrame, region: str, output: str,
 #                 rgb_known: tuple = (0.57, 0.5, 0.120), rgb_novel: tuple = (0.0, 0.5, 0.5)) -> None:
def create_matrix(data: pd.DataFrame, region: str, output: str,
                  rgb_known: tuple = (0.80, 0.00, 0.00), rgb_intermediate: tuple = (0.98, 0.90, 0.55), rgb_novel: tuple = (0.00, 0.55, 1.00)) -> None:
    # Filter data simply by region (Function no longer a part of filtering)
    data_total = data[data['Region'] == region]
    data_total = data_total.sort_values(by=['Sample'])
    data_total = data_total.set_index('Sample')
    filtered_data = data_total.drop(columns=["Function"])
    filtered_data = filtered_data.groupby("Sample").agg(lambda x: ', '.join(filter(None, x)) if len(list(filter(None, x))) > 1 else (
            list(filter(None, x))[0] if len(list(filter(None, x))) == 1 else ''))


    # Get references and presence-absence matrix
    references = filtered_data.columns[1:]
    samples = filtered_data.index.values  # Save original sample names
    presence_absence_matrix = filtered_data[references].copy()

    # Map to categorical codes (1: Known, 2: Novel, 0: Empty)
    def map_to_code(value):
        if isinstance(value, str):
            value_set = set(value.split(','))
            known_count = sum(1 for v in value_set if "_Known" in v)
            novel_count = sum(1 for v in value_set if "_Novel" in v)
            total_count = known_count + novel_count
            if novel_count >= 1:
                return 2 - known_count / total_count, total_count
            elif known_count >= 1:
                return known_count / total_count, total_count
            else:
                return 0, 0
        return 0, 0

    cm_coded = presence_absence_matrix.applymap(lambda x: map_to_code(x)[0])
    cm_coded.index = samples  # Ensure original sample names are preserved in the index
    count_matrix = presence_absence_matrix.applymap(lambda x: map_to_code(x)[1])
    count_matrix.index = samples

    # --- Clean Matrix ---
    row_is_empty = cm_coded.sum(axis=1) == 0
    col_is_empty = cm_coded.sum(axis=0) == 0

    cm_coded = cm_coded.loc[~row_is_empty, ~col_is_empty]
    count_matrix = count_matrix.loc[~row_is_empty, ~col_is_empty]
    presence_absence_matrix = presence_absence_matrix.loc[~row_is_empty, ~col_is_empty]

    if cm_coded.empty:
        print(f"No data remains for region '{region}' after cleaning. Skipping heatmap generation.")
        return

    remaining_samples = cm_coded.index.tolist()
    remaining_references = cm_coded.columns.tolist()
    num_samples = len(remaining_samples)
    num_references = len(remaining_references)
    cm = presence_absence_matrix.values

    # Reorder Matrix
    cm_coded, reordered_samples = reorder_matrix(cm_coded)
    count_matrix = count_matrix.loc[reordered_samples, :].reset_index()
    cm_coded = cm_coded.replace(0.0, np.nan)

    #Make Fig
    fig = plt.figure(figsize=(num_references * 0.35 + 2, num_samples * 0.35 + 2))
    gs = gridspec.GridSpec(2, 2, width_ratios=[5, 0.5], height_ratios=[0.5, 5], hspace=0.02, wspace=0.01)
    colors = [plt.cm.tab20(i / 20) for i in range(20)]
    ax_heatmap = plt.subplot(gs[1, 0])

    # Making upperbar
    columns_to_change = [col for col in data_total.columns if col not in ['Sample', 'Function', 'Region']]
    data_total[columns_to_change] = data_total[columns_to_change].applymap(lambda x: len(x.split(',')) if isinstance(x, str) and x.strip() != '' else 0)

    unique_functionalities = data_total["Function"].unique()
    functionality_colors = dict(zip(unique_functionalities, colors[::-1]))

    ax_bar_x = plt.subplot(gs[0, 0], sharex=ax_heatmap)
    #bar_positions = np.arange(num_references)

    functionality_sums = data_total.groupby('Function')[remaining_references].sum()
    functionalities = functionality_sums.index
    bottom = np.zeros(len(remaining_references))

    for functionality in (functionalities):
        values = functionality_sums.loc[functionality].values
        #bars_x = ax_bar_x.bar(bar_positions, values, bottom=bottom, color=functionality_colors[functionality], label=functionality)
        bars_x = ax_bar_x.bar(range(num_references), values, bottom=bottom, color=functionality_colors[functionality], label=functionality)
        bottom += values
    for i, total in enumerate(bottom):
        ax_bar_x.text(i, total, int(total), ha='center', va='bottom', fontsize=10)
    ax_bar_x.set_xlim(-0.5, num_references - 0.5)
    ax_bar_x.axis('off')

    real_num_samples = sum(len(re.split(r'\s*,\s*', s)) for s in remaining_samples)

    handles, labels = ax_bar_x.get_legend_handles_labels()
    #ax_bar_x.legend(handles, labels, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower right", prop={"size": 10}, ncol=len(labels)).set_title("Functionality")
    ax_bar_x.set_title(f'Region: {region} (N={real_num_samples})',  fontsize=20, y=1.3, weight='bold', loc='left')

    # Side barchard
    ax_bar_y = plt.subplot(gs[1, 1], sharey=plt.subplot(gs[1, 0]))

    reference_labels = presence_absence_matrix.columns.values
    segment_names = list(set([re.sub(r'[A-Z]$', '', re.split(r'[-*]', col)[0]) for col in reference_labels]))
    segment_names = [clean_segment_name(name) for name in segment_names]
    segment_names = sorted(segment_names, key=natural_sort_key)
    segment_names = list(dict.fromkeys(segment_names))
    bottom = np.zeros(num_samples)

    distinct_colors = plt.cm.tab20.colors
    color_cycle = itertools.cycle(distinct_colors)

    for segment_name in segment_names:
        # Identify all columns in reference_labels matching this segment_name
        matching_columns = [col for col in reference_labels if re.split(r'[-*]', col)[0] == segment_name]
        if not matching_columns:  # Skip if no matching columns
            continue
        values = count_matrix[matching_columns].sum(axis=1)
        segment_color = next(color_cycle)
        ax_bar_y.barh(range(num_samples), values, left=bottom, color=segment_color, label=segment_name)
        bottom += values

    # --- Plotting the Heatmap ---
    #custom_cmap = LinearSegmentedColormap.from_list("Known-novel-gradient", [rgb_known, rgb_novel], N=256)
    custom_cmap = LinearSegmentedColormap.from_list("Known-novel-gradient", [rgb_known, rgb_intermediate, rgb_novel], N=256)

    im = ax_heatmap.imshow(cm_coded, interpolation='nearest', cmap=custom_cmap, aspect='auto')
    count_matrix = count_matrix.iloc[:, 1:]

    # Set labels for heatmap axes
    ax_heatmap.set_xticks(np.arange(num_references))
    ax_heatmap.set_xticklabels(remaining_references, rotation=90, ha='center', fontsize=14)
    ax_heatmap.set_yticks(np.arange(num_samples))
    ax_heatmap.set_yticklabels(reordered_samples, fontsize=14)
    ax_heatmap.set_xticks(np.arange(-.5, num_references, 1), minor=True)
    ax_heatmap.set_yticks(np.arange(-.5, num_samples, 1), minor=True)
    #ax_heatmap.tick_params(which='minor', bottom=False, left=False)

    ax_heatmap.grid(which='minor', color='w', linestyle='-', linewidth=4)
    ax_heatmap.tick_params(which='minor', bottom=False, left=False)

    ax_heatmap.set_ylabel('Samples', fontsize=14)
    ax_heatmap.set_xlabel('Segments', fontsize=14)

    for i in range(len(count_matrix.index)):
        for j in range(len(count_matrix.columns)):
            if count_matrix.iloc[i, j] > 1:
                total_count = count_matrix.iloc[i, j]
                if total_count > 1:
                    ax_heatmap.text(j, i, int(total_count), ha='center', color='Black', va='center', fontsize=10)

    # Add Legends
    handles_upper, labels_upper = ax_bar_x.get_legend_handles_labels()  # Upper barplot
    handles_right, labels_right = ax_bar_y.get_legend_handles_labels()  # Right-side barplot

    # First legend for "Functionality" (upper barplot)
    legend_func = fig.legend(
        handles_upper, labels_upper,
        bbox_to_anchor=(1, 1),  # Place in top-right corner
        loc="upper right",
        prop={"size": 10},  # Font size
        title="Functionality",  # Title for legend box
        ncol=3)

    # Draw canvas to update positions
    fig.canvas.draw()
    box_func = legend_func.get_window_extent().transformed(fig.transFigure.inverted())

    # ---- Add the second legend for "Segment Families" ----
    # Define position below the first legend
    offset_func_fam = 0.01  # Offset between functionality and segment legend
    legend_box_x, legend_box_y1 = box_func.x0, box_func.y0
    legend_seg_bottom = legend_box_y1 - offset_func_fam

    # Create the second legend for "Segment Families"
    legend_seg = fig.legend(
        handles_right, labels_right,
        bbox_to_anchor=(legend_box_x, legend_seg_bottom),  # Place directly below first legend
        loc="upper left",
        prop={"size": 10},
        title="Segment Families",
        ncol=3)

    ax_bar_y.spines['top'].set_visible(False)
    ax_bar_y.spines['bottom'].set_visible(False)
    ax_bar_y.spines['left'].set_visible(False)
    ax_bar_y.spines['right'].set_visible(False)

    ax_bar_y.yaxis.set_visible(False)
    ax_bar_y.set_xlabel('Number of segments', fontsize=14, loc='left')

    # Draw canvas again to update positions for the second legend
    fig.canvas.draw()
    # Get bounding box of the Segment Families legend
    box_seg = legend_seg.get_window_extent().transformed(fig.transFigure.inverted())

    # ---- Add the colorbar below both legends ----
    offset_fam_cbar = 0.01  # Offset between segment legend and colorbar
    cbar_height = 0.02  # Height of the colorbar

    # Calculate the bottom of the colorbar
    cbar_bottom = box_seg.y0 - offset_fam_cbar - cbar_height

    # Create an axes for the colorbar and place it below the Segment Families legend
    cbar_ax = fig.add_axes([
        box_seg.x0,  # Align with the left-hand side of the legends
        cbar_bottom,  # Bottom of the colorbar
        box_seg.width,  # Same width as the legends
        cbar_height])  # Height of the colorbar


    # Add colorbar to the figure
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([1, 2])  # Set tick positions for the colorbar
    cbar.set_ticklabels(["Known", "Novel"])  # Label the ticks
    cbar.ax.tick_params(labelsize=12)

    # Save the figure
    plt.savefig(f"{output}/{region}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output}/{region}.svg", bbox_inches='tight')

    plt.close(fig)



def main(path, species, methods, region) -> None:
    data = pd.read_excel(f"{path}/{methods}/{region}/annotation/annotation_report_all.xlsx")
    #Not necessary.
    #data['Sample'] = data['Sample'].str.split('_').str[0]
    completeness = pd.read_excel(f"{path}figure/broken_regions/{species}_{methods}_completeness.xlsx")
    duplicates = find_duplicates_by_region(data)

    data['Short name strip'] = data['Short name'].str.replace(r'-like*', '', regex=True)
    data['Subgroup'] = data['Short name strip'].str.split('*').str[0]
    data['Allele'] = data['Short name strip'].str.split('*').str[-1]
    data = data.sort_values(by=['Status', 'Subgroup', 'Allele']).reset_index(drop=True)

    for idx, row in data.iterrows():
        if row['Status'] == 'Novel':
            key = (row["Subgroup"], row["Allele"])
            subgroup = row["Subgroup"]
            if key not in novel_dict:
                if subgroup not in subgroup_dict:
                    subgroup_dict[subgroup] = 1
                else:
                    subgroup_dict[subgroup] += 1
                novel_dict[key] = str(subgroup_dict[subgroup]) + "_Novel"
            data.at[idx, "Allele_status"] = novel_dict[key]
        else:
            data.at[idx, "Allele_status"] = row["Allele"] + "_Known"

    data['Allele'] = data['Allele_status'].str.split('_').str[0].apply(lambda x: str(int(x)) if x.isdigit() else x)
    data['Allele_status'] = data['Allele'] + "_" + data['Status']

    output = f"{path}figure/Heatmap/"
    os.makedirs(output, exist_ok=True)

    pivot_df = make_pivot_table(data, output)

    #missing_threshold = 0.9
    column_order = ['Sample', 'Region', 'Function'] + [col for col in pivot_df.columns if col not in ['Sample', 'Region', 'Function']]
    pivot_df = pivot_df[column_order]
    pivot_df = remove_missing_data(pivot_df, completeness)
    duplicates = remove_missing_data(duplicates, completeness)
    pivot_df = merge_duplicates(pivot_df, duplicates)

    for region in pivot_df['Region'].unique():
        create_matrix(pivot_df, region=region, output=output)

    print("We made the Figures!!")


if __name__ == '__main__':
    species = "homo_sapiens" #rhesus_macaque is also possible.
    methods = "scaffolds" #contigs or scaffods
    region = "bcr" #bcr or tcr
    path = "/path/to/outcomes/HPGC"
    main(path, species, methods, region)
