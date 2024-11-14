import pandas as pd
import tempfile
from Bio import AlignIO
from io import StringIO
import shutil
from statistics import multimode
import subprocess

def trim_gaps(alignment):
    """
    Trims the alignment to the region where all sequences have a base, excluding gap positions.
    """
    alignment_length = alignment.get_alignment_length()
    start_pos = None
    end_pos = None

    for i in range(alignment_length):
        bases = [record.seq[i] for record in alignment]
        if all(base != "-" for base in bases):
            start_pos = i
            break

    for i in reversed(range(alignment_length)):
        bases = [record.seq[i] for record in alignment]
        if all(base != "-" for base in bases):
            end_pos = i
            break

    if start_pos is not None and end_pos is not None:
        trimmed_seqs = [record.seq[start_pos:end_pos + 1] for record in alignment]
        return trimmed_seqs
    else:
        return ['' for _ in alignment]

def calculate_overlap(trimmed_alignment):
    overlap_seq = []
    mismatch_count = 0

    for bases in zip(*trimmed_alignment):
        if "-" not in bases:
            if len(set(bases)) == 1:
                overlap_seq.append(bases[0].upper())
            else:
                overlap_seq.append("-")
                mismatch_count += 1
        else:
            overlap_seq.append("-")

    overlap_seq_str = "".join(overlap_seq)
    overlap_percentage = (
        (len(overlap_seq_str) - overlap_seq_str.count("-")) / len(overlap_seq_str)) * 100

    return overlap_percentage, mismatch_count

def run_mafft(input_file):
    mafft_exe = shutil.which("mafft")
    if mafft_exe is None:
        raise FileNotFoundError("MAFFT executable not found.")

    result = subprocess.run([mafft_exe, input_file], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed with the following error: {result.stderr}")
    
    return AlignIO.read(StringIO(result.stdout), "fasta")

def filter_on_alignment(group):
    seqs = list(group["Old name-like seq"])
    with tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.fa') as temp_file:
        for num, seq in enumerate(seqs):
            temp_file.write(f">seq{num+1}\n{seq}\n")
    if align(temp_file.name):
        return group.nlargest(1, 'Length')

def check_groups(group):
    known = group.query("Status == 'Known'")
    segments = group["Segment"].to_list()
    if not known.empty and len(known) == 1:
        return known
    elif len(multimode(segments)) > 1:
        print(segments)
        return group.nlargest(1, 'Length')
    else:
        return filter_on_alignment(group)

def align(input_file):
    alignment = run_mafft(input_file)
    trimmed_alignment = trim_gaps(alignment)
    overlap_percentage, mismatch_count = calculate_overlap(trimmed_alignment)
    return True if overlap_percentage == 100 else False

def assign_group_ids(df):
    """
    Assign group IDs to overlapping or adjacent segments.
    """
    df = df.sort_values(by=["Sample", "Haplotype", "Region", "Start coord"]).reset_index(drop=True)

    group_id = 1
    current_sample = df.loc[0, 'Sample']
    current_start = df.loc[0, 'Start coord']
    current_end = df.loc[0, 'End coord']
    current_haplotype = df.loc[0, 'Haplotype']
    current_region = df.loc[0, 'Region']
    df.loc[0, 'Group'] = group_id

    for i in range(1, len(df)):
        sample = df.loc[i, 'Sample']
        start = df.loc[i, 'Start coord']
        end = df.loc[i, 'End coord']
        haplotype = df.loc[i, 'Haplotype']
        region = df.loc[i, 'Region']
        
        if sample == current_sample and haplotype == current_haplotype and region == current_region and start <= current_end:
            current_end = max(current_end, end)
            df.loc[i, 'Group'] = group_id
        else:
            group_id += 1
            current_sample = sample
            current_start = start
            current_end = end
            current_haplotype = haplotype
            current_region = region
            df.loc[i, 'Group'] = group_id

    return df

def find_non_best_rows(df: pd.DataFrame) -> pd.Index:
    """
    Identify and return the indices of rows that are not the 'best' in overlapping segments.
    'Best' is determined based on the row with the highest count of True values in the boolean columns.
    """
    boolean_columns = ["12_heptamer_matched", "12_nonamer_matched",
                       "13_heptamer_matched", "13_nonamer_matched", 
                       "23_heptamer_matched", "23_nonamer_matched", 
                       "24_heptamer_matched", "24_nonamer_matched"]
    df[boolean_columns] = df[boolean_columns].apply(pd.to_numeric, errors='coerce').fillna(0).astype(bool)
    df = df.sort_values(by=["Sample", "Haplotype", "Region", "Start coord"]).reset_index(drop=True)
    df = assign_group_ids(df)
    non_best_indices = df.groupby(['Haplotype', 'Region', 'Group', 'Sample']).apply(lambda group: select_non_best_rows(group, boolean_columns)).explode()
    return pd.Index(non_best_indices.dropna())

def select_non_best_rows(group: pd.DataFrame, boolean_columns: list) -> pd.Index:
    """
    Select the non-best rows in a group based on the sum of True values in the boolean columns.
    """
    if len(group) > 1:
        group['True_Count'] = group[boolean_columns].sum(axis=1)

        max_true_count = group['True_Count'].max()

        best_group = group.query(f"True_Count == {max_true_count}")
        if len(best_group) > 1:
            try:
                best_group["Length"] = best_group["Old name-like seq"].str.len()
                rechecked = check_groups(best_group)
                return group.index.difference(rechecked.index)
            except:
                return group.index
        else:
            non_best_rows = group[group['True_Count'] != max_true_count]
            return non_best_rows.index
    else:
        return pd.Index([])


def remove_overlapping_segments(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove non-best rows based on overlapping segments and boolean conditions.
    """
    non_best_indices = find_non_best_rows(df)
    df = df.sort_values(
        by=["Sample", "Haplotype", "Region", "Start coord"]).reset_index(drop=True)
    df_cleaned = df.drop(
        non_best_indices, errors='ignore').reset_index(drop=True)
    return df_cleaned
