import pandas as pd
import tempfile
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from io import StringIO
import shutil


def trim_gaps(alignment):
    alignment_length = alignment.get_alignment_length()
    start_pos = next(
        (i for i in range(alignment_length) if sum(
            base != "-" for base in [record.seq[i] for record in alignment]) >= 2),
        None
    )
    end_pos = next(
        (i for i in reversed(range(alignment_length)) if sum(
            base != "-" for base in [record.seq[i] for record in alignment]) >= 2),
        None
    )
    return [record.seq[start_pos:end_pos + 1] for record in alignment]


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

    mafft_cline = MafftCommandline(mafft_exe, input=input_file)
    stdout, stderr = mafft_cline()
    return AlignIO.read(StringIO(stdout), "fasta")


def filter_on_alignment(group):
    seqs = list(group["Old name-like seq"])
    with tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.fa') as temp_file:
        for num, seq in enumerate(seqs):
            temp_file.write(f">seq{num+1}\n{seq}\n")
    if align(temp_file.name):
        return group.nlargest(1, 'Length')


def check_groups(group):
    known = group.query("Status == 'Known'")
    if not known.empty and len(known) == 1:
        return known
    else:
        return filter_on_alignment(group)


def align(input_file):
    alignment = run_mafft(input_file)
    trimmed_alignment = trim_gaps(alignment)
    overlap_percentage, mismatch_count = calculate_overlap(trimmed_alignment)
    return True if overlap_percentage == 100 else False


def find_non_best_rows(df: pd.DataFrame) -> pd.Index:
    """
    Identify and return the indices of rows that are not the 'best' in overlapping segments.
    'Best' is determined based on the row with the highest count of True values in the boolean columns.
    """
    boolean_columns = ["12_heptamer_matched", "12_nonamer_matched",
                       "23_heptamer_matched", "23_nonamer_matched"]

    # Ensure boolean columns are correctly handled
    df[boolean_columns] = df[boolean_columns].apply(
        pd.to_numeric, errors='coerce').fillna(0).astype(bool)

    # Sort the DataFrame by Haplotype, Region, and Start coord
    df = df.sort_values(
        by=["Haplotype", "Region", "Start coord", "Sample"]).reset_index(drop=True)

    # Create a column to track overlap groupings
    df['Group'] = (df['Start coord'] > df['End coord'].shift()).cumsum()

    # Select non-best rows for each group
    non_best_indices = df.groupby(['Haplotype', 'Region', 'Group', 'Sample']).apply(
        lambda group: select_non_best_rows(group, boolean_columns)).explode()

    # Return as a pandas Index object
    return pd.Index(non_best_indices.dropna())


def select_non_best_rows(group: pd.DataFrame, boolean_columns: list) -> pd.Index:
    """
    Select the non-best rows in a group based on the sum of True values in the boolean columns.
    """
    # Calculate True count for each row
    group['True_Count'] = group[boolean_columns].sum(axis=1)

    # Determine the maximum True count
    max_true_count = group['True_Count'].max()

    best_group = group.query(f"True_Count == {max_true_count}")
    if len(best_group) > 1:
        best_group["Length"] = best_group["Old name-like seq"].str.len()
        rechecked = check_groups(best_group)
        return group.index.difference(rechecked.index)
    non_best_rows = group[group['True_Count'] != max_true_count]
    return non_best_rows.index


def remove_non_best_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove non-best rows based on overlapping segments and boolean conditions.
    """
    non_best_indices = find_non_best_rows(df)
    df = df.sort_values(
        by=["Haplotype", "Region", "Start coord", "Sample"]).reset_index(drop=True)
    df_cleaned = df.drop(
        non_best_indices, errors='ignore').reset_index(drop=True)

    return df_cleaned


def remove_overlapping_segments(df):
    return remove_non_best_rows(df)
