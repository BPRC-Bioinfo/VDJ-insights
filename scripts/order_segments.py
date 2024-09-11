from mapping import mapping_main
from pathlib import Path
import pandas as pd

from logger import custom_logger

# Method for logger current states of the program.
logger = custom_logger(__name__)


def combine_df(mapping_tools, input_dir, library):
    """
    Calling the mapping script annotation.py. This script is called 
    with either "bowtie", "bowtie2" or "minimap2" as input value. 
    These all return a df.
    Then all df's are concatenated to form a new complete one. 
    Finally dropping all duplicates in the df based on a
    subset of ["start", "stop"], resetting the index of this df 
    and return this new df.

    Returns:
        unique_combinations (DataFrame): A df containing all the unique
        mapping entries from bowtie(2) and minimap2.
    """
    df = pd.DataFrame()
    for tool in mapping_tools:
        mapping_df = mapping_main(tool, input_dir, library, 100, 100)
        df = pd.concat([df, mapping_df])
    unique_combinations = df.drop_duplicates(
        subset=["start", "stop", "haplotype"])
    return unique_combinations.reset_index(drop=True)


def write_report(df, report):
    """
    Creates a excel (xlsx) file of the df and saved as "report.xlsx".

    Args:
        df (DataFrame): df to be saved.
        report (Path): Path of a excel to be saved to.
    """

    df.to_excel(report, index=False)


def get_or_create(annotation_folder, mapping_tool, input_dir, library):
    """
    Verifies if the report.xlsx is present.
    If present it returns the content of the file as a df. Otherwise
    call the combine_df() function to create the df;
    then save it as report.xlsx and return the df.

    Args:
        annotation_folder (Path): Path of the annotation folder.

    Returns:
        df (DataFrame): df containing information from report.xlsx.
    """
    report = annotation_folder / "report.xlsx"
    if not report.exists():
        logger.info("The report.xlsx file does not exist! Creating it!")
        df = combine_df(mapping_tool, input_dir, library)
        write_report(df, report)
        return df
    else:
        return pd.read_excel(report)


def create_annotation_folder(cwd, mapper):
    """Ensure the annotation folder exists."""
    annotation_folder = cwd / f"annotation_{mapper}"
    annotation_folder.mkdir(parents=True, exist_ok=True)
    return annotation_folder


def process_dataframe(mapper, df):
    """Process the DataFrame based on the mapper used."""
    if mapper == "minimap2":
        return df
    elif mapper == "bowtie2":
        return df.query('(stop - start) <= 80')


def prepare_data_frames(library, input_dir, cwd):
    mapping_options = ['minimap2', 'bowtie2']
    dfs = []
    for mapper in mapping_options:
        annotation_folder = create_annotation_folder(cwd, mapper)
        df = get_or_create(annotation_folder, [mapper], input_dir, library)
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['stop'] = pd.to_numeric(df['stop'], errors='coerce')
        processed_df = process_dataframe(mapper, df)
        dfs.append(processed_df)
    return dfs


def concatenate_and_clean(dfs):
    final_df = pd.concat(dfs)
    final_df.drop_duplicates(
        subset=['start', 'stop', 'haplotype'], inplace=True)
    parts = final_df['name'].str.split('_')
    final_df['Status'] = parts.str[-1]
    final_df['Short name'] = parts.str[:-1].apply('_'.join)
    final_df["Sample"] = str(final_df.iloc[0]["reference"]).split("_")[0]
    return final_df


def rename_and_sort_df(df):
    renamed = df.rename(columns={
        'reference': 'Reference',
        'region': "Region",
        'segment': 'Segment',
        'haplotype': 'Haplotype',
        'start': 'Start coord',
        'fasta-file': 'Path'
    })
    renamed.sort_values(by=['Haplotype', 'Region', 'Start coord'], inplace=True)
    return renamed


def filter_data_frame(df, min_distance=5):
    filtered_df = df.groupby(['Haplotype', 'Region']).apply(
        lambda group: group.loc[group['Start coord'].diff().fillna(
            min_distance) >= min_distance]
    ).reset_index(drop=True)
    return filtered_df


def save_data_to_excel(df, file_path):
    df.to_excel(file_path, index=False)
    print(f"DataFrame saved to {file_path}")
    return file_path


def order_main(library, input_dir):
    cwd = Path.cwd()
    dfs = prepare_data_frames(library, input_dir, cwd)
    final_df = concatenate_and_clean(dfs)
    renamed = rename_and_sort_df(final_df)
    filtered_df = filter_data_frame(renamed)
    excel_file = cwd / 'reevaluated.xlsx'
    save_data_to_excel(filtered_df, excel_file)
    return 'reevaluated.xlsx'


if __name__ == "__main__":
    order_main()
