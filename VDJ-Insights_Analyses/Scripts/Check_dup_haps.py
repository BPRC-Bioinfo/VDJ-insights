import pandas as pd

def find_duplicates_by_region(df):
    """
    Finds duplicate 'Short name' values within each region in the provided DataFrame,
    ignoring the order of items in the 'Short name' strings.

    Parameters:
        df (pd.DataFrame): The input DataFrame. It must include the columns 'Region', 'Sample', and 'Short name'.

    Returns:
        pd.DataFrame: A DataFrame containing the duplicate rows within each region.
    """
    # Ensure the required columns are present
    required_columns = ["Region", "Sample", "Short name"]
    if not all(column in df.columns for column in required_columns):
        raise KeyError(f"The required columns {required_columns} are not in the DataFrame.")

    # Group the data by Region and Sample, and gather the "Short name" values as sorted tuples
    grouped = df.groupby(["Region", "Sample"])["Short name"].apply(lambda x: tuple(sorted(set(x)))).reset_index()

    # Find duplicates within each region
    duplicates = grouped.groupby("Region").apply(
        lambda group: group[group.duplicated(subset="Short name", keep=False)]
    ).reset_index(drop=True)

    if not duplicates.empty:
        print("The following Samples share exactly the same segments within the same region:")
        print(duplicates)
    else:
        print("No Samples have the exact same segments within any region.")

    return duplicates
'''
file_path = "/Users/ott/PycharmProjects/VDJ-insights/Digger/output/annotation_report_all_BCR.xlsx"  # Your file path
sheet_name = "Sheet1"  # Replace with the correct sheet name if needed

# Read the Excel file into a DataFrame
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Use the function to find duplicates
duplicates = find_duplicates_by_region(df)

# Use or print the duplicates as needed
if not duplicates.empty:
    print("The following Samples share exactly the same segments within the same region:")
    print(duplicates)
else:
    print("No Samples have the exact same segments within any region.")'''



