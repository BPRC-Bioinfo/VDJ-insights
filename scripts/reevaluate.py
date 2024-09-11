from pathlib import Path
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO
from order_segments import order_main

from util import make_dir


# Define the current working directory
cwd = Path.cwd()


def read_excel_files(filenames):
    """
    Read multiple Excel files and concatenate them into a single DataFrame.

    Parameters:
    filenames (list of str): List of Excel file names to read.

    Returns:
    pd.DataFrame: A DataFrame containing concatenated data from all files.
    """
    df_list = []
    for filename in filenames:
        try:
            df = pd.read_excel(filename)
            df_list.append(df)
        except Exception as e:
            print(f"Failed to read {filename}: {e}")
    return pd.concat(df_list, ignore_index=True)


def write_library_file(df, output_path):
    """
    Write data from a DataFrame to a FASTA format file.

    Parameters:
    df (pd.DataFrame): DataFrame containing the columns 'Short name', 'Old name-like seq', and 'Status'.
    output_path (Path): The file path to write the FASTA file.
    """
    with open(output_path, "w") as f:
        for _, row in df.iterrows():
            f.write(
                f">{row['Short name']}_{row['Status']}\n{row['Old name-like seq']}\n")


def seqio_dict(path):
    return SeqIO.to_dict(SeqIO.parse(path, "fasta"))


def rename(files: Path, sample_directory: Path) -> None:
    combined = sample_directory / "combined"
    sorted_files = sorted(list(files.glob("*.fasta")))
    if combined.exists():
        shutil.rmtree(combined)
    make_dir(combined)
    for file_path in sorted_files:
        path = Path(file_path)
        name = path.stem
        sequences = seqio_dict(path)
        header = list(sequences.keys())[0]
        sequence = sequences[header].seq
        fasta_out = "_".join(name.split("_")[:-1])
        region_file = combined / f"{fasta_out}.fasta"
        with open(region_file, "a") as w:
            w.write(f">{name}\n{sequence}\n")


def write_read_to_file(read_data, header, seperated):
    """Write read data to a file named after the header."""
    filename = seperated / f"{header}.fasta"
    with open(filename, 'w') as output_file:
        output_file.write(read_data)
    print(f"Written {filename}")


def seperate_fasta(aligned_path, seperated):
    make_dir(seperated)
    """Splits a FASTA file into multiple files based on each read, naming files after the read headers."""
    try:
        with open(aligned_path, 'r') as file:
            read_data = ""
            current_header = ""
            for line in file:
                if line.startswith('>'):
                    if read_data:
                        # Write the current read data to a new file named after the header
                        write_read_to_file(read_data, current_header, seperated)
                        read_data = ""
                    # Remove '>' and any trailing newline/whitespace
                    current_header = line[1:].strip()
                    read_data += line
                else:
                    read_data += line
            # Write the last read to a file
            if read_data:
                write_read_to_file(read_data, current_header, seperated)
    except IOError as e:
        print(f"Error processing file: {e}")


def MAFFT(sample_directory: Path):
    combined = sample_directory / "combined"
    aligned = sample_directory / "aligned"
    seperated = sample_directory / "seperated"
    make_dir(aligned)
    for file in combined.glob("*.fasta"):
        aligened_path = aligned / f"{file.stem}_aligned.fasta"
        if not aligened_path.is_file():
            mafft_command = f"mafft --auto --thread 20 {file} > {aligened_path}"
            subprocess.run(mafft_command, shell=True)
        seperate_fasta(aligened_path, seperated)
    return seperated


def re_evaluate_main(excel_files):
    library_path = cwd / "library.fasta"
    df = read_excel_files(excel_files)
    write_library_file(df, library_path)
    directory = cwd / "region"
    sample_directory = Path(cwd / "reevaluate")
    rename(directory, sample_directory)
    seperated = MAFFT(sample_directory)
    reevaluated_file = order_main(library_path, seperated)
    return str(reevaluated_file)


if __name__ == "__main__":
    re_evaluate_main()
