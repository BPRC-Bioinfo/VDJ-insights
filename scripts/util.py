from pathlib import Path
import pandas as pd
import yaml
import zipfile
import argparse

from logger import custom_logger


logger = custom_logger(__name__)


def make_dir(path: str | Path) -> Path:
    """
    Ensures the specified directory exists by creating it if necessary. Checks if the directory exists,
    and if it doesn't, creates the directory along with any necessary parent directories. Logs the creation
    or existence of the directory. If an error occurs during the directory creation, raises an `OSError`.

    Args:
        path (str): The path of the directory to be created.

    Returns:
        Path: The path of the created or existing directory.

    Raises:
        OSError: If the directory cannot be created due to a system-related error.
    """
    try:
        Path(path).mkdir(parents=True, exist_ok=True)
        logger.info(f"Directory created or already exists: {path}")
    except Exception as e:
        logger.error(f"Failed to create directory {path}: {e}")
        raise OSError(f"Failed to create directory {path}") from e
    return Path(path)


def validate_directory(path: str) -> str:
    """
    Validates that the specified directory exists. Checks if the provided path corresponds to an existing directory.
    If the directory does not exist, raises an `argparse.ArgumentTypeError`.

    Args:
        path (str): Path to the directory to validate.

    Returns:
        str: The validated directory path.

    Raises:
        argparse.ArgumentTypeError: If the directory does not exist.
    """
    if not Path(path).is_dir():
        raise argparse.ArgumentTypeError(
            f"The directory {path} does not exist. "
            "Try another directory!"
        )
    return path


def validate_file(path: str) -> str:
    """
    Validates that the specified file exists. Checks if the provided path corresponds to an existing file.
    If the file does not exist, raises an `argparse.ArgumentTypeError`.

    Args:
        path (str): Path to the file to validate.

    Returns:
        str: The validated file path.

    Raises:
        argparse.ArgumentTypeError: If the file does not exist.
    """
    if not Path(path).is_file():
        raise argparse.ArgumentTypeError(
            f"The file {path} does not exist. Try another file, please!"
        )
    return path


def validate_input(path: str) -> str:
    """
    Validates the input directory, ensuring it exists and contains FASTA files. Checks that the specified directory
    exists and contains at least one FASTA file (with extensions .fasta, .fa, or .fna). If the directory is empty
    or contains no FASTA files, raises an `argparse.ArgumentTypeError`.

    Args:
        path (str): Path to the input directory.

    Returns:
        str: The validated input directory path.

    Raises:
        argparse.ArgumentTypeError: If the directory is empty or contains no FASTA files.
    """

    input_path = Path(path).resolve()
    validate_directory(str(input_path))
    if not any(entry.is_file() for ext in ["*.fasta", "*.fa", "*.fna"] for entry in input_path.glob(ext)):
        raise argparse.ArgumentTypeError(
            f"""The directory {
                input_path} is empty or does not contain any FASTA files!"""
        )
    return str(input_path)


def load_config(path: str | Path) -> dict:
    """
    Loads a configuration file (config.yaml) located in the 'config' directory of the current working directory.

    Args:
        path (str): The current working directory.

    Returns:
        dict: Dictionary containing the configuration settings, including chromosomes of interest and their respective flanking genes.

    Raises:
        Exception: If the configuration file cannot be loaded, logs the error and raises an exception.
    """
    try:
        with open(path) as f:
            config = yaml.safe_load(f)
            logger.info("Config file loaded successfully")
            return config
    except Exception as e:
        logger.error(f"Failed to load config file: {e}")
        raise


def unzip_file(file_path: str | Path, unzip_path: str | Path) -> None:
    """
    Extracts a zip file to the specified directory.

    Args:
        file_path (str or Path): Path to the zip file to be extracted.
        unzip_path (str or Path): Directory where the zip file's contents will be extracted.

    Raises:
        Exception: If the extraction fails, logs the error and raises an exception.
    """
    extract_to_path = Path(unzip_path)
    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to_path)
        logger.info(f'Extracted {file_path} to {extract_to_path}')
    except Exception as e:
        logger.error(f"Failed to extract {file_path} to {extract_to_path}: {e}")
        raise


def seperate_annotation(sample_df: pd.DataFrame, annotation_folder: Path, filename: str):
    """
    Creates separate annotation files for individual samples.

    Args:
        sample_df (pd.DataFrame): Data for a single sample.
        annotation_folder (Path): The directory where the individual report will be saved.
        filename (str): The base file name for saving the report.
    """
    sample = sample_df.iloc[0]['Sample']
    path = annotation_folder / "individual" / sample
    make_dir(path)
    new_file = path / f"{sample}_{filename}"
    sample_df.to_excel(new_file, index=False)
