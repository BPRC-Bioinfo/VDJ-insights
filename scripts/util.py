from pathlib import Path
import pandas as pd
import yaml

from logger import custom_logger


logger = custom_logger(__name__)


def make_dir(dir: str) -> Path:
    """
    Ensures the specified directory exists by creating it if necessary. Checks if the directory exists,
    and if it doesn't, creates the directory along with any necessary parent directories. Logs the creation
    or existence of the directory. If an error occurs during the directory creation, raises an `OSError`.

    Args:
        dir (str): The path of the directory to be created.

    Returns:
        Path: The path of the created or existing directory.

    Raises:
        OSError: If the directory cannot be created due to a system-related error.
    """
    try:
        Path(dir).mkdir(parents=True, exist_ok=True)
        logger.info(f"Directory created or already exists: {dir}")
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")
        raise OSError(f"Failed to create directory {dir}") from e
    return Path(dir)


def load_config(cwd):
    """
    Loads a configuration file (config.yaml) located in the 'config' directory of the current working directory.

    Args:
        cwd (Path): The current working directory.

    Returns:
        dict: Dictionary containing the configuration settings, including chromosomes of interest and their respective flanking genes.

    Raises:
        Exception: If the configuration file cannot be loaded, logs the error and raises an exception.
    """
    try:
        with open(cwd / "config" / "config.yaml") as f:
            config = yaml.safe_load(f)
            logger.info("Config file loaded successfully")
            return config
    except Exception as e:
        logger.error(f"Failed to load config file: {e}")
        raise


def load_config2(cwd):
    """
    Loads a configuration file (config.yaml) located in the 'config' directory of the current working directory.

    Args:
        cwd (Path): The current working directory.

    Returns:
        dict: Dictionary containing the configuration settings, including chromosomes of interest and their respective flanking genes.

    Raises:
        Exception: If the configuration file cannot be loaded, logs the error and raises an exception.
    """
    try:
        with open(cwd) as f:
            config = yaml.safe_load(f)
            logger.info("Config file loaded successfully")
            return config
    except Exception as e:
        logger.error(f"Failed to load config file: {e}")
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
