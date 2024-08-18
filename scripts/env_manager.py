import subprocess
from pathlib import Path
import os
import shutil
import filecmp
from logger import custom_logger

# Initialize the logger
logger = custom_logger(__name__)

# Default directory to save the environments
DEFAULT_ENV_ROOT_DIR = Path('.tool/conda')
ARCHIVE_DIR = Path('.tool/conda_archive')
SAVED_ENV_YAML_DIR = Path('.tool/saved_envs')


def create_and_activate_env(env_file, env_root_dir=DEFAULT_ENV_ROOT_DIR, saved_env_yaml_dir=SAVED_ENV_YAML_DIR):
    """
    Creates and activates a conda environment from a given environment YAML file.
    If the environment already exists and matches the YAML file, it is simply activated.
    If the environment has changed, the old one is archived, and a new environment is created.

    Args:
        env_file (Path): The path to the environment YAML file.
        env_root_dir (Path, optional): The root directory where environments are stored. Defaults to DEFAULT_ENV_ROOT_DIR.
        saved_env_yaml_dir (Path, optional): The directory where environment YAML files are saved for future comparison. Defaults to SAVED_ENV_YAML_DIR.

    Returns:
        Path: The full path to the created or activated conda environment, or None if creation failed.
    """
    env_name = env_file.stem
    env_dir = env_root_dir / env_name
    saved_env_yaml_file = saved_env_yaml_dir / f"{env_name}.yaml"

    env_root_dir.mkdir(parents=True, exist_ok=True)
    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    saved_env_yaml_dir.mkdir(parents=True, exist_ok=True)

    if env_dir.exists() and saved_env_yaml_file.exists():
        if filecmp.cmp(env_file, saved_env_yaml_file, shallow=False):
            logger.environment(
                f"Environment {env_name} is up to date. Activating it.")
            activate_env(env_dir, env_name)
            return env_dir
        else:
            logger.environment(
                f"Environment {env_name} has changed. Recreating it.")
            archive_path = ARCHIVE_DIR / env_name
            if archive_path.exists():
                shutil.rmtree(archive_path)
            shutil.move(env_dir, archive_path)

    logger.environment(f"Creating environment {env_name}.")
    result = subprocess.run(
        ["conda", "env", "create", "--file",
            str(env_file), "--prefix", str(env_dir)],
        text=True
    )
    if result.returncode != 0:
        logger.error(f"Error creating environment: {result.stderr}")
        return None

    shutil.copy(env_file, saved_env_yaml_file)

    logger.environment(f"Activating environment {env_name}.")
    activate_env(env_dir, env_name)

    return env_dir


def activate_env(env_path, env_name):
    """
    Activates the specified conda environment by modifying the PATH and CONDA_PREFIX environment variables.

    Args:
        env_path (Path): The full path to the conda environment to be activated.
        env_name (str): The name of the conda environment to be activated.
    """
    conda_prefix = str(env_path)
    os.environ["CONDA_PREFIX"] = conda_prefix
    os.environ["PATH"] = f"{conda_prefix}/bin:" + os.environ["PATH"]
    logger.environment(f"Activated environment {env_name} at {env_path}")


def deactivate_env():
    """
    Deactivates the current conda environment by removing its paths from the PATH environment variable and unsetting CONDA_PREFIX.
    """
    conda_prefix = os.environ.pop("CONDA_PREFIX", None)
    if conda_prefix:
        env_name = Path(conda_prefix).name
        os.environ["PATH"] = ":".join(p for p in os.environ["PATH"].split(
            ":") if not p.startswith(conda_prefix))
        logger.environment(
            f"Deactivated environment {env_name} at {conda_prefix}")
