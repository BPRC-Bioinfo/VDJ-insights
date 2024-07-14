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
    Create and activate a conda environment from an environment YAML file.

    Parameters:
    env_file (Path): Path to the environment YAML file.
    env_root_dir (Path): Root directory to save the environments.
    saved_env_yaml_dir (Path): Directory to save the environment YAML files.

    Returns:
    Path: Full path to the created or activated environment.
    """
    env_name = env_file.stem  # Get the environment name from the file name
    env_dir = env_root_dir / env_name  # Full path to the environment
    saved_env_yaml_file = saved_env_yaml_dir / f"{env_name}.yaml"

    # Ensure the necessary directories exist
    env_root_dir.mkdir(parents=True, exist_ok=True)
    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    saved_env_yaml_dir.mkdir(parents=True, exist_ok=True)

    if env_dir.exists() and saved_env_yaml_file.exists():
        # Compare the saved YAML with the current one
        if filecmp.cmp(env_file, saved_env_yaml_file, shallow=False):
            logger.environment(
                f"Environment {env_name} is up to date. Activating it.")
            # Environment is the same, just activate it
            activate_env(env_dir, env_name)
            return env_dir
        else:
            logger.environment(
                f"Environment {env_name} has changed. Recreating it.")
            # Archive the existing environment
            archive_path = ARCHIVE_DIR / env_name
            if archive_path.exists():
                shutil.rmtree(archive_path)
            shutil.move(env_dir, archive_path)

    # Create the environment if it does not exist or has changed
    logger.environment(f"Creating environment {env_name}.")
    result = subprocess.run(
        ["conda", "env", "create", "--file",
            str(env_file), "--prefix", str(env_dir)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if result.returncode != 0:
        logger.error(f"Error creating environment: {result.stderr}")
        return None

    # Save a copy of the current YAML file
    shutil.copy(env_file, saved_env_yaml_file)

    logger.environment(f"Activating environment {env_name}.")
    # Activate the environment
    activate_env(env_dir)

    return env_dir


def activate_env(env_path, env_name):
    """
    Activate the specified conda environment.

    Parameters:
    env_path (Path): Full path to the conda environment.
    """
    conda_prefix = str(env_path)
    os.environ["CONDA_PREFIX"] = conda_prefix
    os.environ["PATH"] = f"{conda_prefix}/bin:" + os.environ["PATH"]
    logger.environment(f"Activated environment {env_name} at {env_path}")


def deactivate_env():
    """
    Deactivate the current conda environment.
    """
    conda_prefix = os.environ.pop("CONDA_PREFIX", None)
    env_name = Path(conda_prefix).name
    if conda_prefix:
        # Reset the PATH to its original state by removing the conda path
        os.environ["PATH"] = ":".join(p for p in os.environ["PATH"].split(
            ":") if not p.startswith(conda_prefix))
        logger.environment(
            f"Deactivated environment {env_name} at {conda_prefix}")
