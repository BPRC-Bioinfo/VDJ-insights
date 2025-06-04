"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

import subprocess
from pathlib import Path
import os
import shutil
import filecmp

from .logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def get_archive_dir():
    return Path.cwd() / '.tool' / 'conda_archive'


def create_and_activate_env(env_file, env_root_dir=None, saved_env_yaml_dir=None):
    """
    Creates and activates a conda environment from a given environment YAML file.
    If the environment already exists and matches the YAML file, it is simply activated.
    If the environment has changed, the old one is archived, and a new environment is created.

    Args:
        env_file (Path): The path to the environment YAML file.
        env_root_dir (Path, optional): The root directory where environments are stored.
                                       Defaults to the CWD's DEFAULT_ENV_ROOT_DIR.
        saved_env_yaml_dir (Path, optional): The directory where environment YAML files are saved 
                                             for future comparison. Defaults to SAVED_ENV_YAML_DIR.

    Returns:
        Path: The full path to the created or activated conda environment, or None if creation failed.
    """

    env_root_dir = env_root_dir or Path.cwd() / '.tool' / 'conda'
    saved_env_yaml_dir = saved_env_yaml_dir or Path.cwd() / '.tool' / 'saved_envs'

    env_root_dir.mkdir(parents=True, exist_ok=True)
    get_archive_dir().mkdir(parents=True, exist_ok=True)
    saved_env_yaml_dir.mkdir(parents=True, exist_ok=True)

    env_name = env_file.stem
    env_dir = env_root_dir / env_name
    saved_env_yaml_file = saved_env_yaml_dir / f"{env_name}.yaml"

    if env_dir.exists():
        if saved_env_yaml_file.exists() and filecmp.cmp(env_file, saved_env_yaml_file, shallow=False):
            file_log.environment(f"Environment {env_name} is up to date. Activating it.")
            activate_env(env_dir, env_name)
            return env_dir
        elif not saved_env_yaml_file.exists():
            file_log.environment(f"Environment {env_name} exists but no saved YAML found. Assuming up-to-date.")
            activate_env(env_dir, env_name)
            return env_dir

        file_log.environment(f"Environment {env_name} has changed. Recreating it.")
        archive_path = get_archive_dir() / env_name
        if archive_path.exists():
            shutil.rmtree(archive_path)
        shutil.move(env_dir, archive_path)

    console_log.environment(f"Creating environment {env_name}. Install mamba for faster environment creation.")
    if shutil.which("mamba") is not None:
        cmd = ["mamba", "env", "create", "-f", str(env_file), "--prefix", str(env_dir)]
    else:
        cmd = ["conda", "env", "create", "-f", str(env_file), "--prefix", str(env_dir)]

    result = subprocess.run(
        args=cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        file_log.error(f"Error creating environment: {result.stderr}")
        return None

    shutil.copy(env_file, saved_env_yaml_file)

    file_log.environment(f"Activating environment {env_name}.")
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
    file_log.environment(f"Activated environment {env_name} at {env_path}")


def deactivate_env():
    """
    Deactivates the current conda environment by removing its paths from the PATH environment variable and unsetting CONDA_PREFIX.
    """
    conda_prefix = os.environ.pop("CONDA_PREFIX", None)
    if conda_prefix:
        env_name = Path(conda_prefix).name
        os.environ["PATH"] = ":".join(p for p in os.environ["PATH"].split(":") if not p.startswith(conda_prefix))
        file_log.environment(f"Deactivated environment {env_name} at {conda_prefix}")
