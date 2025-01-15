"""
Copyright (c) 2023-2025 Biomedical Primate Research Centre, the Netherlands.
All rights reserved.
"""

from pathlib import Path
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO
from order_segments import order_main

from .util import make_dir
from .logger import console_logger, file_logger

cwd = Path.cwd()

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def re_evaluate_main(excel_files):
    library_path = cwd / "library.fasta"



if __name__ == "__main__":
    re_evaluate_main()
