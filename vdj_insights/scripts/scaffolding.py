from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import subprocess
import shutil
import json
import pandas as pd

from .util import make_dir, calculate_available_resources
from .logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)


def make_json_scaffolds():
    cwd = Path.cwd()
    ragtag_path = cwd / Path("tmp/ragtag")

    scaffold_agp_files = [Path(agp_file) for agp_file in ragtag_path.rglob('*') if agp_file.suffix.lower() == '.agp']
    all_scaffolds = {}
    for agp_file in scaffold_agp_files:
        try:
            assembly = agp_file.parent.name
            agp_df = pd.read_csv(agp_file, sep="\t", header=None, skiprows=2)
            agp_df = agp_df[agp_df[4] == "W"]
            scaffold_contigs = agp_df.groupby(0)[5].apply(list)
            scaffold_dict = scaffold_contigs.to_dict()
            all_scaffolds[assembly] = scaffold_dict
        except Exception as e:
            file_log.error(f"Error processing {agp_file}: {e}")

    json_path = ragtag_path / "ragtag_scaffolds.json"
    with open(json_path, 'w') as f:
        json.dump(all_scaffolds, f, indent=2)


def run_scaffolding(assembly_file: str, reference: str, scaffolding_dir: str):
    cwd = Path.cwd()

    assembly_name = assembly_file.stem
    ragtag_output = cwd / Path("tmp/ragtag") / assembly_name
    make_dir(ragtag_output)

    original_scaffold = cwd / ragtag_output / "ragtag.scaffold.fasta"
    scaffolding_dir = cwd /scaffolding_dir
    make_dir(scaffolding_dir)
    renamed_scaffold = scaffolding_dir / f"{assembly_name}.fasta"
    if not renamed_scaffold.is_file():
        if ragtag_output.exists():
            shutil.rmtree(ragtag_output)
        rag_command = f"ragtag.py scaffold -t 4 {reference} {assembly_file} -o {ragtag_output}/"
        subprocess.run(rag_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        shutil.move(str(original_scaffold), str(renamed_scaffold))


def scaffolding_main(reference: str, assembly_dir: Path, scaffolding_dir: str,  threads: int) -> Path:
    cwd = Path.cwd()
    assembly_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in (cwd / assembly_dir).glob(ext)]
    scaffold_assembly_files = [file for ext in ["*.fna", "*.fasta", "*.fa"] for file in (cwd / scaffolding_dir).glob(ext)]
    if not len(assembly_files) == len(scaffold_assembly_files):
        total_tasks = len(assembly_files)
        max_jobs = calculate_available_resources(max_cores=threads, threads=4, memory_per_process=12)
        with ThreadPoolExecutor(max_workers=max_jobs) as executor:
            futures = [executor.submit(run_scaffolding, assembly_file, reference, scaffolding_dir)for assembly_file in assembly_files]
            with tqdm(total=total_tasks, desc='Scaffolding:', unit="Assemblies") as pbar:
                for future in as_completed(futures):
                    pbar.update(1)

    make_json_scaffolds()
    return Path(scaffolding_dir)

