from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import subprocess
import shutil

from .util import make_dir, calculate_available_resources


def run_scaffolding(assembly_file: str, reference: str, scaffolding_dir: str):
    cwd = Path.cwd()

    assembly_name = assembly_file.stem
    #if not assembly_name in reference_name:
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
    return Path(scaffolding_dir)

