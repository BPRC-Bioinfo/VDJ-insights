import datetime
import json
import sys
import shutil
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
import yaml
from Bio import SeqIO
import pandas as pd

# Global CONFIG
CONFIG = None
IMAGE_DIR = Path('source/images')


def load_config(cwd: Path):
    """
    Load the configuration file.

    Args:
        cwd (Path): Current working directory.
    """
    global CONFIG
    config_file = cwd / 'config' / 'config.yaml'
    if config_file.exists():
        with open(config_file, 'r') as file:
            CONFIG = yaml.safe_load(file)
    else:
        sys.exit("No configuration file provided, closing application!")


def make_dir(directory: Path):
    """
    Create a directory if it does not exist.

    Args:
        directory (Path): Path of the directory to create.
    """
    directory.mkdir(parents=True, exist_ok=True)


def move_dir(src: Path, dst: Path):
    """
    Copy a directory to a new location if it does not already exist.

    Args:
        src (Path): Source directory.
        dst (Path): Destination directory.
    """
    if not dst.exists():
        shutil.copytree(src, dst)


def parse_qc_files(qc_dir: Path) -> dict:
    """
    Parse quality control (QC) files.

    Args:
        qc_dir (Path): Directory containing QC files.

    Returns:
        dict: Parsed QC data.
    """
    data = {}
    for subfolder in qc_dir.iterdir():
        if subfolder.is_dir():
            subfolder_data = []
            for file in subfolder.glob("*.stats"):
                with open(file, 'r') as f:
                    lines = f.readlines()
                    header = lines[0].strip().split()
                    values = lines[1].strip().split()
                    file_data = dict(zip(header, values))
                    file_data['file'] = str(file.relative_to(qc_dir))
                    subfolder_data.append(file_data)
            data[subfolder.name] = subfolder_data
    return data


def parse_quast_data(quast_dir: Path) -> dict:
    """
    Parse QUAST data.

    Args:
        quast_dir (Path): Directory containing QUAST data.

    Returns:
        dict: Parsed QUAST data.
    """
    data = {}
    for subfolder in sorted(quast_dir.iterdir(), reverse=True):
        if subfolder.is_dir():
            subfolder_data = {'reports': [],
                              'pdfs': [], 'transposed_report': None}
            for report in subfolder.glob("**/report.*"):
                subfolder_data['reports'].append(
                    str(report.relative_to(quast_dir.parent)))
            basic_stats_dir = subfolder / 'basic_stats'
            moved_location = IMAGE_DIR / subfolder.name
            if basic_stats_dir.exists():
                move_dir(basic_stats_dir, moved_location)
                for pdf in moved_location.glob("*.pdf"):
                    subfolder_data['pdfs'].append(
                        str(Path('..', *pdf.parts[1:])))
            transposed_report_file = subfolder / 'transposed_report.tsv'
            if transposed_report_file.exists():
                subfolder_data['transposed_report'] = parse_transposed_report(
                    transposed_report_file)
            data[subfolder.name.split("_")[0]] = subfolder_data
    return data


def parse_transposed_report(report_path: Path) -> list:
    """
    Parse a transposed report file.

    Args:
        report_path (Path): Path to the transposed report file.

    Returns:
        list: Parsed report data.
    """
    df = pd.read_csv(report_path, sep='\t')
    return df.to_dict(orient='records')


def parse_region_files(region_dir: Path) -> list:
    """
    Parse region files.

    Args:
        region_dir (Path): Directory containing region files.

    Returns:
        list: Parsed region data.
    """
    data = []
    for fasta_file in region_dir.glob("*.fasta"):
        file_info = {}
        parts = fasta_file.stem.split('_')
        if len(parts) == 4:
            file_info['sample_name'] = parts[0]
            file_info['region_name'] = "/".join(parts[1:3])
            file_info['haplotype'] = parts[-1]
        length = sum(len(record.seq)
                     for record in SeqIO.parse(fasta_file, "fasta"))
        file_info['length'] = length
        file_info['file'] = str(fasta_file.relative_to(region_dir))
        data.append(file_info)
    return data


def parse_rss_meme_data(rss_dir: Path) -> dict:
    """
    Parse RSS meme data.

    Args:
        rss_dir (Path): Directory containing RSS meme data.

    Returns:
        dict: Parsed RSS meme data.
    """
    data = {}
    for main_folder in rss_dir.glob('*_meme'):
        if main_folder.is_dir():
            data[main_folder.name] = {}
            for subfolder in sorted(main_folder.iterdir()):
                if subfolder.is_dir():
                    type_key = subfolder.name[:3]
                    if type_key not in data[main_folder.name]:
                        data[main_folder.name][type_key] = []
                    logo_file = subfolder / 'logo1.png'
                    if logo_file.exists():
                        destination = IMAGE_DIR / main_folder.name / subfolder.name
                        make_dir(destination)
                        shutil.copy(logo_file, destination / 'logo1.png')
                        data[main_folder.name][type_key].append({
                            'folder': subfolder.name,
                            'logo': str(Path('..', *destination.parts[1:]) / 'logo1.png')
                        })
    return data


def summarize_annotation_data(file_path: Path) -> tuple:
    """
    Summarize annotation data from an Excel file.

    Args:
        file_path (Path): Path to the Excel file.

    Returns:
        tuple: Summary data and detailed annotation data.
    """
    df = pd.read_excel(file_path)
    summary = df.groupby(['Region', 'Haplotype', 'Function',
                         'Segment']).size().reset_index(name='Count')
    return summary.to_dict(orient='records'), df.to_dict(orient='records')


def parse_config_to_html(cwd: Path, data: dict):
    """
    Parse configuration and render HTML report.

    Args:
        cwd (Path): Current working directory.
        data (dict): Data to be passed to the template for rendering.
    """
    env = Environment(loader=FileSystemLoader(str(cwd / 'source' / 'template')))
    template = env.get_template('report_template.html')

    # Define the current date and time
    data['date_time'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Render the template with the provided data
    html_content = template.render(data)

    output_dir = cwd / 'source' / 'html'
    make_dir(output_dir)
    with open(output_dir / 'report.html', 'w') as html_file:
        html_file.write(html_content)

    # Render the full annotation details in separate HTML files
    details_template = env.get_template('annotation_details_template.html')

    known_details_html_content = details_template.render(
        date_time=data['date_time'],
        annotation_data_100_plus=data['annotation_data_100_plus'],
        annotation_data_plus=[],
        show_known=True
    )

    with open(output_dir / 'annotation_details_known.html', 'w') as html_file:
        html_file.write(known_details_html_content)

    novel_details_html_content = details_template.render(
        date_time=data['date_time'],
        annotation_data_100_plus=[],
        annotation_data_plus=data['annotation_data_plus'],
        show_known=False
    )

    with open(output_dir / 'annotation_details_novel.html', 'w') as html_file:
        html_file.write(novel_details_html_content)


def load_BUSCO_files(busco_dir: Path) -> dict:
    """
    Load BUSCO files and extract relevant metrics.

    Args:
        busco_dir (Path): Directory containing BUSCO files.

    Returns:
        dict: Extracted BUSCO metrics for each chromosome and haplotype.
    """
    busco_files = {}
    for subdirc in busco_dir.glob('*'):
        for chromosome in subdirc.glob('*'):
            file_structure = chromosome.stem.split("_")
            number = file_structure[0]
            haplotype = file_structure[2] if len(
                file_structure) > 2 and file_structure[2] else "reference"
            busco_file = next(chromosome.glob("*.json"))
            with open(busco_file, 'r') as js:
                content = json.load(js)
                results = content.get('results', {})
                metrics = {
                    'Single BUSCOs': results.get('Single copy BUSCOs', 'N/A'),
                    'Fragmented BUSCOs': results.get('Fragmented BUSCOs', 'N/A'),
                    'Missing BUSCOs': results.get('Missing BUSCOs', 'N/A'),
                    'Duplicates BUSCOs': results.get('Multi copy BUSCOs', 'N/A')
                }
                busco_files.setdefault(
                    number, {}).setdefault(haplotype, metrics)
    return busco_files


def html_main():
    cwd = Path.cwd()
    make_dir(cwd / IMAGE_DIR)
    load_config(cwd)

    busco_data = load_BUSCO_files(cwd / 'BUSCO')
    qc_data = parse_qc_files(cwd / 'QC')
    quast_data = parse_quast_data(cwd / 'quast' / 'hifiasm')
    region_data = parse_region_files(cwd / 'region')
    rss_meme_data = parse_rss_meme_data(cwd / 'RSS')
    annotation_summary_100_plus, annotation_data_100_plus = summarize_annotation_data(
        cwd / 'annotation' / 'annotation_report_100%_plus.xlsx')
    annotation_summary_plus, annotation_data_plus = summarize_annotation_data(
        cwd / 'annotation' / 'annotation_report_plus.xlsx')

    data = {
        'config': CONFIG,
        'busco_data': busco_data,
        'qc_data': qc_data,
        'quast_data': quast_data,
        'region_data': region_data,
        'rss_meme_data': rss_meme_data,
        'annotation_summary_100_plus': annotation_summary_100_plus,
        'annotation_summary_plus': annotation_summary_plus,
        'annotation_data_100_plus': annotation_data_100_plus,
        'annotation_data_plus': annotation_data_plus
    }

    parse_config_to_html(cwd, data)


if __name__ == '__main__':
    html_main()
