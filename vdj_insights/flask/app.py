import re
from flask import Flask, render_template, request, jsonify
from pathlib import Path
import datetime
import yaml
import shutil
import pandas as pd
from Bio import SeqIO
import json
import sys

# Initialize the Flask application
app = Flask(__name__)

# Global CONFIG
CONFIG = None
IMAGE_DIR = Path('flask/static/images')
MODE = None


def load_config(cwd: Path):
    global CONFIG
    config_file = cwd / 'config' / 'config.yaml'
    if config_file.exists():
        with open(config_file, 'r') as file:
            CONFIG = yaml.safe_load(file)
    else:
        sys.exit("No configuration file provided, closing application!")


def make_dir(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)


def move_dir(src: Path, dst: Path):
    if not dst.exists():
        shutil.copytree(src, dst)


def parse_qc_files(qc_dir: Path) -> dict:
    data = {}
    if qc_dir.is_dir():
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
    data = {}
    if quast_dir.is_dir():
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


def extract_sample(path):
    filename = path.split("/")[-1]
    sample_pattern = re.compile(r'(GCA|GCF|DRR|ERR)_?\d{6,9}(\.\d+)?')
    match = sample_pattern.search(filename)
    if match:
        return match.group(0)
    else:
        return filename.split("_")[0]


def parse_transposed_report(report_path: Path) -> list:
    df = pd.read_csv(report_path, sep='\t')
    return df.to_dict(orient='records')


def parse_region_files(region_dir: Path) -> list:
    data = []
    if region_dir.is_dir():
        for fasta_file in region_dir.glob("*.fasta"):

            file_info = {}
            sample = extract_sample(str(fasta_file))
            file_info['sample_name'] = sample
            parts = fasta_file.stem
            parts = parts.removeprefix(f"{sample}_").split('_')
            file_info['region_name'] = "/".join(parts[0:2])
            file_info['haplotype'] = parts[-1]
            length = sum(len(record.seq)
                         for record in SeqIO.parse(fasta_file, "fasta"))
            file_info['length'] = length
            file_info['file'] = str(fasta_file.relative_to(region_dir))
            data.append(file_info)
    return data


def parse_rss_meme_data(rss_dir: Path) -> dict:
    data = {}
    if rss_dir.is_dir():
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
                            if not destination.is_file():
                                make_dir(destination)
                                shutil.copy(
                                    logo_file, destination / 'logo1.png')
                            data[main_folder.name][type_key].append({
                                'folder': subfolder.name,
                                'logo': str(Path('..', *destination.parts[1:]) / 'logo1.png')
                            })
    return data


def summarize_annotation_data(file_path: Path) -> tuple:
    if file_path.is_file():
        df = pd.read_excel(file_path)
        summary = df.groupby(
            ['Region', 'Haplotype', 'Function', 'Segment']).size().reset_index(name='Count')
        annotation = df.to_dict(orient='records')
        annotation_summary = summary.to_dict(orient='records')
    else:
        annotation, annotation_summary = [], []
    return annotation_summary, annotation


def load_BUSCO_files(busco_dir: Path) -> dict:
    data = {}
    if busco_dir.is_dir():
        for subdir in busco_dir.glob('*'):
            for chromosome in subdir.glob('*'):
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
                    data.setdefault(number, {}).setdefault(haplotype, metrics)
    return data


@app.errorhandler(404)
def page_not_found(e):
    # You can return a custom template or plain text response
    return render_template('404.html'), 404

# Route for the home page


@app.route('/')
def home():
    context = {'date_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
    return render_template('index.html', **context)


# Route for the report page
@app.route('/report')
def report():
    cwd = Path.cwd()
    load_config(cwd)
    make_dir(cwd / IMAGE_DIR)

    # Loading different types of data
    busco_data = load_BUSCO_files(cwd / 'BUSCO')
    qc_data = parse_qc_files(cwd / 'QC')
    quast_data = parse_quast_data(cwd / 'quast' / 'hifiasm')
    region_data = parse_region_files(cwd / 'region')
    rss_meme_data = parse_rss_meme_data(cwd / 'RSS')

    # Parsing annotation data
    annotation_summary_100_plus, annotation_data_100_plus = summarize_annotation_data(
        cwd / 'annotation' / 'annotation_report_known_rss.xlsx')
    annotation_summary_plus, annotation_data_plus = summarize_annotation_data(
        cwd / 'annotation' / 'annotation_report_novel_rss.xlsx')

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
        'annotation_data_plus': annotation_data_plus,
        'date_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'mode': 'Pipeline' if qc_data else 'Annotation'
    }

    return render_template('report.html', **data)

def filter_dict(broken: dict, partial_key: str) -> dict:
    return {key: value for key, value in broken.items() if partial_key in key}


@app.route('/<string:sample_id>')
def sample_details(sample_id):
    cwd = Path.cwd()
    individual_path = cwd / 'annotation' / 'individual' / sample_id

    if not individual_path.exists():
        return render_template('404.html'), 404

    region_data = parse_region_files(cwd / 'region')
    
    sample_region_data = [item for item in region_data if item['sample_name'] == sample_id]
    with open(cwd / "broken_regions.json", "r") as f:
        broken = list(filter_dict(json.load(f), sample_id).values())[0]
    region_data_dict = {
        element["region_name"].replace("/", "-"): {
            "haplotype": element["haplotype"],
            "length": element["length"],
            "file": element["file"],
            "broken": broken.get(element["region_name"].replace("/", "-")),
        }
        for element in sample_region_data
    }


    known_report_file = individual_path / 'annotation_report_known_rss.xlsx'
    novel_report_file = individual_path / 'annotation_report_novel_rss.xlsx'

    known_annotation_summary = {}
    novel_annotation_summary = {}

    if known_report_file.exists():
        known_df = pd.read_excel(known_report_file)

        known_annotation_summary = {
            'total_annotations': len(known_df),
            'unique_segments': known_df['Segment'].nunique(),
            'segments': known_df['Segment'].value_counts().to_dict(),
            'functions': known_df['Function'].value_counts().to_dict(),
            'regions': known_df['Region'].value_counts().to_dict(),
            'haplotypes': known_df['Haplotype'].value_counts().to_dict(),
            'average_mismatches': known_df['Mismatches'].mean(),
            'max_mismatches': known_df['Mismatches'].max(),
        }
    else:
        known_annotation_summary = None

    if novel_report_file.exists():
        novel_df = pd.read_excel(novel_report_file)

        # Compute statistics
        novel_annotation_summary = {
            'total_annotations': len(novel_df),
            'unique_segments': novel_df['Segment'].nunique(),
            'segments': novel_df['Segment'].value_counts().to_dict(),
            'functions': novel_df['Function'].value_counts().to_dict(),
            'regions': novel_df['Region'].value_counts().to_dict(),
            'haplotypes': novel_df['Haplotype'].value_counts().to_dict(),
            'average_mismatches': novel_df['Mismatches'].mean(),
            'max_mismatches': novel_df['Mismatches'].max(),
        }
    else:
        novel_annotation_summary = None

    # Sample metadata (replace with actual data retrieval)
    sample_metadata = {
        'source': 'Blood',
        'collection_date': '2023-10-01',
        'description': 'Sample collected from patient XYZ.'
    }

    context = {
        'sample_id': sample_id,
        'regions': region_data_dict,
        'known_annotation_summary': known_annotation_summary,
        'novel_annotation_summary': novel_annotation_summary,
        'date_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'sample_metadata': sample_metadata
    }

    return render_template('sample_details.html', **context)


@app.route('/annotation/<string:annotation_type>', defaults={'sample_id': 'all'})
@app.route('/annotation/<string:annotation_type>/<string:sample_id>')
def annotation(annotation_type, sample_id):
    cwd = Path.cwd()
    individual_path = cwd / 'annotation' / 'individual'
    if sample_id == "all" or sample_id == None:
        report_path = cwd / 'annotation'
    elif sample_id:
        report_path = individual_path / sample_id

    if annotation_type == 'known':
        report_file = report_path / 'annotation_report_known_rss.xlsx'
        show_known = True
    else:
        report_file = report_path / 'annotation_report_novel_rss.xlsx'
        show_known = False

    annotation_summary, annotation_data = summarize_annotation_data(report_file)
    all_samples = [i.name for i in individual_path.glob(
        "*") if not i.name.startswith(".")]
    data = {
        'annotation_summary': annotation_summary,
        'annotation_data_100_plus': annotation_data if show_known else None,
        'annotation_data_plus': None if show_known else annotation_data,
        'show_known': show_known,
        'date_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'sample_ids': all_samples,
        'current_sample': sample_id
    }

    return render_template('annotation.html', **data)


def read_fasta(file_path):
    sequences = {}
    header = None
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:
            sequences[header] = ''.join(sequence)  # Add the last sequence

    return sequences

# Function to append new sequences to FASTA file


def append_to_fasta(file_path, new_sequences):
    existing_sequences = read_fasta(file_path)
    total_new = 0
    with open(file_path, 'a') as file:
        for header, sequence in new_sequences.items():
            if header not in existing_sequences and sequence not in existing_sequences.values():
                total_new += 1
                file.write(f"{header}\n{sequence}\n")
    return total_new


@app.route('/save_sequences', methods=['POST'])
def save_sequences():
    new_sequences = request.json.get('sequences', {})
    base_dir = Path.cwd()
    fasta_file_path = base_dir / ".tool/library/library.fasta"

    if not new_sequences:
        return jsonify({"error": "No sequences provided"}), 400

    total_new = append_to_fasta(fasta_file_path, new_sequences)
    if total_new > 0:
        message = f"Added a total of {total_new} new sequences to the library"
    else:
        message = "No new sequences found to add"
    return jsonify({"message": message}), 200


@app.route('/remove_sequences', methods=['POST'])
def remove_sequences():
    sequences_to_remove = request.json.get('sequences', [])
    base_dir = Path.cwd()
    fasta_file_path = base_dir / ".tool/library/library.fasta"

    if not sequences_to_remove:
        return jsonify({"error": "No sequences provided"}), 400

    # Read existing sequences
    existing_sequences = read_fasta(fasta_file_path)

    # Filter out the sequences to remove
    remaining_sequences = {header: seq for header, seq in existing_sequences.items(
    ) if header not in sequences_to_remove}

    # Overwrite the FASTA file with the remaining sequences
    with open(fasta_file_path, 'w') as fasta_file:
        for header, sequence in remaining_sequences.items():
            fasta_file.write(f"{header}\n{sequence}\n")

    return jsonify({"message": f"Removed {len(sequences_to_remove)} sequences from the library"}), 200


@app.route('/download_fasta', methods=['GET'])
def download_fasta():
    base_dir = Path.cwd()
    fasta_file_path = base_dir / ".tool/library/library.fasta"

    if not fasta_file_path.exists():
        return jsonify({"error": "FASTA file not found"}), 404

    with open(fasta_file_path, 'r') as file:
        fasta_content = file.read()

    # Set the response headers to trigger a download
    response = app.response_class(
        fasta_content,
        mimetype='text/plain',
        headers={"Content-Disposition": "attachment;filename=library.fasta"}
    )

    return response


@app.route('/imgt_report')
def imgt_report():
    base_dir = Path.cwd()
    context = json.load(open(base_dir / "library" / "library_info.json"))
    context["date_time"] = datetime.datetime.now().strftime(
        "%Y-%m-%d %H:%M:%S")
    return render_template('imgt_report.html', **context)


@app.route('/library')
def view_library():
    base_dir = Path.cwd()
    fasta_file_path = base_dir / ".tool/library/library.fasta"
    sequences = read_fasta(fasta_file_path)
    date_time = datetime.datetime.now().strftime(
        "%Y-%m-%d %H:%M:%S")  # Get current date and time
    context = {
        'sequences': sequences,
        'date_time': date_time
    }
    return render_template('library.html', **context)


@app.route('/sequence_library')
def sequence_library():
    base_dir = Path.cwd()
    sequence_json = base_dir / ".tool/library/library.json"
    with open(sequence_json) as f:
        data = json.load(f)
    return render_template('sequence_library.html', data=data)


@app.route('/about')
def about():
    return render_template('about.html')


# Run the application
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
