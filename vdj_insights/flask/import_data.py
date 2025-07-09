import json
import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def open_json(file_name: str) -> dict:
    if os.path.exists(file_name):
        with open(file_name) as file:
            data = json.load(file)
        return data


def get_sequences(file_name: str) -> list:
    sequences = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append({
            "id": record.id,
            "sequence_length": len(record.seq),
            "sequence": str(record.seq)
        })
    return sequences


def get_region_data(path: Path):
    try:
        region_json = open_json(path / "broken_regions.json")
        records = []
        if region_json:
            for file, segments in region_json.items():
                for segment, details in segments.items():
                    for region in details:
                        record = {"File": file, "Segment": segment, **region}
                        records.append(record)
        return pd.DataFrame(records)
    except FileNotFoundError:
        return pd.DataFrame()


def get_commando_data(path: Path):
    try:
        commando_json = open_json(path / "used_commando.json")
        records = []
        if commando_json:
            for command, values in commando_json.items():
                record = {"Argument": command, "Given argument": values}
                records.append(record)
        return pd.DataFrame(records)
    except FileNotFoundError:
        return pd.DataFrame()


def get_scaffold_data(path: Path):
    try:
        scaffold_json = open_json(path / "ragtag_scaffolds.json")
        records = []
        if scaffold_json:
            for file, scaffold_dict in scaffold_json.items():
                for scaffold_name, contigs in scaffold_dict.items():
                    for contig in contigs:
                        records.append({
                            "File": file,
                            "Scaffold": scaffold_name,
                            "Contig": contig
                        })
        return pd.DataFrame(records)
    except FileNotFoundError:
        return pd.DataFrame()