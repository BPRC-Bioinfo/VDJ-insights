import subprocess
import tempfile
from bowtie2 import bowtie2_main
from minimap2 import minimap2_main
from write_annotation_report import write_report_report
from pathlib import Path
import pandas as pd


def combine_df():
    minimap2_df = minimap2_main()
    bowtie2_df = bowtie2_main()
    return pd.concat([minimap2_df, bowtie2_df])


def write_report(df, report):
    df.to_excel(report, index=False)


def select_row(df):
    if 'V' in df['segment'].values:
        return df[df['tool'] == 'minimap2']
    else:
        return df[df['tool'] == 'bowtie2']


def get_or_create(annotation_folder):
    report = annotation_folder / "report.xlsx"
    if not report.exists():
        df = combine_df()
        write_report(df, report)
        return df
    else:
        return pd.read_excel(report)


def make_dir(dir):
    Path(dir).mkdir(parents=True, exist_ok=True)


def make_blast_db(cwd):
    db = cwd / "blast_db"
    if not db.exists():
        reference = cwd / "library" / "retained.fasta"
        make_dir(db)
        command = f"makeblastdb -in {reference} -dbtype nucl -out {db}/blast_db"
        subprocess.run(command, shell=True)
    return db


def run_blast(row, db, cut_off) -> str:
    header, sequence, start, stop, path, strand = row['name'], row[
        "sequence"], row["start"], row["stop"], row["fasta-file"], row["strand"]
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as fasta_file, \
            tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as result_file:
        fasta_file.write(
            f">{header}:{start}:{stop}:{strand}:{path}\n{sequence}\n")
        fasta_file.flush()
        blast_columns = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
        command = f"blastn -query {fasta_file.name} -db {db}/blast_db -outfmt '{blast_columns}' -perc_identity {cut_off} -out {result_file.name}"
        subprocess.run(command, shell=True)
        result_file.flush()
        with open(result_file.name, 'r') as file:
            blast_result = file.readlines()
            if blast_result:
                all_hits = []
                for number, hit in enumerate(blast_result):
                    hits = hit.strip().split("\t")
                    hits.extend([number + 1, cut_off, start, stop])
                    all_hits.append(hits)
                return all_hits


def make_blast_df(filtered_df: pd.DataFrame, db):
    df = pd.DataFrame()
    for cut_off in [100, 75, 50]:
        blast_results = filtered_df.apply(
            lambda row: run_blast(row, db, cut_off), axis=1)
        blast_results = blast_results.dropna().reset_index(drop=True)
        blast_columns = [
            'query', 'subject', '% identity',
            'alignment length', 'mismatches', 'gap opens',
            'q. start', 'q. end', 's. start',
            's. end', 'evalue', 'bit score',
            'query seq', 'subject seq', 'number',
            'cut-off', 'start', 'stop'
        ]
        flattened_hits = [hit for sublist in blast_results for hit in sublist]
        blast_results_df = pd.DataFrame(flattened_hits, columns=blast_columns)
        df = pd.concat([df, blast_results_df], ignore_index=True)
    return df


def main():
    cwd = Path.cwd()
    annotation_folder = cwd / "annotation"
    make_dir(annotation_folder)
    db = make_blast_db(cwd)
    df = get_or_create(annotation_folder)
    filtered_df = df.groupby(["name", "start", "stop"]).apply(select_row)
    filtered_df: pd.DataFrame = filtered_df.reset_index(drop=True)
    filtered_df = filtered_df.query("region != 'LOC'")
    blast_df = make_blast_df(filtered_df, db)
    blast_df = blast_df.drop_duplicates(
        subset=blast_df.columns[:12]).reset_index(drop=True)
    blast_df.to_excel(annotation_folder / "blast_results.xlsx", index=False)
    write_report_report(annotation_folder)


if __name__ == '__main__':
    main()
