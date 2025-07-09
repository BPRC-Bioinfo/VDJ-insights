import datetime
import os
import json
import shutil
import base64
from pathlib import Path
from datetime import datetime
import io

import jsonify
import pandas as pd
import numpy as np
from flask import Flask, render_template, request, make_response, flash, redirect, url_for, send_file, Response, session, abort

from celery import Celery
from flask_caching import Cache
import base64
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import uuid
import subprocess
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.express as px
from plotly.offline import plot

from scipy.cluster.hierarchy import linkage
import plotly.figure_factory as ff

from task import merge_sequences

from import_data import open_json, get_sequences, get_region_data, get_scaffold_data, get_commando_data

from figures.figures import get_known_novel_plot
from figures.figures import get_vdj_plot
from figures.figures import get_shared_shortname_heatmap_plotly_simple
from figures.figures import get_sample_rss_plot
from figures.figures import get_rss_plot
from figures.figures import get_library_violin_plot
from figures.figures import get_library_plot
from figures.figures import plot_count_distribution_bar_chart
from figures.figures import get_dna_vieuwer_plot
from figures.figures import get_pi_plot
from figures.figures import get_venn_diagram
from figures.figures import generate_bokeh_segment_div_from_df

from HaplotypeAligner import HaplotypeAligner
from NetworkBuilder import  NetworkBuilder


app = Flask(__name__)
app.secret_key = os.urandom(12)

app.config["CACHE_TYPE"] = "SimpleCache"
app.config["CACHE_DEFAULT_TIMEOUT"] = 3600
app.config["CELERY_BROKER_URL"] = "redis://localhost:6379/0"
app.config["CELERY_RESULT_BACKEND"] = "redis://localhost:6379/0"

celery = Celery(app.name, broker=app.config["CELERY_BROKER_URL"])
celery.conf.update(app.config)

cache = Cache(app)

BASE_PATH = Path("")
LIBRARY_PATH =  BASE_PATH / "library"
RSS_PATH = BASE_PATH / "tmp/RSS"


@app.route('/', methods=['POST', 'GET'])
def home():
    commando_data = get_commando_data(BASE_PATH)
    input_row = commando_data.loc[commando_data["Argument"] == "input", "Given argument"].iloc[0]
    if input_row == None:
        commando_data = commando_data[commando_data["Argument"] != "input"]

    desired_order = [
        "species",
        "receptor_type",

        "assembly",
        "input",

        "flanking_genes",
        "library",
        "mapping_tool",

        "scaffolding",
        "metadata",

        "output",

        "threads",
        "verbose",
        "command line",
    ]
    present_args = [arg for arg in desired_order if arg in commando_data["Argument"].values]
    df_sorted = commando_data.set_index("Argument").loc[present_args].reset_index()
    df_sorted["Argument"] = df_sorted["Argument"].str.title()
    return render_template("home.html", commando_data=df_sorted.to_dict(orient='records'))


@app.route('/help_annotation', methods=['POST', 'GET'])
def help_annotation():
    return render_template("help_annotation.html")

@app.route('/get_compare', methods=['POST', 'GET'])
def get_compare():
    df = get_annotation_data()

    samples = df["Sample"].unique()
    regions = df["Region"].unique()

    selected_sample_one = session.get("selected_sample_one", None)
    selected_sample_two = session.get("selected_sample_two", None)
    selected_region = session.get("selected_region", None)

    data_table = pd.DataFrame()
    venn_svg_dict = {}

    if request.method == "POST":
        selected_sample_one = request.form.get("sample_one", None)
        selected_sample_two = request.form.get("sample_two", None)
        selected_region = request.form.get("region", None)

        session["selected_sample_one"] = selected_sample_one
        session["selected_sample_two"] = selected_sample_two
        session["selected_region"] = selected_region

        data_table = df[(df["Sample"].isin([selected_sample_one, selected_sample_two])) & (df["Region"] == selected_region)].copy()
        presence_list = []

        for status in data_table["Status"].unique():
            short_names_sample_one = set(data_table[(data_table["Sample"] == selected_sample_one) & (data_table["Status"] == status)]["Short name"].unique())
            short_names_sample_two = set(data_table[(data_table["Sample"] == selected_sample_two) & (data_table["Status"] == status)]["Short name"].unique())

            if not short_names_sample_one or not short_names_sample_two:
                flash("One of the samples has no data1 in this region.", "warning")
                return redirect(url_for("get_compare"))

            venn_svg = get_venn_diagram(short_names_sample_one, short_names_sample_two, selected_sample_one, selected_sample_two, status)
            venn_svg_dict[f"venn_svg_dict_{status}"] = venn_svg

            for short_name in short_names_sample_one | short_names_sample_two:
                if short_name in short_names_sample_one and short_name in short_names_sample_two:
                    presence = "Shared"
                elif short_name in short_names_sample_one:
                    presence = f"Only in {selected_sample_one}"
                elif short_name in short_names_sample_two:
                    presence = f"Only in {selected_sample_two}"
                else:
                    presence = "Unknown"

                presence_list.append({"Short name": short_name, "Status": status, "Presence": presence})
        if len(presence_list) > 0:
            presence_df = pd.DataFrame(presence_list)
            data_table = data_table.merge(presence_df, on=["Short name", "Status"], how="left")

        return render_template("sample_compare.html",
                                   samples=samples,
                                   regions=regions,
                                   selected_sample_one=selected_sample_one,
                                   selected_sample_two=selected_sample_two,
                                   selected_region=selected_region,
                                   venn_svg=venn_svg_dict,
                                   data_table=data_table
                                   )

    return render_template("sample_compare.html",
                           samples=samples,
                           regions=regions,
                           data_table=data_table,
                           selected_sample_one=selected_sample_one,
                           selected_sample_two=selected_sample_two,
                           selected_region=selected_region,
                           )


def get_fasta(file_names):
    sequences = []
    for file_name in file_names:
        for record in SeqIO.parse(file_name, "fasta"):
            sequences.append({"id": record.id, "sequence": str(record.seq)})
    return sequences


@app.route('/run_mafft', methods=['POST', 'GET'])
def run_mafft():
    df = get_annotation_data()

    samples = df["Sample"].unique()
    regions = df["Region"].unique()

    selected_sample_one = None
    selected_sample_two = None
    selected_region = None

    if request.method == "POST":
        selected_sample_one = request.form.get("sample_one", None)
        selected_sample_two = request.form.get("sample_two", None)
        selected_region = request.form.get("region", None)

        if selected_sample_one == selected_sample_two:
            flash("Same sample given! Please choose two different samples.", "warning")
            return redirect(url_for("run_mafft"))
        if not selected_region:
            flash("Select region of interest", "warning")
            return redirect(url_for("run_mafft"))

        path_selected_sample_one = df[(df["Sample"] == selected_sample_one) & (df["Region"] == selected_region)]["Path"].unique()[0]
        path_selected_sample_two = df[(df["Sample"] == selected_sample_two) & (df["Region"] == selected_region)]["Path"].unique()[0]

        if not path_selected_sample_one or not path_selected_sample_two:
            flash("One of the samples does not have data for the selected region", "warning")
            return redirect(url_for("run_mafft"))

        UUID = uuid.uuid4().hex[:8]
        output_dir = os.path.join( BASE_PATH, f"tmp/mafft/mafft_{UUID}/")
        merged_dir = os.path.join(output_dir, "merged")
        aligned_dir = os.path.join(output_dir, "aligned")
        split_dir = os.path.join(output_dir, "split_regions")

        os.makedirs(merged_dir, exist_ok=True)
        os.makedirs(aligned_dir, exist_ok=True)
        os.makedirs(split_dir, exist_ok=True)

        merged_path = os.path.join(merged_dir, "merged_sequences.fasta")
        aligned_path = os.path.join(aligned_dir, "aligned_sequences.fasta")
        log_file = os.path.join(output_dir, "mafft.log")

        sequences = get_fasta([path_selected_sample_one, path_selected_sample_two])
        records = [SeqRecord(Seq(seq["sequence"]), id=seq["id"], description="") for seq in sequences]
        with open(merged_path, "w") as f:
            SeqIO.write(records, f, "fasta")


        conda_base = subprocess.run(["conda", "info", "--base"], stdout=subprocess.PIPE, text=True, check=True).stdout.strip()
        cmd = f'''nohup bash -c "source {conda_base}/etc/profile.d/conda.sh && conda activate .tool/conda/vdj-insights_env && mafft '{merged_path}' > '{aligned_path}'" > '{log_file}' 2>&1 &'''

        process = subprocess.Popen(cmd, shell=True)

        #if result.returncode != 0:
           # flash(f"{result.stderr}", "error")
           # return redirect(url_for("run_mafft"))
        """
        aligned_sequences = get_fasta([aligned_path])
        records = [SeqRecord(Seq(seq["sequence"]), id=seq["id"], description="") for seq in aligned_sequences]
        for record in records:
            fasta_file = os.path.join(split_dir, f"{record.id}.fasta")
            with open(fasta_file, "w") as f:
                SeqIO.write(record, f, "fasta")
        """

        return redirect(url_for("run_mafft"))

    return render_template("mafft.html",
                           samples=samples,
                           regions=regions,
                           selected_sample_one=selected_sample_one,
                           selected_sample_two=selected_sample_two,
                           selected_region=selected_region)


def load_tooltip_data(data: pd.DataFrame, accessions: list, region: str, set_segments: list) -> dict:
    tooltip_data = {}

    filtered_df = data[(data["Sample"].isin(accessions)) & (data["Short name"].isin(set_segments)) & (data["Region"] == region)].copy()

    filtered_df["rss_sum"] = filtered_df[["5'-RSS", "3'-RSS"]].sum(axis=1)
    conditions_rss = [
        (filtered_df["Segment"] == "D") & (filtered_df["rss_sum"] == 2),
        (filtered_df["Segment"] == "D") & (filtered_df["rss_sum"] == 1),
        (filtered_df["Segment"].isin(["V", "J"])) & (filtered_df["rss_sum"] == 1),
    ]
    rss_choices = ["", "", ""]
    filtered_df["found_rss_count"] = np.select(conditions_rss, rss_choices, default="/")
    filtered_df.drop(columns=["rss_sum"], inplace=True)

    for _, row in filtered_df.iterrows():
        gene = row["Short name"]
        sample = row["Sample"]
        tooltip_data[(sample, gene)] = {
            "full_name": row.get("Short name", ""),
            "function_segment": row.get("Function", ""),
            "region": row.get("Region", ""),
            "segment": row.get("Segment", ""),
            "strand": row.get("Strand", ""),
            "status": row.get("Status", ""),
            "SNPs": row.get("SNPs", ""),
            "Insertions": row.get("Insertions", ""),
            "Deletions": row.get("Deletions", ""),
            "found_rss_count": row.get("found_rss_count", ""),
            "5_RSS_seq": row.get("5'-RSS seq", ""),
            "3_RSS_seq": row.get("3'-RSS seq", ""),
        }

    return tooltip_data


def haplo_network_main(accessions: list, region: str, segment_types: list, data: pd.DataFrame, color_theme) -> dict:
    #filtered_df = data1[(data1["Sample"].isin(accessions)) & (data1["Region"] == region) & (data1["Segment"].isin(segment_types))].copy()
    filtered_df = data[(data["Sample"].isin(accessions)) & (data["Region"] == region)].copy()

    if not filtered_df.empty:
        hap1 = tuple(filtered_df[filtered_df["Sample"] == accessions[0]].sort_values(by="Start coord")[["Short name", "Start coord", "End coord"]].itertuples(index=False, name=None))
        hap2 = tuple(filtered_df[filtered_df["Sample"] == accessions[1]].sort_values(by="Start coord")[["Short name", "Start coord", "End coord"]].itertuples(index=False, name=None))

        set_segments = list({item[0] for item in hap1}.union({item[0] for item in hap2}))
        tooltip_data = load_tooltip_data(data=data, accessions=accessions, region=region, set_segments=set_segments)

        aligner = HaplotypeAligner(hap1, hap2)
        alignment = aligner.align()

        left_flanking_gene = "5'- flaking gene"
        right_flanking_gene = "3'- flaking gene"
        network_builder = NetworkBuilder(alignment, tooltip_data, left_flanking_gene, right_flanking_gene, accessions, color_theme)
        data_dict = network_builder.build_network()
        return data_dict


@app.route('/get_haplotype_alignment', methods=['POST', 'GET'])
def get_haplotype_alignment():
    accessions = []
    selected_accession_one = None
    selected_accession_two = None
    region = None
    color_theme = None

    accessions = ["GCA_018466835.1", "GCA_009914755.4"]

    selected_accession_one = "GCA_018466835.1"
    selected_accession_two = "GCA_009914755.4"
    region = "IGH"
    color_theme = "default"
    segment_types = ["V", "D", "J"]

    if request.method == "POST":
        region = request.form.get("region")
        selected_accession_one = request.form.get("accession_one")
        selected_accession_two = request.form.get("accession_two")
        accessions = [selected_accession_one, selected_accession_two]
        color_theme = request.form.get("color_theme", "default")
        segment_types = request.form.getlist("segment_type", None)

    data = get_annotation_data()
    data_dict = haplo_network_main(accessions, region, segment_types, data, color_theme)
    print(data_dict)
    return render_template(template_name_or_list="graph.html",
                           data=data_dict,
                           selected_accession_one=selected_accession_one,
                           selected_accession_two=selected_accession_two,
                           selected_region=region,
                           selected_color_theme=color_theme,
                           selected_segment_types=segment_types,
                           all_regions=data["Region"].unique(),
                           all_accession_groups=data["Sample"].unique())



def get_unique_sequences(reference):
    df = get_annotation_data()

    data = {}
    unique_shortname_rows = []

    for sample_naam in df['Sample'].unique():
        data[sample_naam] = {"Total_known": 0, "Total_novel": 0, "Total": 0}

        for status in ["Known", "Novel"]:
            df_status = df[(df["Sample"] == sample_naam) & (df["Status"] == status)]

            if not df_status.empty:
                short_names = df_status['Short name'].tolist()
                alle_short_names_lijst = df['Short name'].tolist()

                exclusieve_short_names = []
                for short_name in short_names:
                    if alle_short_names_lijst.count(short_name) == 1:
                        exclusieve_short_names.append(short_name)
                        if sample_naam == reference:
                            unique_shortname_rows.append(df_status[df_status["Short name"] == short_name])

                data[sample_naam][f"Total_{status.lower()}"] = len(exclusieve_short_names)
                data[sample_naam]["Total"] += len(exclusieve_short_names)

    df = pd.DataFrame.from_dict(data, orient="index").reset_index()

    data_table = df[["index", "Total_known", "Total_novel", "Total"]]
    data_table.rename(columns={"index": "Sample", "Total_known": "Known", "Total_novel": "Novel"}, inplace=True)
    data_table.sort_values(by=["Total"], ascending=False, inplace=True)

    if unique_shortname_rows:
        df_unique_shortnames = pd.concat(unique_shortname_rows, ignore_index=True)
    else:
        df_unique_shortnames = pd.DataFrame(columns=df.columns)
    print(data_table[data_table['Known'] != 0])
    return data_table, df_unique_shortnames


@app.route('/unique_sequences', methods=['POST', 'GET'])
def unique_sequences():
    selected_accession = request.form.get('reference_id')
    data_table, df_unieke_shortnames = get_unique_sequences(selected_accession)
    if selected_accession:
        region_data = get_region_data(BASE_PATH)
        if not region_data.empty:
            region_data = region_data[region_data["File"].str.contains(selected_accession)]

        return render_template("unique_sequences.html",
                               data_table=df_unieke_shortnames,
                               selected_accession=selected_accession,
                               region_data=region_data)
    plot = get_known_novel_plot(data_table)
    return render_template("unique_sequences.html",
                           data_table=data_table,
                           plot=plot)


@app.route('/sequence_list', methods=['POST', 'GET'])
def sequence_list():
    df = get_annotation_data()
    regions = df["Region"].unique()
    selected_region = request.form.get("region", regions[0])

    filterd_df = df[(df["Region"] == selected_region) & (df["Segment"].isin(["V", "D", "J"]))]
    segments = filterd_df["Segment"].unique()
    selected_segment = request.form.get("segment", segments[0])
    filterd_df = df[(df["Region"] == selected_region) & (df["Segment"] == selected_segment)]

    grouped_data = filterd_df.groupby(["Target name", "Status"]).size().reset_index(name="Count")
    pivot_data = grouped_data.pivot_table(
        index="Target name",
        columns="Status",
        values="Count",
        fill_value=0
    ).reset_index()
    pivot_data = pivot_data.astype(int, errors="ignore")

    pivot_data["Known"] = pivot_data.get("Known", 0)
    pivot_data["Novel"] = pivot_data.get("Novel", 0)
    pivot_data["Total"] = pivot_data["Known"] + pivot_data["Novel"]
    pivot_data = pivot_data.sort_values(by="Total", ascending=False)
    return render_template("sequence.html",
                           regions=regions,
                           selected_region=selected_region,
                           segments=segments,
                           selected_segment=selected_segment,
                           pivot_data=pivot_data.to_dict('records')
                           )


@app.route('/get_sequence_table', methods=['POST', 'GET'])
def get_sequence_table():
    reference = request.form.get('reference_id')
    df = get_annotation_data()

    filtered_df = df[df["Target name"] == reference][
        ["Sample", "Short name", "Start coord", "End coord", "Status", "SNPs", "Insertions", "Deletions", "Library sequence", "Target sequence", "Library name"]
    ].copy()
    return render_template("sequence_list.html",reference=reference, data_table=filtered_df, zip=zip)


def make_dendrogram(df):
    """
    Build interactive dendrograms per Region and return a dict mapping region -> HTML div.
    Requires:
      - df contains columns: 'Region', 'Sample', 'Target sequence', 'Continent'
      - at least 50 observations per Sample to include in clustering
    """
    region_divs = {}
    for region, df_reg in df.groupby('Region'):

        ctab = pd.crosstab(df['Sample'], df['Target sequence'])
        Z = linkage(ctab.values, method='average', metric='jaccard')
        labels = ctab.index.tolist()
        fig = ff.create_dendrogram(
            ctab.values,
            orientation='left',
            labels=labels,
            linkagefun=lambda x: Z
        )

        for trace in fig.data:
            if trace.mode == 'lines':
                trace.hoverinfo = 'x'
                trace.hovertemplate = 'Jaccard distance: %{x:.2f}<extra></extra>'

        height = max(400, 20 * len(labels))
        fig.update_layout(
            width=1200,
            height=height,
            margin=dict(l=150, r=50, t=50, b=50),
            template="simple_white",
            xaxis_title = "Jaccard distance",
            yaxis_title = "Sample"
        )

        div = plot(
            fig,
            config={
                'displaylogo': False,
                'modeBarButtonsToRemove': [
                    'zoom2d','pan2d','select2d','lasso2d',
                    'resetScale2d','zoomIn2d','zoomOut2d','autoScale2d'
                ],
                'modeBarButtonsToAdd': ['toImage'],
                'scrollZoom': False,
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': f'pca_{region}',
                    'scale': 1
                }
            },
            output_type='div',
            include_plotlyjs=True
        )
        region_divs[region] = div

    return region_divs


@app.route('/get_dendrogram', methods=['POST', 'GET'])
def get_dendrogram():
    df = get_annotation_data()
    plot_divs = make_dendrogram(df)
    return render_template("phylogenetic_tree.html", dendrogram_divs=plot_divs)


def make_pca_figures(df):
    """
    Return a dict mapping each Region to its Plotly PCA <div> string.
    If 'Population' is missing, colors are by 'Sample'.
    """
    has_pop = 'Population' in df.columns
    region_divs = {}

    for region, df_reg in df.groupby('Region'):
        feat = (
            df_reg
            .groupby(['Sample','Target name'])
            .size()
            .unstack(fill_value=0)
        )

        labels = pd.DataFrame({'Sample': feat.index})
        if has_pop:
            pop = df_reg[['Sample','Population']].drop_duplicates()
            labels = labels.merge(pop, on='Sample', how='left')
        else:
            labels['Population'] = labels['Sample']
        labels.set_index('Sample', inplace=True)

        X = StandardScaler().fit_transform(feat)
        pca = PCA(n_components=2, random_state=42)
        pcs = pca.fit_transform(X)

        plot_df = pd.DataFrame({
            'PC1': pcs[:,0],
            'PC2': pcs[:,1],
            'Sample': feat.index,
            'Population': labels['Population']
        })

        color_col = 'Population' if has_pop else 'Sample'
        hover = ['Sample','Population'] if has_pop else ['Sample']
        fig = px.scatter(
            plot_df,
            x='PC1', y='PC2',
            color=color_col,
            hover_data=hover,
            #title=f"PCA — Region: {region}",
            template='simple_white'
        )
        fig.update_traces(marker_size=10)

        div = plot(
            fig,
            config={
                'displaylogo': False,
                'modeBarButtonsToRemove': [
                    'zoom2d','pan2d','select2d','lasso2d',
                    'resetScale2d','zoomIn2d','zoomOut2d','autoScale2d'
                ],
                'modeBarButtonsToAdd': ['toImage'],
                'scrollZoom': False,
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': f'pca_{region}',
                    'scale': 1
                }
            },
            output_type='div',
            include_plotlyjs=True
        )

        region_divs[region] = div

    return region_divs

@app.route('/get_pca', methods=['POST', 'GET'])
def get_pca():
    df = get_annotation_data()
    plot_divs = make_pca_figures(df)
    return render_template("pca.html", pca_divs=plot_divs)


def chunk_sequence(sequence, chunk_size=100):
    return [sequence[i:i + chunk_size] for i in range(0, len(sequence), chunk_size)]


def format_fixed_width(text, width=100):
    return text.ljust(width)[:width]


@app.route('/msa_sequence', methods=['POST'])
def msa_sequence():
    reference =request.form.get('reference')
    allele =request.form.get('allele')

    width  = max(len(reference), len(allele)) + 5
    reference = format_fixed_width(reference, width)
    allele = format_fixed_width(allele, width)

    query = request.form.get('query')
    subject = request.form.get('subject')

    query_chunks = chunk_sequence(query, 80)
    subject_chunks = chunk_sequence(subject, 80)

    formatted_query_lines = []
    match_lines = []
    class_chunks = []
    for q_chunk, s_chunk in zip(query_chunks, subject_chunks):
        formatted_q_chunk = ""
        match_line = ""
        class_chunk = []
        for index, (q, s) in enumerate(zip(q_chunk, s_chunk)):
            if ((index % 10) == 0) and (index != 0):
                formatted_q_chunk += " "
                match_line += " "
                class_chunk.append(" ")
            if q == s:
                match_line += "-"
                class_chunk.append("")
            elif q == "-":
                match_line += s
                class_chunk.append("deletion")
            elif s == "-":
                match_line += "*"
                class_chunk.append("insertion")
            else:
                match_line += s
                class_chunk.append("snp")

            formatted_q_chunk += q

        formatted_query_lines.append(formatted_q_chunk)
        match_lines.append(match_line)
        class_chunks.append(class_chunk)

    query_chunks = [query_chunk for query_chunk in query_chunks]
    return render_template(
            "msa.html",
            reference=reference,
            allele=allele,
            query_chunks=formatted_query_lines,
            match_lines=match_lines,
            class_chunks=class_chunks,
            zip=zip
        )


def get_pivot_table(df: pd.DataFrame) -> pd.DataFrame:
    pivot_table = df.pivot_table(
        index="Sample",
        columns="Region",
        aggfunc="size",
        fill_value=0
    )
    return pivot_table


@app.route('/get_report', methods=['POST', 'GET'])
def get_report():
    df = get_annotation_data()

    pivot_table_known = get_pivot_table(df[df["Status"] == "Known"])
    vdj_plot_known = get_vdj_plot(pivot_table_known, "known")

    pivot_table_novel = get_pivot_table(df[df["Status"] == "Novel"])
    vdj_plot_novel = get_vdj_plot(pivot_table_novel, "novel")
    return render_template("report.html", vdj_plot_known=vdj_plot_known, vdj_plot_novel=vdj_plot_novel)


@app.route('/get_count_contigs', methods=['POST', 'GET'])
def get_count_contigs():
    region_data = get_region_data(BASE_PATH)
    region_data["Sample"] = (region_data["File"].str.replace(".fasta", "", regex=False).str.replace(".fa", "", regex=False))
    region_data_pivot = (
        region_data
        .pivot_table(
            index='Sample',
            columns='Segment',
            values='Contig counts',
            aggfunc='sum',
            fill_value=0
        )
        .reset_index()
    )
    return render_template("contigs.html", pivot_data=region_data_pivot.to_dict('records'))


@app.route('/get_scaffold_figure', methods=['POST', 'GET'])
def get_scaffold_figure():
    selected_sample = request.form.get('selected_sample')
    selected_contig = request.form.get('selected_contig')

    scaffold_data = get_scaffold_data(BASE_PATH)

    contigs = scaffold_data[
        (scaffold_data['File'] == selected_sample) &
        (scaffold_data['Scaffold'] == selected_contig)
        ]["Contig"].to_list()

    return send_file(
        contigs,
        mimetype="text/plain",
        as_attachment=True,
        download_name=f"{selected_sample}.bed"
    )



@app.route('/get_sample_table', methods=['POST', 'GET'])
def get_sample_table():
    selected_sample = request.form.get('selected_sample')
    selected_region = request.form.get('region')
    df = get_annotation_data()
    regions = df["Region"].unique().tolist()

    vdj_plot = None
    if selected_sample:
        df_sample = df[df["Sample"] == selected_sample].copy()
        regions = df_sample["Region"].unique().tolist()
        selected_region = request.form.get('region', request.form.get('selected_region'))

        if selected_region is None or selected_region not in regions:
            selected_region = regions[0]

        df_filtered = df_sample[df_sample["Region"] == selected_region].copy()

        region_data = get_region_data(BASE_PATH)
        plot_dict = {}

        if not region_data.empty:
            region_data = region_data[region_data["File"].str.contains(selected_sample)]

            region_data["5 coords"] = pd.to_numeric(region_data["5 coords"], errors="coerce").fillna(0).astype(int)
            region_data["3 coords"] = pd.to_numeric(region_data["3 coords"], errors="coerce").fillna(0).astype(int)
            region_data = region_data.iloc[:, 2:]

            for (region_name, contig), row_df in region_data.groupby(["Region", "Contig"]):
                plot_df_filtered = df[
                    (df["Sample"] == selected_sample) &
                    (df["Region"] == region_name) &
                    (df["Contig"] == contig)
                    ].copy()

                if not plot_df_filtered.empty:
                    script, div = generate_bokeh_segment_div_from_df(plot_df_filtered)

                    plot_dict[(region_name, contig)] = {
                        "script": script,
                        "div": div
                    }

        segments_plot = get_pi_plot(df_filtered, "Segment", "Segment")
        known_novel_plot = get_pi_plot(df_filtered, "Status", "Status")
        function_plot = get_pi_plot(df_filtered, "Function", "Function")
        function_messenger_plot = get_pi_plot(df_filtered[df_filtered["Function"] == "ORF"], "Function_messenger", "ORF reason")
        function_messenger_plot2 = get_pi_plot(df_filtered[df_filtered["Function"] == "pseudo"], "Function_messenger", "Pseudo reason")


        region_viewer = get_dna_vieuwer_plot(df_filtered)

        return render_template("annotation_sample.html",
                               data_table=df_filtered,
                               region_data=region_data,
                               selected_sample=selected_sample,
                               plots=plot_dict,
                               regions=regions,
                               selected_region=selected_region,
                               segments_plot=segments_plot,
                               known_novel_plot=known_novel_plot,
                               function_plot=function_plot,
                               function_messenger_plot=function_messenger_plot,
                               function_messenger_plot2=function_messenger_plot2,
                               region_viewer=region_viewer
                               )
    elif selected_region:
        df = df[df["Region"] == selected_region].copy()
        vdj_plot = get_vdj_plot(df.pivot_table(index="Sample",columns=["Segment"], aggfunc="size", fill_value=0), f"{selected_region}")

    grouped_data = df.groupby(["Sample", "Status"]).size().reset_index(name="Count")

    pivot_data = grouped_data.pivot_table(
        index="Sample",
        columns="Status",
        values="Count",
        fill_value=0
    ).reset_index()
    pivot_data = pivot_data.astype(int, errors="ignore")

    pivot_data["Known"] = pivot_data.get("Known", 0)
    pivot_data["Novel"] = pivot_data.get("Novel", 0)
    pivot_data["Total"] = pivot_data["Known"] + pivot_data["Novel"]
    pivot_data = pivot_data.sort_values(by="Total", ascending=False)


    return render_template("samples.html", pivot_data=pivot_data.to_dict('records'), regions=regions, selected_region=selected_region, vdj_plot=vdj_plot)



@app.route('/get_rss', methods=['POST', 'GET'])
def get_rss():
    df = get_annotation_data()

    regions = df["Region"].unique()
    selected_region = request.form.get("region", regions[0])

    filterd_df = df[(df["Region"] == selected_region) & (df["Segment"].isin(["V", "D", "J"]))]
    segments = filterd_df["Segment"].unique()
    selected_segment = request.form.get("segment", segments[0])

    count_df = filterd_df[filterd_df["Segment"] == selected_segment]
    rss_plot = get_rss_plot(count_df)

    rss_data = {}
    if RSS_PATH.is_dir():
        meme_pattern = f"meme_output/{selected_region}{selected_segment}_*"
        fimo_pattern = f"fimo_output/{selected_region}{selected_segment}_*"

        meme_submappen = list(RSS_PATH.glob(meme_pattern))
        fimo_submappen = list(RSS_PATH.glob(fimo_pattern))

        for meme_submap in meme_submappen:
            meme_logo = meme_submap / "logo1.png"
            if meme_logo.exists():
                with open(meme_logo, "rb") as img_file:
                    encoded_image = base64.b64encode(img_file.read()).decode('utf-8')
                rss_data[meme_submap.name] = {
                    'meme_logo': encoded_image
                }

        for fimo_submap in fimo_submappen:
            fimo_txt = fimo_submap / "fimo.tsv"
            if fimo_txt.exists():
                fimo_data = []
                with open(fimo_txt, "r", encoding="utf-8") as txt_file:
                    headers = txt_file.readline().strip().split("\t")
                    selected_headers = [headers[i] for i in [2, 6, 7, 8, 9]]
                    for line in txt_file:
                        lines = line.strip().split("\t")
                        if len(lines) < 2:
                            continue
                        selected_lines = [lines[i] for i in [2, 6, 7, 8, 9]]

                        fimo_data.append({
                            header: value for header, value in zip(selected_headers, selected_lines)
                        })

                if fimo_submap.name in rss_data:
                    rss_data[fimo_submap.name]['fimo_data'] = fimo_data
                else:
                    rss_data[fimo_submap.name] = {
                        'fimo_data': fimo_data
                    }

    return render_template("rss.html",
                           segments=segments,
                           selected_segment=selected_segment,
                           regions=regions,
                           selected_region=selected_region,
                           rss_plot=rss_plot,
                           rss_data=rss_data
                           )


@app.route("/download_xlsx")
def download_xlsx():
    df = get_annotation_data()
    selected_sample = request.args.get("reference_id", "").strip()
    region = request.args.get("region", "").strip()

    if selected_sample and region:
        df = df[(df["Sample"] == selected_sample) & (df["Region"] == region)]

    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        df.to_excel(writer, index=False, sheet_name="Sample Data")
    output.seek(0)

    return send_file(
        output,
        mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        as_attachment=True,
        download_name="sample_data.xlsx"
    )


@app.route("/download_gtf")
def download_gtf():
    df = get_annotation_data()
    selected_sample = request.args.get("reference_id", "").strip()
    region = request.args.get("region", "").strip()

    df = df[(df["Sample"] == selected_sample) & (df["Region"] == region)]
    gtf_path = Path(f"/tmp/{selected_sample}_{region}.gtf")

    with open(gtf_path, "w") as f:
        for _, row in df.iterrows():
            attribute_field = f'gene_id "{row["Target name"]}"'
            line = [
                row["Contig"],              #seqname
                "VDJ-Insights",             #source
                "Gene",                     #feature
                str(row["Start coord"]),    #start
                str(row["Start coord"]),    #end
                "0",                        #score
                row["Strand"],              #strand
                ".",                        #frame
                attribute_field             #attribute
            ]
            f.write("\t".join(line) + "\n")

    return send_file(
        gtf_path,
        mimetype="text/plain",
        as_attachment=True,
        download_name=gtf_path.name
    )


@app.route("/download_fasta")
def download_fasta():
    df = get_annotation_data()
    selected_sample = request.args.get("reference_id", "").strip()
    region = request.args.get("region", "").strip()

    fasta_path = Path(df[(df["Sample"] == selected_sample) & (df["Region"] == region)]["Path"].unique()[0])
    if not fasta_path.exists():
        abort(404, description="Fasta error")

    return send_file(
        fasta_path,
        mimetype="text/plain",
        as_attachment=True,
        download_name=fasta_path.name
    )


@app.route("/download_bed")
def download_bed():
    df = get_annotation_data()
    selected_sample = request.args.get("reference_id", "").strip()
    region = request.args.get("region", "").strip()

    df = df[(df["Sample"] == selected_sample) & (df["Region"] == region)]

    bed_columns = ["Contig", "Start coord", "End coord", "Short name", "Strand"]

    output = io.BytesIO()
    df[bed_columns].to_csv(output, sep="\t", index=False, header=False)
    output.seek(0)


    return send_file(
            output,
            mimetype="text/plain",
            as_attachment=True,
            download_name=f"{selected_sample}_{region}.bed"
            )



@app.route('/download_sequences', methods=['POST'])
def download_sequences(selected_ids):
    sequences = get_sequences(LIBRARY_PATH / 'library.fasta')

    selected_sequences = [seq for seq in sequences if seq['id'] in selected_ids]

    fasta_content = ""
    for seq in selected_sequences:
        fasta_content += f">{seq['id']}\n{seq['sequence']}\n"

    response = make_response(fasta_content)
    response.headers["Content-Disposition"] = "attachment; filename=selected_sequences.fasta"
    response.headers["Content-Type"] = "text/plain"
    return response


@app.route('/remove_sequences', methods=['POST'])
def remove_sequences(selected_ids):
    sequences = get_sequences(LIBRARY_PATH / 'library.fasta')
    selected_sequences = [seq for seq in sequences if seq['id'] not in selected_ids]

    library_file = Path(LIBRARY_PATH / "library.fasta")

    old_library_dir = Path(LIBRARY_PATH / "old")
    old_library_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archived_file = old_library_dir / f"library_{timestamp}.fasta"
    shutil.copy(library_file, archived_file)

    with open(library_file, "w") as file:
        for seq in selected_sequences:
            file.write(f">{seq['id']}\n{seq['sequence']}\n")

    flash("Sequences removed from library.", "info")
    clear_cache()
    return redirect(request.referrer)


def add_to_library(selected_ids: list):
    df = get_annotation_data()
    mask = df["Short name"].isin(selected_ids) & (df["Status"] == "Novel")
    df.loc[mask, "Status"] = "Known"

    df_known = df[df["Status"] == "Known"]
    df_novel = df[df["Status"] == "Novel"]

    file_known = Path(BASE_PATH / "annotation/annotation_report_known_rss.xlsx")
    file_novel = Path(BASE_PATH / "annotation/annotation_report_novel_rss.xlsx")

    df_known.to_excel(file_known, index=False)
    df_novel.to_excel(file_novel, index=False)

    selected_sequences = df[df["Short name"].isin(selected_ids)][["Short name", "Target sequence"]].drop_duplicates(subset=["Short name"])
    selected_sequences["Target sequence"] = selected_sequences["Target sequence"].str.replace("_", "", regex=False)

    library_file = Path(LIBRARY_PATH / "library.fasta")

    old_library_dir = Path(LIBRARY_PATH / "old")
    old_library_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archived_file = old_library_dir / f"library_{timestamp}.fasta"
    shutil.copy(library_file, archived_file)

    with open(library_file, "a") as f:
        for _, row in selected_sequences.iterrows():
            f.write(f">{row['Short name']}\n{row['Target sequence']}\n")


@app.route('/add_novel', methods=['POST', 'GET'])
def add_novel():
    selected_ids = request.form.getlist('selected_sequences')
    if len(selected_ids) >= 1:
        add_to_library(selected_ids)
        clear_cache()
        flash(f"Successfully added {len(selected_ids)} segment(s) to the library!", "success")

    selected_filter_amount = request.form.get("filter_amount", type=int)

    df = get_annotation_data()
    df = df[df["Status"] == "Novel"]
    if not df.empty:
        df_unique = df.drop_duplicates(subset=["Sample", "Short name"])

        grouped_data = df_unique.groupby(["Short name"]).size().reset_index(name="Count")
        extra_columns = ["% identity", "SNPs", "Insertions", "Deletions"]
        grouped_extra = df_unique.groupby("Short name")[extra_columns].mean().reset_index()
        grouped_extra[["SNPs", "Insertions", "Deletions"]] = grouped_extra[["SNPs", "Insertions", "Deletions"]].astype(int)
        grouped_extra["% identity"] = grouped_extra["% identity"].round(2)

        pivot_data = pd.merge(grouped_data, grouped_extra, on="Short name")
        pivot_data = pivot_data.sort_values(by="Count", ascending=False)

        if selected_filter_amount:
                pivot_data = pivot_data[pivot_data["Count"] >= selected_filter_amount]

        plot_count_distribution = plot_count_distribution_bar_chart(pivot_data)
        return render_template("add_to_library.html",
                               plot_count_distribution=plot_count_distribution,
                               pivot_data=pivot_data.to_dict('records'),
                               selected_filter_amount=selected_filter_amount
                               )
    return render_template("add_to_library.html")


@cache.cached(timeout=600, key_prefix="data1")
@app.route('/get_library', methods=['POST', 'GET'])
def get_library():
    action = request.form.get("action")
    selected_ids = request.form.getlist("selected_sequences")

    if action == "remove":
        return remove_sequences(selected_ids)
    elif action == "download":
        return download_sequences(selected_ids)

    sequences = get_sequences(LIBRARY_PATH / 'library.fasta')
    library_info = open_json(LIBRARY_PATH / 'library_info.json')

    if library_info:
        tilte = f"""Library {library_info["type"]} {library_info["species"]} ({library_info["set_release"]})"""
    else:
        tilte = f"""Own library"""

    regions = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ", "TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", "TRGV", "TRGJ", "TRDV", "TRDD", "TRDJ"]

    for seq in sequences:
        for region in regions:
            if region in seq["id"]:
                seq["region"] = region
                break
        else:
            seq["region"] = "Unknown"

    result = [{"name": region, "entries": count} for region in regions if (count := sum(region in seq["id"] for seq in sequences)) > 0]

    names = [item["name"] for item in result]
    entries = [item["entries"] for item in result]

    library_plot = get_library_plot(names, entries, tilte)
    library_violin = get_library_violin_plot(sequences)

    return render_template("library.html",
                           sequences=sequences,
                           library_plot=library_plot,
                           library_violin=library_violin
                           )


@app.route('/get_about', methods=['POST', 'GET'])
def get_about():
    return render_template("about.html")


@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html")


def clear_cache():
    cache.delete("data1")


#moet nog naar import_data maar geeft foutmerling met betrekking tot cache
@cache.cached(timeout=600, key_prefix="data1")
def get_annotation_data() -> pd.DataFrame:
    df = pd.read_excel(BASE_PATH / "annotation/annotation_report_all.xlsx")
    return df

