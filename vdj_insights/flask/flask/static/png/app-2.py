import configparser
import datetime
from pathlib import Path
from flask_caching import Cache
from dna_features_viewer import GraphicFeature, GraphicRecord
from itertools import cycle

from bokeh.models import Legend, LegendItem, NumeralTickFormatter, HoverTool
from bokeh.embed import file_html
from bokeh.resources import CDN
from Bio import SeqIO
import os
import json
import plotly.express as px
import matplotlib
matplotlib.use('agg')
from plotly.offline import plot
import pandas as pd
import base64

from flask import Flask, render_template, request, redirect, make_response, session


config = configparser.ConfigParser()
config.read('config.ini')

app = Flask(__name__)
app.secret_key = os.urandom(12)
app.config["CACHE_TYPE"] = "SimpleCache"
app.config["CACHE_DEFAULT_TIMEOUT"] = 3600
cache = Cache(app)

@app.route('/', methods=['POST', 'GET'])
def start():
    return render_template("home.html")


def get_pivot_table(df: pd.DataFrame) -> pd.DataFrame:
    pivot_table = df.pivot_table(
        index="Sample",
        columns="Region",
        aggfunc="size",
        fill_value=0
    )
    return pivot_table


def get_vdj_plot(pivot_table: pd.DataFrame) -> dict:
    fig = px.bar(
        pivot_table,
        barmode="stack",
        labels={"value": "Aantal sequenties", "Sample": "Samples", "Region": "Regio's"},
        title="Aantal sequenties per sample en regio",
        color_discrete_sequence=px.colors.qualitative.Set2
    )
    fig.update_traces(
        textfont_size=12,
        textangle=0,
        textposition="outside",
        cliponaxis=False
    )
    fig.update_traces(
        hoverinfo="y",
        hovertemplate="<b>Total:</b> %{y}<extra></extra>"
    )
    fig.update_layout(
        dragmode=False,
        template="plotly_white",
        yaxis=dict(showgrid=False),
    )

    plot_div = plot(
        fig,
        config={
            'displayModeBar': False,
            'scrollZoom': False
        },
        output_type='div',
        include_plotlyjs=False
    )
    return plot_div


@app.route('/get_report', methods=['POST', 'GET'])
def get_report():
    df = get_annotation_data("data1")
    pivot_table_known = get_pivot_table(df[df["Status"] == "Known"])
    vdj_plot_known = get_vdj_plot(pivot_table_known)

    pivot_table_novel = get_pivot_table(df[df["Status"] == "Novel"])
    vdj_plot_novel = get_vdj_plot(pivot_table_novel)

    vdj_plot = get_vdj_plot(df.pivot_table(
        index="Sample",
        columns=["Segment"],
        aggfunc="size",
        fill_value=0
    ))
    return render_template("report.html", vdj_plot_known=vdj_plot_known, vdj_plot_novel=vdj_plot_novel, vdj_plot=vdj_plot)


@cache.cached(timeout=600, key_prefix="data1")
def get_annotation_data(file_path: str) -> pd.DataFrame:
    file_known = Path("data1/annotation/annotation_report_known.xlsx")
    file_novel = Path("data1/annotation/annotation_report_novel.xlsx")

    dataframes = []
    if file_known.exists():
        dataframes.append(pd.read_excel(file_known))
    if file_novel.exists():
        dataframes.append(pd.read_excel(file_novel))
    df = pd.concat(dataframes, ignore_index=True)
    return df


def get_dna_vieuwer_plot(df_filtered: pd.DataFrame, selected_accession: str, selected_region: str) -> dict:
    if not df_filtered.empty:
        named_colors = [
            "lightblue", "lightcoral", "lightgreen", "lightpink", "lightyellow",
            "lightgray", "lightcyan", "lavender", "peachpuff", "khaki"
        ]
        unique_segments = df_filtered["Segment"].unique()
        segments_colors = {seg: color for seg, color in zip(unique_segments, cycle(named_colors))}

        features = [
            GraphicFeature(
                start=row["Start coord"],
                end=row["End coord"],
                strand=int(str(row["Strand"]) + "1"),
                color=segments_colors[row["Segment"]],
                label=f"{row['Short name']}"
            )
            for _, row in df_filtered.iterrows()
        ]

        start_coord = df_filtered["Start coord"].min() - 100
        end_coord = df_filtered["End coord"].max() - start_coord + 100

        record = GraphicRecord(
            sequence="ATCG",
            first_index=start_coord,
            sequence_length=end_coord,
            features=features
        )
        bokeh_plot = record.plot_with_bokeh(figure_width=80, figure_height=40)
        bokeh_plot.title.text = f"{selected_accession} - {selected_region}"
        bokeh_plot.title.text_font_size = "20pt"
        bokeh_plot.title.align = "center"
        from bokeh.layouts import gridplot

        bokeh_plot.xaxis[0].formatter = NumeralTickFormatter(format="0")
        hover = HoverTool(tooltips=[("Short name", "@label"), ("Sample", selected_accession)])
        bokeh_plot.add_tools(hover)

        grid = gridplot([[bokeh_plot]], sizing_mode="stretch_both")
        graph = file_html(grid, CDN, "Region Viewer")
        #graph = file_html(bokeh_plot, CDN, "Region Viewer")
        return graph

@app.route('/get_samples', methods=['POST', 'GET'])
def get_samples():
    df = get_annotation_data("data1")
    samples = df["Sample"].unique()
    regions = df["Region"].unique()

    selected_accession = request.form.get("sample", samples[0])
    selected_region = request.form.get("region", regions[0])

    df_filtered = df[(df["Sample"] == selected_accession) & (df["Region"] == selected_region)].copy()

    df_filtered.sort_values(by="Start coord", inplace=True)
    region_viewer = get_dna_vieuwer_plot(df_filtered, selected_accession, selected_region)

    print(df_filtered.groupby("Segment").size())
    fig = px.pie(df_filtered.groupby("Segment").size())
    print(fig)
    return render_template("samples_report.html",
                           samples=samples,
                           selected_sample=selected_accession,
                           regions=regions,
                           selected_region=selected_region,
                           region_viewer=region_viewer,
                           data_table=df_filtered)


@app.route('/get_rss', methods=['POST', 'GET'])
def get_rss():
    df = get_annotation_data("data1")
    regions = df["Region"].unique()
    segments = df["Segment"].unique()

    selected_region = request.form.get("region", regions[0])
    selected_segment = request.form.get("segment", segments[0])
    cwd = Path.cwd()
    rss_dir = cwd / "data1" / "RSS"

    rss_data = {}
    if rss_dir.is_dir():
        meme_pattern = f"meme_output/{selected_region}{selected_segment}_*"
        fimo_pattern = f"fimo_output/{selected_region}{selected_segment}_*"

        meme_submappen = list(rss_dir.glob(meme_pattern))
        fimo_submappen = list(rss_dir.glob(fimo_pattern))

        for meme_submap in meme_submappen:
            meme_logo = meme_submap / "logo1.png"
            if meme_logo.exists():
                with open(meme_logo, "rb") as img_file:
                    encoded_image = base64.b64encode(img_file.read()).decode('utf-8')
                rss_data[meme_submap.name] = {
                    'meme_logo': encoded_image
                }

        for fimo_submap in fimo_submappen:
            fimo_txt = fimo_submap / "fimo.txt"
            if fimo_txt.exists():
                fimo_data = []
                with open(fimo_txt, "r", encoding="utf-8") as txt_file:
                    headers = txt_file.readline().strip().split("\t")
                    selected_headers = [headers[i] for i in [1, 5, 6, 7]]

                    for line in txt_file:
                        lines = line.strip().split("\t")
                        selected_lines = [lines[i] for i in [1, 5, 6, 7]]

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
                           rss_data=rss_data
                           )


def get_sequences(file_name: str) -> list:
    sequences = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append({
            "id": record.id,
            "sequence_length": len(record.seq),
            "sequence": str(record.seq)
        })
    return sequences


def get_library_info(file_name: str) -> dict:
    if os.path.exists(file_name):
        with open(file_name) as file:
            data = json.load(file)
        return data


def get_library_violin_plot(sequences: list) -> dict:
    df = pd.DataFrame(sequences)

    fig = px.violin(
        df,
        y="sequence_length",
        box=True,
        points="all",
        labels={"sequence_length": "Sequence length"},
        title="Distribution of Sequence Lengths",
        template="simple_white",
        color_discrete_sequence=px.colors.qualitative.Set2
    )

    fig.update_layout(
        dragmode=False,
        yaxis=dict(showgrid=True),
        xaxis=dict(showgrid=False),
    )

    plot_div = plot(
        fig,
        config={
            'displayModeBar': False,
            'scrollZoom': False
        },
        output_type='div',
        include_plotlyjs=False
    )

    return plot_div


def get_library_plot(library_info: dict) -> dict:
    names = [item["name"].strip(".fasta") for item in library_info.get("fasta_files", [])]
    entries = [item["entries"] for item in library_info.get("fasta_files", [])]

    fig = px.bar(
        x=names,
        y=entries,
        labels={"x": "Region", "y": "Number of sequences"},
        title=f"""Library {library_info["type"]} {library_info["species"]} ({library_info["set_release"]})""",
        text=entries,
        template = 'simple_white',
        color_discrete_sequence=px.colors.qualitative.Set2

    )

    fig.update_traces(
        textfont_size=12,
        textangle=0,
        textposition="outside",
        cliponaxis=False
    )

    fig.update_traces(
        hoverinfo="y",
        hovertemplate="<b>Total:</b> %{y}<extra></extra>"
    )
    fig.update_layout(
        dragmode=False,
        template="plotly_white",
        yaxis=dict(showgrid=False),
    )

    plot_div = plot(
        fig,
        config={
            'displayModeBar': False,
            'scrollZoom': False
        },
        output_type='div',
        include_plotlyjs=False
    )

    return plot_div


@app.route('/download_sequences', methods=['POST'])
def download_sequences():
    selected_ids = request.form.getlist('selected_sequences')
    print(selected_ids)
    sequences = get_sequences('data1/library/library.fasta')

    selected_sequences = [seq for seq in sequences if seq['id'] in selected_ids]

    fasta_content = ""
    for seq in selected_sequences:
        fasta_content += f">{seq['id']}\n{seq['sequence']}\n"

    response = make_response(fasta_content)
    response.headers["Content-Disposition"] = "attachment; filename=selected_sequences.fasta"
    response.headers["Content-Type"] = "text/plain"
    return response


@cache.cached(timeout=600, key_prefix="data1")
@app.route('/get_library', methods=['POST', 'GET'])
def get_library():
    sequences = get_sequences('data1/library/library.fasta')

    library_info = get_library_info('data1/library/library_info.json')
    library_plot = get_library_plot(library_info)
    library_violin = get_library_violin_plot(sequences)


    return render_template("library.html", sequences=sequences, library_plot=library_plot, library_violin=library_violin, **library_info)


@app.route('/get_about', methods=['POST', 'GET'])
def get_about():
    return render_template("about.html")


@app.errorhandler(401)
def internal_server_error(error):
    return redirect("/login")


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5001, debug=True)
