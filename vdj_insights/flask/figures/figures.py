from itertools import cycle

from bokeh.models import Legend, LegendItem, NumeralTickFormatter, HoverTool
from bokeh.embed import file_html
from bokeh.resources import CDN

import plotly.express as px
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import io
import base64

from dna_features_viewer import GraphicFeature, GraphicRecord
from plotly.offline import plot

matplotlib.use('agg')


def get_known_novel_plot(df: dict) -> str:
    fig = px.bar(
        df[df["Sample"] != "Total"],
        x="Sample",
        y=["Known", "Novel"],
        barmode="stack",
        labels={"value": "Count", "variable": "Type", "Sample": "Samples"},
        title="Unique sequences",
        color_discrete_sequence=px.colors.qualitative.Set2
    )

    fig.update_traces(
        textfont_size=12,
        textangle=0,
        textposition="outside",
        cliponaxis=False,
        hoverinfo="y",
        hovertemplate="<b>Total:</b> %{y}<extra></extra>"
    )

    fig.update_layout(
        xaxis_tickangle=-45,
        template="plotly_white",
        yaxis=dict(title="Aantal segmenten", showgrid=True),
        xaxis=dict(title="Samples"),
        legend_title="Segment Type",
        margin=dict(l=40, r=40, t=40, b=120),
        dragmode=False
    )

    plot_div = plot(
        fig,
        config={
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                       'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )
    return plot_div


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
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                       'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )
    return plot_div


def get_sample_rss_plot(df_filtered: pd.DataFrame) -> dict:
    if not df_filtered.empty:
        rss_counts = df_filtered.groupby(['Segment', 'Status']).agg({
            "3'-RSS": ['sum', lambda x: len(x) - sum(x)],
            "5'-RSS": ['sum', lambda x: len(x) - sum(x)]
        }).reset_index()

        rss_counts.columns = ['Segment', 'Status', "3'-RSS True", "3'-RSS False", "5'-RSS True", "5'-RSS False"]
        rss_counts_melted = rss_counts.melt(
            id_vars=['Segment', 'Status'],
            value_vars=["3'-RSS True", "3'-RSS False", "5'-RSS True", "5'-RSS False"],
            var_name='RSS Type',
            value_name='Count'
        )

        rss_counts_melted[['RSS Type', 'Value']] = rss_counts_melted['RSS Type'].str.split(' ', expand=True)
        rss_counts_melted = rss_counts_melted[rss_counts_melted['Count'].notna() & (rss_counts_melted['Count'] > 0)]
        fig = px.bar(
            rss_counts_melted,
            x='Segment',
            y='Count',
            color='Value',
            facet_col='RSS Type',
            facet_row='Status',
            text='Count',
            barmode='group',
            title=f"RSS",
            labels={'Segment': 'Segment', 'Count': 'Count', 'Value': 'RSS', 'Status': 'Status'},
            color_discrete_sequence=px.colors.qualitative.Set2
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

        fig.update_yaxes(showgrid=False)
        fig.update_xaxes(showgrid=False)

        fig.update_traces(
            textfont_size=12,
            textangle=0,
            textposition="outside",
            cliponaxis=False,
            hoverinfo="y",
            hovertemplate="<b>Total:</b> %{y}<br><b>Segment:</b> %{x}<extra></extra>"
        )

        fig.update_layout(
            dragmode=False,
            title=dict(
                text=f"RSS",
                font=dict(size=18, family="Arial", color="black"),
                x=0.5
            ),
            template="plotly_white",
            yaxis=dict(showgrid=False),
        )
        plot_div = plot(
            fig,
            config={
                'displaylogo': False,
                'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                           'zoomOut2d', 'autoScale2d'],
                'modeBarButtonsToAdd': ['toImage'],
                'scrollZoom': False,
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': 'custom_image',
                    'scale': 1
                }
            },
            output_type='div',
            include_plotlyjs=True
        )
        return plot_div


def get_rss_plot(count_df: pd.DataFrame) -> dict:
    if not count_df.empty:
        counts_3rss_status = count_df.groupby("Status")["3'-RSS"].value_counts(dropna=False).reset_index(name="Count")
        counts_3rss_status["Column"] = "3'-RSS"
        counts_3rss_status.rename(columns={"3'-RSS": "Value"}, inplace=True)

        counts_5rss_status = count_df.groupby("Status")["5'-RSS"].value_counts(dropna=False).reset_index(name="Count")
        counts_5rss_status["Column"] = "5'-RSS"
        counts_5rss_status.rename(columns={"5'-RSS": "Value"}, inplace=True)

        combined_counts_status = pd.concat([counts_3rss_status, counts_5rss_status], ignore_index=True)
        combined_counts_status["Value"] = combined_counts_status["Value"].replace({1.0: True, 0.0: False})

        fig = px.bar(
            combined_counts_status,
            x="Column",
            y="Count",
            color="Value",
            barmode="group",
            text="Count",
            facet_col="Status",
            facet_col_wrap = 2,
            title="Found RSS sequence motifs",
            labels={"Column": "RSS Type", "Count": "Count"},
            color_discrete_sequence = px.colors.qualitative.Set2
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

        fig.update_traces(
            textfont_size=12,
            textangle=0,
            textposition="outside",
            cliponaxis=False,
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
                'displaylogo': False,
                'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                           'zoomOut2d', 'autoScale2d'],
                'modeBarButtonsToAdd': ['toImage'],
                'scrollZoom': False,
                'toImageButtonOptions': {
                    'format': 'svg',
                    'filename': 'custom_image',
                    'scale': 1
                }
            },
            output_type='div',
            include_plotlyjs=True
        )
        return plot_div

def get_shared_shortname_heatmap_plotly_simple(df: pd.DataFrame):
    # Maak een matrix waarin de rijen 'Short name' zijn en de kolommen 'Sample' bevatten (transpose matrix)
    sample_shortname_matrix = df.pivot_table(index="Sample", columns="Short name", aggfunc="size", fill_value=0)

    fig = px.imshow(
        sample_shortname_matrix.values,
        labels=dict(x="Short name", y="Sample", color="Aantal gedeelde Short names"),
        x=sample_shortname_matrix.columns.tolist(),
        y=sample_shortname_matrix.index.tolist(),
        color_continuous_scale="Blues",
    )

    fig.update_layout(
        title="Aantal gedeelde 'Short name' tussen samples",
        xaxis_title="Short name",
        yaxis_title="Sample",
        template="plotly_white",
        coloraxis_showscale=False,
    )

    plot_div = plot(
        fig,
        config={
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                       'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )

    return plot_div

def get_shared_shortname_heatmap_plotly_simple2(df: pd.DataFrame):
    # Maak een matrix waarin de rijen Samples zijn en kolommen 'Short name' bevatten
    sample_shortname_matrix = df.pivot_table(index="Sample", columns="Short name", aggfunc="size", fill_value=0)

    shared_matrix = sample_shortname_matrix @ sample_shortname_matrix.T

    fig = px.imshow(
        shared_matrix.values,
        labels=dict(x="Sample", y="Sample", color="Aantal gedeelde Short names"),
        x=shared_matrix.columns.tolist(),
        y=shared_matrix.index.tolist(),
        color_continuous_scale="Blues",
        #text_auto=True
    )


    fig.update_layout(
        title="Aantal gedeelde 'Short name' tussen samples",
        xaxis_title="Sample",
        yaxis_title="Sample",
        template="plotly_white",
        autosize=True,
        width=2000,
        height=2000
    )

    plot_div = plot(
        fig,
        config={
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                       'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )

    return plot_div


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
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                       'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )
    return plot_div


def get_library_plot(names: list, entries: list, title: str) -> dict:
    fig = px.bar(
        x=names,
        y=entries,
        labels={"x": "Region", "y": "Number of sequences"},
        title=title,
        text=entries,
        template = 'simple_white',
        color_discrete_sequence=px.colors.qualitative.Set2

    )

    fig.update_traces(
        textfont_size=12,
        textangle=0,
        textposition="outside",
        cliponaxis=False,
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
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )
    return plot_div


def plot_count_distribution_bar_chart(df: pd.DataFrame):
    count_distribution = df["Count"].value_counts().reset_index()
    count_distribution.columns = ["Number of occurrences", "Number of genes"]

    count_distribution = count_distribution.sort_values(by="Number of occurrences", ascending=False)

    fig = px.bar(
        count_distribution,
        x="Number of occurrences",
        y="Number of genes",
        labels={"Number of occurrences": "Number of occurrences", "Number of genes": "Number of genes"},
        title="Distribution of Gene Occurrences",
        color_discrete_sequence=px.colors.qualitative.Set2
    )

    fig.update_traces(
        text=count_distribution["Number of genes"],  # Zet de waarden als tekst in de balken
        textposition="outside",
        hoverinfo="y",
        hovertemplate="<b>Number of occurrences:</b> %{x}<br><b>Number of genes:</b> %{y}<extra></extra>"
    )

    fig.update_layout(
        xaxis_title="Number of occurrences",
        yaxis_title="Number of genes",
        xaxis=dict(tickangle=-45, showgrid=False),
        yaxis=dict(showgrid=True),
        dragmode=False,
        template="plotly_white"
    )

    plot_div = plot(
        fig,
        config={
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d',
                                       'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )

    return plot_div


def get_pi_plot(df_filtered: pd.DataFrame, groupby_status: str, title: str) -> dict:
    segment_counts = df_filtered.groupby(groupby_status).size().reset_index(name="Count")

    fig = px.pie(
        segment_counts,
        names=groupby_status,
        values="Count",
        color=groupby_status,
        color_discrete_sequence=px.colors.qualitative.Set2,
        hole=0.4
    )

    fig.update_traces(
        textinfo="label+value",
        hoverinfo="label+value+percent"
    )

    fig.update_layout(
        title=dict(
            text=f"{title}",
            font=dict(size=18, family="Arial", color="black"),
            x=0.5
        ),
        template="plotly_white"
    )

    plot_div = plot(
        fig,
        config={
            'displaylogo': False,
            'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'resetScale2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d'],
            'modeBarButtonsToAdd': ['toImage'],
            'scrollZoom': False,
            'toImageButtonOptions': {
                'format': 'svg',
                'filename': 'custom_image',
                'scale': 1
            }
        },
        output_type='div',
        include_plotlyjs=True
    )

    return plot_div


def get_dna_vieuwer_plot(df_filtered: pd.DataFrame) -> dict:
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
        bokeh_plot = record.plot_with_bokeh(figure_width=20, figure_height=5)
        bokeh_plot.title.text_font_size = "20pt"
        bokeh_plot.title.align = "center"

        bokeh_plot.xaxis.formatter = NumeralTickFormatter(format="0")
        graph = file_html(bokeh_plot, CDN, "Region Viewer")
        return graph


def get_venn_diagram(set1, set2, label1, label2, status):
    fig, ax = plt.subplots()
    venn2([set1, set2], set_labels=(label1, label2),  set_colors=("lightblue", "salmon"))
    plt.title(f"{status}")
    img = io.BytesIO()
    plt.savefig(img, format='svg')
    img.seek(0)
    plt.close(fig)

    encoded_svg = base64.b64encode(img.getvalue()).decode('utf-8')
    return f"data:image/svg+xml;base64,{encoded_svg}"
