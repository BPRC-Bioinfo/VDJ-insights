import argparse
from collections import deque
from itertools import zip_longest
import itertools
from pathlib import Path
from typing import List, Tuple, Dict, Any
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from collections import Counter
from reevaluate import re_evaluate_main
from logger import custom_logger


# Method for logger current states of the program.
logger = custom_logger(__name__)

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['svg.fonttype'] = 'none'

class PlotGenerator:
    """
    A class to generate plots for genomic segments.

    Attributes:
        colors (dict): A dictionary mapping segment types to color codes.
        legend_labels (dict): A dictionary mapping segment types to labels.
        style (str): The plotting style. Currently supports only "combined".
        footer (bool): Flag to add footer to plots.
        cwd (Path): Current working directory.
        output_dir_name (str): Name of the output directory.
        output_dir (Path): Path to the output directory.
    """

    def __init__(self) -> None:
        """Initializes the PlotGenerator with default settings."""
        self.legend_labels = {
            "V": "V", "D": "D", "J": "J","C":"C", "Known_a1": "Known allele","Known_a2": "Known allele", "Novel_a1": "Novel allele","Novel_a2": "Novel allele"
                                                                                                                      }
        self.cwd = Path.cwd()

    def add_filter(self, elements: List[str]) -> List[str]:
        """Filters elements by removing characters after '*' and '-'."""
        return [elem.split('*')[0].split("-")[0] for elem in elements]

    def align_lists_by_matching(self, list1: List[str], list2: List[str]) -> Tuple[List[str], List[str]]:
        """Aligns two lists by matching elements, inserting blanks for mismatches."""
        aligned_list1, aligned_list2 = [], []
        i, j = 0, 0
        while i < len(list1) or j < len(list2):
            if i < len(list1) and j < len(list2):
                if list1[i] == list2[j]:
                    aligned_list1.append(list1[i])
                    aligned_list2.append(list2[j])
                    i += 1
                    j += 1
                else:
                    if list1[i] in list2[j:]:
                        aligned_list1.append("")
                        aligned_list2.append(list2[j])
                        j += 1
                    elif list2[j] in list1[i:]:
                        aligned_list1.append(list1[i])
                        aligned_list2.append("")
                        i += 1
                    else:
                        aligned_list1.append(list1[i])
                        aligned_list2.append(list2[j])
                        i += 1
                        j += 1
            elif i < len(list1):
                aligned_list1.append(list1[i])
                aligned_list2.append("")
                i += 1
            elif j < len(list2):
                aligned_list1.append("")
                aligned_list2.append(list2[j])
                j += 1
        return aligned_list1, aligned_list2

    def update_df(self, df: pd.DataFrame, haplotype_list: List[str]) -> pd.DataFrame:
        """Updates a DataFrame based on the provided haplotype list."""
        new_df = pd.DataFrame(columns=df.columns)
        row_counter = 0
        for value in haplotype_list:
            if value == "":
                new_df = pd.concat([new_df, pd.DataFrame(
                    {col: "" for col in df.columns}, index=[0])], ignore_index=True)
            else:
                if row_counter < len(df):
                    new_df = pd.concat([new_df, pd.DataFrame(
                        df.iloc[row_counter]).T], ignore_index=True)
                    row_counter += 1
                else:
                    print(
                        f"Warning: Attempted to access row {row_counter} but only {len(df)} rows are available.")
        return new_df

    def add_to_dict(self, d, total_dict):
        for key, value in d.items():
            total_dict.setdefault(key, list()).append(value)
        return total_dict

    def generate_footer_lists(self, aligned_hap1, aligned_hap2, df_hap1, df_hap2):
        main, secondary = list(), list()
        dict_hap1 = df_hap1.set_index('Start coord')['Short name'].to_dict()
        dict_hap2 = df_hap2.set_index('Start coord')['Short name'].to_dict()
        for hap1, hap2 in zip_longest(aligned_hap1, aligned_hap2):
            hap1_value, hap2_value = dict_hap1.get(hap1), dict_hap2.get(hap2)
            if hap1_value == hap2_value:
                main.append(hap1_value)
                secondary.append("")
            elif hap1_value != hap2_value:
                main.append(hap2_value)
                secondary.append(hap1_value)

        return secondary, main

    def run_alignment(self, df_hap1: pd.DataFrame, df_hap2: pd.DataFrame, region: str) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
        """Runs alignment on DataFrame based on the region and returns aligned DataFrames and merged haplotypes."""
        hap1, hap2 = [i["Start coord"].to_list() for i in [df_hap1, df_hap2]]
        aligned_hap1, aligned_hap2 = self.align_lists_by_matching(hap1, hap2)
        footer = self.generate_footer_lists(
            aligned_hap1, aligned_hap2, df_hap1, df_hap2)
        df_alignment_hap1, df_alignment_hap2 = self.update_df(
            df_hap1, aligned_hap1), self.update_df(df_hap2, aligned_hap2)
        return df_alignment_hap1, df_alignment_hap2, footer

    def load_dataframe(self, filepath: str) -> pd.DataFrame:
        """Loads a DataFrame from an Excel file located at the given filepath."""
        read_path = Path(self.cwd / filepath)
        return pd.read_excel(read_path)

    def generate_product(self, regions: List[str], haplotypes: List[str]) -> List[Tuple[str, str]]:
        """Generates a Cartesian product of regions and haplotypes."""
        return list(itertools.product(regions, haplotypes))

    def calculate_plot_size(self, num_segments: int) -> Tuple[int, float, float, int]:
        """Calculates and returns plot size parameters based on the number of segments."""
        if num_segments < 20:
            return 1, 1.5, 0.2, 10  # block size, block spacing, width per block, font size
        elif num_segments < 40:
            return 8, 3, 0.2, 8
        else:
            return 8, 1, 0.2, 6

    def create_dual_color_features(self, ax: plt.Axes, start_pos: float, block_size: float,
                                   segment_color: str, fcolor: str) -> None:
        """Creates a dual-colored block feature on the plot."""
        half_block_size = block_size / 10
        bottom_block_y = -half_block_size / 2
        top_block_y = half_block_size / 2
        ax.add_patch(patches.Rectangle((start_pos, bottom_block_y), block_size,
                     half_block_size, facecolor=fcolor, edgecolor='black', linewidth=0))
        ax.add_patch(patches.Rectangle((start_pos, top_block_y), block_size,
                     half_block_size, facecolor=segment_color, edgecolor='black', linewidth=0))

    def calculate_total_width(self, df: pd.DataFrame, block_size: float, block_spacing: float) -> float:
        """Calculates the total width of the plot based on DataFrame size, block size, and spacing."""
        return len(df) * (block_size + block_spacing) - block_spacing

    def initialize_plot(self, total_width: float, block_size: float) -> Tuple[plt.Figure, plt.Axes]:
        """Initializes and returns a matplotlib plot with the specified dimensions."""
        fig, ax = plt.subplots(figsize=(20, 3))
        y_padding = 2
        ax.set_xlim(0, total_width)
        ax.set_ylim(-block_size, block_size + y_padding)
        ax.axis('off')
        ax.grid(False)
        ax.figure.set_tight_layout(True)
        return fig, ax

    def add_segment_blocks(self, ax: plt.Axes, df: pd.DataFrame, block_size: float, block_spacing: float) -> None:
        """Adds segment blocks to the plot based on DataFrame data."""
        current_position = 0
        for _, row in df.iterrows():
            self.create_dual_color_features(
                ax, current_position, block_size, self.colors[row['Segment']], self.colors[row['Allele']])
        #        ax, current_position, block_size, self.colors[row['Segment']], self.colors[row['Status']])
            current_position += block_size + block_spacing


    def segment_legend(self, ax: plt.Axes) -> None:
        """Adds a segment legend to the plot."""
        color_legend_handles = [patches.Patch(
            color=self.colors[key], label=label) for key, label in self.legend_labels.items()]
        color_legend = ax.legend(handles=color_legend_handles, loc='best', bbox_to_anchor=(
            0.5, 1), ncol=len(color_legend_handles), handletextpad=1)
        ax.add_artist(color_legend)

    def sample_haplotype_legend(self, ax: plt.Axes, region: str, haplotype: str, sample: str) -> plt.Axes:
        """Adds a sample haplotype legend to the plot."""
        region_haplotype_legend = ax.legend(handles=[patches.Patch(
            color='none', label=f"{region}-{haplotype} ({sample})")], loc='best', bbox_to_anchor=(1, 1), fontsize=10)
        region_haplotype_legend.get_title().set_fontsize('12')
        return ax

    def add_footer(self, ax: plt.Axes, merged: List[str], step: float, block_size: float, block_spacing: float, haplotype: str) -> None:
        """Adds a footer to the plot if specified in the style."""
        total_segments = len(merged)
        start_pos = 0
        for index, item in enumerate(merged):
            x_position = start_pos + block_size / 2
            y_position = -1 if total_segments < 20 else -6
            font_size = 8 if total_segments < 20 else 5
            if item:
                ax.text(x_position, y_position, item, ha='center',
                #        va='bottom', fontsize=font_size, rotation=90)
                        va='top', fontsize=font_size, rotation=90)
            start_pos += block_size + block_spacing

    def add_legends(self, ax: plt.Axes, region: str, haplotype: str, seg_type: str, sample: str) -> None:
        """Adds legends to the plot."""
        if self.style == "combined" and haplotype == "hap1" or self.style == "single":
            self.segment_legend(ax)
        self.sample_haplotype_legend(ax, region, haplotype, sample)

    def add_region_indicators(self, ax: plt.Axes, region: str, region_indicators: Dict[str, deque], block_size: float, block_spacing: float) -> None:
        """Adds region indicators to the plot."""
        for key, deq in region_indicators.items():
            if deq:
                x_min = min(deq) * (block_size + block_spacing)
                x_max = max(deq) * (block_size + block_spacing) + block_size
                single_color = self.colors[key]
                y_position_line = 2 if block_size != 1 else 0.5
                ax.hlines(y_position_line, x_min, x_max,
                          colors=single_color, linestyles='solid', linewidth=5)
                x_mid = (x_min + x_max) / 2
                text_y_position = 3 if block_size != 1 else 0.7
                ax.text(x_mid, text_y_position, f'{region}{key}', horizontalalignment='center', verticalalignment='bottom', fontsize=10, color='black', bbox=dict(
                    facecolor='white', edgecolor='black', boxstyle='round,pad=0.3', linewidth=0.2))

    def get_index(self, segments: List[str]) -> Dict[str, deque]:
        """Generates indices for segment blocks."""
        index_counter = {"V": deque(), "D": deque(), "J": deque(),"C": deque()}
        counter = Counter([item for item in segments[0:3] if item != ""])

        most_common_element, count = counter.most_common(1)[0]

        collect = deque(most_common_element)
        for x, seg in enumerate(segments):
            if seg not in collect and seg != '':
                seq_check = segments[x+1] if x+1 < len(segments) else None
                if seq_check not in collect:
                    collect.append(seg)
            index_counter[collect[-1]].append(x)
        return index_counter

    def configure_plot(self, df: pd.DataFrame, block_size: float, block_spacing: float,
                       block_width: float, region: str, haplotype: str, seg_type: str,
                       sample: str, merged: List[str]) -> Tuple[plt.Figure, plt.Axes]:
        """
        Configures the plot with segments, legends, and region indicators.

        Parameters:
            df: DataFrame containing the segments and their statuses.
            block_size: The size of each block in the plot.
            block_spacing: The spacing between blocks in the plot.
            block_width: The width per block in the plot.
            region: The genomic region being plotted.
            haplotype: The haplotype identifier.
            seg_type: The type of segment.
            sample: The sample identifier.
            merged: The list of merged haplotypes.

        Returns:
            A tuple containing the Figure and Axes objects for the plot.
        """
        all_list = list(df["Segment"])
        region_indicators = self.get_index(all_list)
        total_width = self.calculate_total_width(df, block_size, block_spacing)
        fig, ax = self.initialize_plot(total_width, block_size)
        self.add_segment_blocks(ax, df, block_size, block_spacing)
        self.add_legends(ax, region, haplotype, seg_type, sample)
        self.add_region_indicators(
            ax, region, region_indicators, block_size, block_spacing)
        if self.footer:
            self.add_footer(ax, merged, block_size + block_spacing,
                            block_size, block_spacing, haplotype)
        return fig, ax

    def save_plot(self, fig: plt.Figure, output_dir: Path, region: str, haplotype: str) -> None:
        """
        Saves the generated plot to the specified directory.

        Parameters:
            fig: The Figure object of the plot to save.
            output_dir: The directory path where the plot will be saved.
            region: The genomic region of the plot.
            haplotype: The haplotype identifier.
        """
        filename = f"{region}-{haplotype}.svg"
        filepath = output_dir / filename
        logger.info(f"Generation plot {filename}, saving to {filepath}")
        #fig.savefig(filepath, dpi=300, bbox_inches='tight', format='svg')
        fig.savefig(filepath, bbox_inches='tight', format='svg')

    def create_plot(self, df: pd.DataFrame, output_dir: Path, region: str,
                    haplotype: str, sample: str, merged: List[str]) -> Tuple[plt.Figure, plt.Axes]:
        """
        Creates a plot for the given DataFrame and parameters.

        Parameters:
            df: DataFrame containing the segments and their statuses.
            output_dir: The directory path where the plot will be saved.
            region: The genomic region of the plot.
            haplotype: The haplotype identifier.
            sample: The sample identifier.
            merged: The list of merged haplotypes.

        Returns:
            A tuple containing the Figure and Axes objects for the plot.
        """
        if not df.empty:
            block_size, block_spacing, base_width_per_block, font_size = self.calculate_plot_size(
                len(df))
            return self.configure_plot(df, block_size, block_spacing, base_width_per_block, region, haplotype, 'test', sample, merged)


    def load_and_combine_dataframes(self, files: List[str]) -> pd.DataFrame:
        """
        Loads multiple Excel files as DataFrames, assigns a status to each, and combines them into a single DataFrame.

        Parameters:
            files: A list of tuples, each containing the file path and the status to assign to its DataFrame.

        Returns:
            A combined DataFrame containing all data from the files with assigned statuses.
        """
        dataframes = [self.load_dataframe(filename) for filename in files]
        combined_df = pd.concat(dataframes, ignore_index=True)

        def determine_alleles(row):
            # Check if Start coord and End coord appear double
            if combined_df[(combined_df['Start coord'] == row['Start coord']) & (combined_df['stop'] == row['stop'])].shape[0] > 1:
                # Find rows with same Start coord and End coord
                same_coords = combined_df[(combined_df['Start coord'] == row['Start coord']) & (combined_df['stop'] == row['stop'])]
                # Check if Old name-like seq is identical
                if same_coords['sequence'].nunique() == 1:
                    if row['Status'] == 'Known':
                        # Alleles become Known_a1
                        return 'Known_a1'
                    else:
                        return 'Novel_a1'
                else:
                    if same_coords['Status'].nunique() > 1:
                        # Check if Status is known or novel
                        if row['Status'] == 'Known':
                            return 'Known_a1'
                        else:
                            return 'Novel_a1'
                    else:
                        if row['Haplotype'] == 'hap1':
                            if row['Status'] == 'Known':
                                return 'Known_a1'
                            else:
                                return 'Novel_a1'
                        else:
                            if row['Status'] == 'Known':
                                return 'Known_a2'
                            else:
                                return 'Novel_a2'
            else:
                if row['Status'] == 'Known':
                    return 'Known_a1'
                else:
                    return 'Novel_a1'
        combined_df['Allele'] = combined_df.apply(determine_alleles, axis=1)

        return combined_df

    def combine_and_save_plots_vertically(self, fig1: plt.Figure, fig2: plt.Figure, output_path: Path, dpi: int = 300) -> None:
        """
        Combines two matplotlib figures vertically and saves the combined figure to the specified path.

        Parameters:
            fig1: The first Figure object.
            fig2: The second Figure object.
            output_path: The path where the combined figure will be saved.
            dpi: The resolution of the output image.
        """
        def fig_to_numpy(fig: plt.Figure) -> np.ndarray:
            """Converts a Figure to a NumPy array."""
            canvas = FigureCanvas(fig)
            canvas.draw()
            buf = canvas.buffer_rgba()
            ncols, nrows = canvas.get_width_height()
            return np.frombuffer(buf, dtype=np.uint8).reshape(nrows, ncols, 4)

        arr1, arr2 = fig_to_numpy(fig1), fig_to_numpy(fig2)
        max_width = max(arr1.shape[1], arr2.shape[1])
        combined_arr = np.vstack([np.pad(arr, ((0, 0), (0, max_width - arr.shape[1]),
                                               (0, 0)), 'constant', constant_values=255) for arr in [arr1, arr2]])
        combined_fig, combined_ax = plt.subplots(
            figsize=(max_width / dpi, combined_arr.shape[0] / dpi), dpi=dpi)
        combined_ax.imshow(combined_arr, aspect='auto')
        combined_ax.axis('off')
        plt.subplots_adjust(left=0, right=1, top=1,
                            bottom=0, wspace=0, hspace=0)
        logger.info(
            f"Creating {output_path.name} and saving it to {self.cwd /output_path}!")
        plt.savefig(self.cwd / output_path, dpi=dpi, bbox_inches='tight',
                    pad_inches=0, transparent=True)
        plt.close(combined_fig)

    def orientation_calculation(self, df):
        segment = df.sort_values(by="Start coord").head(n=3)
        single_segment = segment["Segment"].mode().iloc[0]
        return False if single_segment in ("J", "D") else True

    def run_seperate_haplotype(self, df: pd.DataFrame, region, sample):
        for hap in ["hap1", "hap2"]:
            hap_df = df.query(f"Haplotype == '{hap}'")
            if not hap_df.empty:
                orientation_check = self.orientation_calculation(
                    hap_df) if self.rotate else True
                hap_df = hap_df.sort_values(
                    by="Start coord", ascending=orientation_check)
                merged = hap_df["Filtered"]
                plot = self.create_plot(hap_df, self.output_dir,
                                        region, hap, sample, merged)
                self.save_plot(plot[0], self.output_dir, region, hap)

    def get_haplotype(self, df):
        return df["Haplotype"].mode()[0]

    def run_combined_haplotype(self, new_df: pd.DataFrame, region, sample):
        orientation_check = self.orientation_calculation(
            new_df) if self.rotate else True
        df_haplotypes = [new_df.query(f"Haplotype == '{i}'") for i in [
            "hap1", "hap2"]]
        df_hap1, df_hap2 = [hap.sort_values(
            by="Start coord", ascending=orientation_check) for hap in df_haplotypes]
        if not df_hap1.empty and not df_hap2.empty:
            hap1, hap2, merged = self.run_alignment(df_hap1, df_hap2, region)
            plots = [self.create_plot(
                hap, self.output_dir, region, self.get_haplotype(hap),
                sample, footer) for hap, footer in zip([hap1, hap2], merged)]
            self.combine_and_save_plots_vertically(
                plots[0][0], plots[1][0], self.output_dir / f"{region}.svg")

    def process_region(self, df, references):
        """Process references within a region, aligning haplotypes if conditions are met."""
        if len(references) > 2:
            return df[df['Path'].isin(references)]
        return df

    def process_region_data(self, df: pd.DataFrame, region: str, sample: str) -> None:
        """Process data for a specific region."""
        region_df = df.query(f"Region == '{region}'")
        references = region_df['Path'].value_counts()
        needed_references = references.head(2).index.tolist()
        new_df = self.process_region(region_df, needed_references)
        print(len(new_df))
        filtered = new_df.copy()
        filtered['Filtered'] = filtered['Short name'].apply(
            lambda x: self.add_filter([x])[0])
        if self.style == "single":
            self.run_seperate_haplotype(filtered, region, sample)
        else:
            self.run_combined_haplotype(filtered, region, sample)


def parse_arguments():
    cwd = Path.cwd()
    parser = argparse.ArgumentParser(
        description="Generate plots for VDJ segments.")
    # Setting default file paths and statuses
    default_files = [
        "annotation_report_100%_plus.xlsx",
        "annotation_report_plus.xlsx"
    ]
    parser.add_argument('-f', '--files', nargs='+', default=default_files,
                        help="List of file paths to process.")
    parser.add_argument('-o', '--output_dir', type=str,
                        help="Path to the output directory.", default=cwd / "VDJ_display")
    parser.add_argument('-s', '--style', choices=['single', 'combined'],
                        help="The plotting style: 'single' for separate haplotypes, 'combined' for combined haplotypes.", default='combined')
    parser.add_argument('--no-footer', action='store_false',
                        help="Flag to not add footer to plots.", dest='footer')
    parser.add_argument('-r', '--rotate', action='store_true', default=False,
                        help="Set orientation to True, defaults to False if not provided.")
    parser.add_argument('-c', '--color_theme', choices=[
        'default', 'easy_eye', 'IBM', 'Wong', 'Tol','Susan'
    ], help="Choose a color theme for the plots.", default='default')
    parser.add_argument('--re-evaluate', action='store_true', default=False,
                        help="Reevaluate the VDJ segments to better determine the order")

    return parser.parse_args()


def main():
    # Define color themes
    color_themes = {
        'default': {
            "V": '#ffdfba',
            "D": '#baffc9',
            "J": '#bae1ff',
            "C": '#F3F29B',
            "Known_a1": '#8287D4',
            "Known_a2": '#CFD2FB',
            "Novel_a1": '#FF8D98',
            "Novel_a2": '#FFDDE0',
            "": '#ffffff'
        },
        'easy_eye': {
            'V': '#E69F00',
            'D': '#56B4E9',
            'J': '#009E73',
            'C': '#008F00',
            'Known_a1': '#DED227',
            'Known_a2': '#F0E989',
            'Novel_a1': '#0072B2',
            'Novel_a2': '#77A8C3',
            '': '#FFFFFF'
        },
        'IBM': {
            'V': '#785EF0',
            'D': '#DC267F',
            'J': '#FE6100',
            'C': '#F0E442',
            'Known_a1': '#FFB000',
            'Known_a2': '#E69F00',
            'Novel_a1': '#648FFF',
            'Novel_a2': '#56B4E9',
            '': '#FFFFFF'
        },
        'Wong': {
            'V': '#CC79A7',
            'D': '#009E73',
            'J': '#D55E00',
            'C': '#000000',
            'Known_a1': '#E69F00',
            'Known_a2': '#F0E442',
            'Novel_a1': '#0072B2',
            'Novel_a2': '#56B4E9',
            '': '#FFFFFF'
        },
        'Tol': {
            'V': '#88CCEE',
            'D': '#44AA99',
            'J': '#117733',
            'C': '#332288',
            'Known_a1': '#882255',
            'Known_a2': '#AA4499',
            'Novel_a1': '#CC6677',
            'Novel_a2': '#DDCC77',
            '': '#FFFFFF'
        },
        'Susan': {
            'V': '#BDB2CE',
            'D': '#C8A588',
            'J': '#C9E79F',
            'C': '#B2D4E0',
            'Known_a1': '#36753B',
            'Known_a2': '#61A899',
            'Novel_a1': '#C66526',
            'Novel_a2': '#DCA237',
            '': '#FFFFFF'
        },

        # Add more themes as needed
    }
    args = parse_arguments()
    input_files = [file.rstrip(' ,;.') for file in args.files]
    plot_generator = PlotGenerator()
    # Apply selected color theme
    plot_generator.colors = color_themes[args.color_theme]
    plot_generator.output_dir = Path(args.output_dir)
    plot_generator.style = args.style
    plot_generator.footer = args.footer
    plot_generator.rotate = args.rotate
    reevaluated_excel = Path(Path.cwd() / "reevaluated.xlsx")
    if args.re_evaluate:
        logger.info(
            f"Reevaluating {' '.join(input_files)}!"
            f"See {reevaluated_excel} for the results!")
        if not reevaluated_excel.exists():
            input_files = [re_evaluate_main(input_files)]
        else:
            input_files = [reevaluated_excel]

    plot_generator.output_dir.mkdir(parents=True, exist_ok=True)

    df = plot_generator.load_and_combine_dataframes(input_files)
    regions = list(df["Region"].unique())
    sample = df.iloc[0]["Sample"]
    for region in regions:
        plot_generator.process_region_data(df, region, sample)


if __name__ == "__main__":
    main()
