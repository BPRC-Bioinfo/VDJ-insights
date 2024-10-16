import argparse
from collections import deque
from itertools import zip_longest
import itertools
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from collections import Counter
from reevaluate import re_evaluate_main
from logger import console_logger, file_logger

logger = console_logger(__name__)
file_log = file_logger(__name__)


plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


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

    def __init__(self):
        """
        Initializes the PlotGenerator with default settings. Sets the current working directory (cwd).
        Initializes default attributes for color themes and legend labels. Other attributes are defined
        during runtime, such as the output directory, style, and footer flag.
        """
        self.cwd = Path.cwd()

    def add_filter(self, elements):
        """
        Filters elements by removing characters after '*' and '-'. Strips extra identifiers
        from the input list of strings to retain the core segment names.

        Args:
            elements (list): List of segment strings to filter.

        Returns:
            filtered_elements (list): Filtered list of segment names.
        """
        return [elem.split('*')[0].split("-")[0] for elem in elements]

    def align_lists_by_matching(self, list1, list2):
        """
        Aligns two lists by matching elements, inserting blanks for mismatches. This function compares two lists,
        aligns matching elements, and inserts empty strings ("") where elements do not match.

        Args:
            list1 (list): First list to align.
            list2 (list): Second list to align.

        Returns:
            aligned_list1 (list): Aligned versions of list 1.
            aligned_list2 (list): Aligned versions of list 2.
        """
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

    def update_df(self, df, haplotype_list):
        """
        Updates a DataFrame by aligning it to a given haplotype list. This function updates
        rows in the DataFrame based on the elements in the haplotype list, aligning rows 
        or adding blank entries where needed.

        Args:
            df (pd.DataFrame): The DataFrame to update.
            haplotype_list (list): List of haplotypes used for alignment.

        Returns:
            new_df (pd.DataFrame): Updated DataFrame aligned with the haplotype list.
        """
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
                    file_log.warning(
                        f"Warning: Attempted to access row {row_counter} but only {len(df)} rows are available.")
        return new_df

    def add_to_dict(self, d, total_dict):
        """
        Adds key-value pairs from one dictionary to another, appending values to lists.

        Args:
            d (dict): The source dictionary.
            total_dict (dict): The target dictionary where key-value pairs are appended.

        Returns:
            total_dict (dict): Updated total dictionary with appended values.
        """
        for key, value in d.items():
            total_dict.setdefault(key, list()).append(value)
        return total_dict

    def generate_footer_lists(self, aligned_hap1, aligned_hap2, df_hap1, df_hap2):
        """
        Generates the footer lists for aligned haplotypes. Compares the two haplotypes and generates
        two lists: one for the main haplotype and another for secondary elements where haplotypes differ.

        Args:
            aligned_hap1 (list): Aligned list for haplotype 1.
            aligned_hap2 (list): Aligned list for haplotype 2.
            df_hap1 (pd.DataFrame): DataFrame for haplotype 1.
            df_hap2 (pd.DataFrame): DataFrame for haplotype 2.

        Returns:
            secondary (list): List of secondary elements from haplotype 1.
            main (list): List of main elements from haplotype 2.
        """
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

    def run_alignment(self, df_hap1, df_hap2, region):
        """
        Runs alignment on DataFrame based on the region and returns aligned DataFrames and merged haplotypes. 
        Aligns haplotypes using their starting coordinates and generates aligned DataFrames for each haplotype.

        Args:
            df_hap1 (pd.DataFrame): DataFrame for haplotype 1.
            df_hap2 (pd.DataFrame): DataFrame for haplotype 2.
            region (str): The genomic region being processed.

        Returns:
            df_alignment_hap1 (pd.DataFrame): Aligned DataFrame for haplotype 1.
            df_alignment_hap2 (pd.DataFrame): Aligned DataFrame for haplotype 2.
            footer (list): List of merged haplotypes.
        """
        hap1, hap2 = [i["Start coord"].to_list() for i in [df_hap1, df_hap2]]
        aligned_hap1, aligned_hap2 = self.align_lists_by_matching(hap1, hap2)
        footer = self.generate_footer_lists(
            aligned_hap1, aligned_hap2, df_hap1, df_hap2)
        df_alignment_hap1, df_alignment_hap2 = self.update_df(
            df_hap1, aligned_hap1), self.update_df(df_hap2, aligned_hap2)
        return df_alignment_hap1, df_alignment_hap2, footer

    def load_dataframe(self, filepath):
        """
        Loads a DataFrame from an Excel file located at the given filepath.

        Args:
            filepath (str): Path to the Excel file.

        Returns:
            df (pd.DataFrame): Loaded DataFrame.
        """
        read_path = Path(self.cwd / filepath)
        return pd.read_excel(read_path)

    def generate_product(self, regions, haplotypes):
        """
        Generates a Cartesian product of regions and haplotypes.

        Args:
            regions (list): List of regions.
            haplotypes (list): List of haplotypes.

        Returns:
            product (list): Cartesian product of regions and haplotypes as tuples.
        """
        return list(itertools.product(regions, haplotypes))

    def calculate_plot_size(self, num_segments):
        """
        Calculates and returns plot size parameters based on the number of segments.

        Args:
            num_segments (int): Number of segments to plot.

        Returns:
            block_size (int): Block size for the plot.
            block_spacing (float): Spacing between blocks.
            width_per_block (float): Width of each block.
            font_size (int): Font size for the plot.
        """
        if num_segments < 20:
            return 1, 1.5, 0.2, 10
        elif num_segments < 40:
            return 8, 3, 0.2, 8
        else:
            return 8, 1, 0.2, 6

    def create_dual_color_features(self, ax, start_pos, block_size, segment_color, fcolor):
        """
        Creates a dual-colored block feature on the plot.

        Args:
            ax (plt.Axes): The Axes object of the plot.
            start_pos (float): Starting position of the block.
            block_size (float): Size of the block.
            segment_color (str): Color for the upper half of the block.
            fcolor (str): Color for the lower half of the block.
        """
        half_block_size = block_size / 10
        bottom_block_y = -half_block_size / 2
        top_block_y = half_block_size / 2
        ax.add_patch(patches.Rectangle((start_pos, bottom_block_y), block_size,
                     half_block_size, facecolor=fcolor, edgecolor='black', linewidth=0))
        ax.add_patch(patches.Rectangle((start_pos, top_block_y), block_size,
                     half_block_size, facecolor=segment_color, edgecolor='black', linewidth=0))

    def calculate_total_width(self, df, block_size, block_spacing):
        """
        Calculates the total width of the plot based on DataFrame size, block size, and spacing.

        Args:
            df (pd.DataFrame): DataFrame containing the segments.
            block_size (float): Size of each block in the plot.
            block_spacing (float): Spacing between blocks.

        Returns:
            total_width (float): Total width of the plot.
        """
        return len(df) * (block_size + block_spacing) - block_spacing

    def initialize_plot(self, total_width, block_size):
        """
        Initializes and returns a matplotlib plot with the specified dimensions.

        Args:
            total_width (float): Total width of the plot.
            block_size (float): Size of each block in the plot.

        Returns:
            fig (plt.Figure): The Figure object for the plot.
            ax (plt.Axes): The Axes object for the plot.
        """
        fig, ax = plt.subplots(figsize=(20, 3))
        y_padding = 2
        ax.set_xlim(0, total_width)
        ax.set_ylim(-block_size, block_size + y_padding)
        ax.axis('off')
        ax.grid(False)
        ax.figure.set_tight_layout(True)
        return fig, ax

    def add_segment_blocks(self, ax, df, block_size, block_spacing):
        """
        Adds segment blocks to the plot based on DataFrame data.

        Args:
            ax (plt.Axes): The Axes object of the plot.
            df (pd.DataFrame): DataFrame containing the segments and their statuses.
            block_size (float): Size of each block.
            block_spacing (float): Spacing between blocks.
        """
        current_position = 0
        for _, row in df.iterrows():
            self.create_dual_color_features(
                ax, current_position, block_size, self.colors[row['Segment']], self.colors[row['Status']])
            current_position += block_size + block_spacing

    def segment_legend(self, ax):
        """
        Adds a segment legend to the plot.

        Args:
            ax (plt.Axes): The Axes object of the plot.
        """
        color_legend_handles = [patches.Patch(
            color=self.colors[key], label=label) for key, label in self.legend_labels.items()]
        color_legend = ax.legend(handles=color_legend_handles, loc='best', bbox_to_anchor=(
            0.5, 1), ncol=len(color_legend_handles), handletextpad=1)
        ax.add_artist(color_legend)

    def sample_haplotype_legend(self, ax, region, haplotype, sample):
        """
        Adds a sample haplotype legend to the plot.

        Args:
            ax (plt.Axes): The Axes object of the plot.
            region (str): The genomic region of the plot.
            haplotype (str): The haplotype identifier.
            sample (str): The sample identifier.

        Returns:
            ax (plt.Axes): The updated Axes object with the added legend.
        """
        region_haplotype_legend = ax.legend(handles=[patches.Patch(
            color='none', label=f"{region}-{haplotype} ({sample})")], loc='best', bbox_to_anchor=(1, 1), fontsize=10)
        region_haplotype_legend.get_title().set_fontsize('12')
        return ax

    def add_footer(self, ax, merged, step, block_size, block_spacing, haplotype):
        """
        Adds a footer to the plot if specified in the style.

        Args:
            ax (plt.Axes): The Axes object of the plot.
            merged (list): List of merged haplotypes.
            step (float): Step value for the footer.
            block_size (float): Size of each block.
            block_spacing (float): Spacing between blocks.
            haplotype (str): The haplotype identifier.
        """
        total_segments = len(merged)
        start_pos = 0
        for index, item in enumerate(merged):
            x_position = start_pos + block_size / 2
            y_position = -1 if total_segments < 20 else -6
            font_size = 8 if total_segments < 20 else 5
            if item:
                ax.text(x_position, y_position, item, ha='center',
                        va='bottom', fontsize=font_size, rotation=90)
            start_pos += block_size + block_spacing

    def add_legends(self, ax, region, haplotype, seg_type, sample):
        """
        Adds legends to the plot.

        Args:
            ax (plt.Axes): The Axes object of the plot.
            region (str): The genomic region of the plot.
            haplotype (str): The haplotype identifier.
            seg_type (str): The segment type.
            sample (str): The sample identifier.
        """
        if self.style == "combined" and haplotype == "hap1" or self.style == "single":
            self.segment_legend(ax)
        self.sample_haplotype_legend(ax, region, haplotype, sample)

    def add_region_indicators(self, ax, region, region_indicators, block_size, block_spacing):
        """
        Adds region indicators to the plot.

        Args:
            ax (plt.Axes): The Axes object of the plot.
            region (str): The genomic region of the plot.
            region_indicators (dict): Dictionary mapping region indicators to positions.
            block_size (float): Size of each block.
            block_spacing (float): Spacing between blocks.
        """
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

    def get_index(self, segments):
        """
        Generates indices for segment blocks.

        Args:
            segments (list): List of segments in the plot.

        Returns:
            index_counter (dict): Dictionary containing the segment block indices.
        """
        index_counter = {"V": deque(), "D": deque(), "J": deque()}
        counter = Counter([item for item in segments[0:3] if item != ""])

        most_common_element, count = counter.most_common(1)[0]

        collect = deque(most_common_element)
        for x, seg in enumerate(segments):
            if seg not in collect and seg != '' and segments[x+1] not in collect:
                collect.append(seg)
            index_counter[collect[-1]].append(x)
        return index_counter

    def configure_plot(self, df, block_size, block_spacing,
                       block_width, region, haplotype, seg_type,
                       sample, merged):
        """
        Configures the plot with segments, legends, and region indicators.

        Args:
            df (pd.DataFrame): DataFrame containing the segments and their statuses.
            block_size (float): The size of each block in the plot.
            block_spacing (float): The spacing between blocks in the plot.
            block_width (float): The width per block in the plot.
            region (str): The genomic region being plotted.
            haplotype (str): The haplotype identifier.
            seg_type (str): The type of segment.
            sample (str): The sample identifier.
            merged (list): The list of merged haplotypes.

        Returns:
            fig (plt.Figure): The Figure object for the plot.
            ax (plt.Axes): The Axes object for the plot.
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

    def save_plot(self, fig, output_dir, region, haplotype):
        """
        Saves the generated plot to the specified directory.

        Args:
            fig (plt.Figure): The Figure object of the plot to save.
            output_dir (Path): The directory path where the plot will be saved.
            region (str): The genomic region of the plot.
            haplotype (str): The haplotype identifier.
        """
        filename = f"{region}-{haplotype}.png"
        filepath = output_dir / filename
        file_log.info(f"Generation plot {filename}, saving to {filepath}")
        fig.savefig(filepath, dpi=300, bbox_inches='tight', format='png')

    def create_plot(self, df, output_dir, region,
                    haplotype, sample, merged):
        """
        Creates a plot for the given DataFrame and parameters.

        Args:
            df (pd.DataFrame): DataFrame containing the segments and their statuses.
            output_dir (Path): The directory path where the plot will be saved.
            region (str): The genomic region of the plot.
            haplotype (str): The haplotype identifier.
            sample (str): The sample identifier.
            merged (list): The list of merged haplotypes.

        Returns:
            fig (plt.Figure): The Figure object for the plot.
            ax (plt.Axes): The Axes object for the plot.
        """
        if not df.empty:
            block_size, block_spacing, base_width_per_block, font_size = self.calculate_plot_size(
                len(df))
            return self.configure_plot(df, block_size, block_spacing, base_width_per_block, region, haplotype, 'test', sample, merged)

    def load_and_combine_dataframes(self, files):
        """
        Loads multiple Excel files as DataFrames, assigns a status to each, and combines them into a single DataFrame.

        Args:
            files (list): A list of file paths to process.

        Returns:
            combined_df (pd.DataFrame): A combined DataFrame containing all data from the files with assigned statuses.
        """
        dataframes = [self.load_dataframe(filename) for filename in files]
        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df

    def combine_and_save_plots_vertically(self, fig1, fig2, output_path, dpi=300):
        """
        Combines two matplotlib figures vertically and saves the combined figure to the specified path.

        Args:
            fig1 (plt.Figure): The first Figure object.
            fig2 (plt.Figure): The second Figure object.
            output_path (Path): The path where the combined figure will be saved.
            dpi (int): The resolution of the output image.
        """
        def fig_to_numpy(fig):
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
        file_log.info(
            f"Creating {output_path.name} and saving it to {self.cwd /output_path}!")
        plt.savefig(self.cwd / output_path, dpi=dpi, bbox_inches='tight',
                    pad_inches=0, transparent=True)
        plt.close(combined_fig)

    def orientation_calculation(self, df):
        """
        Determines the orientation of the plot based on the first few segments.

        Args:
            df (pd.DataFrame): DataFrame containing segment information.

        Returns:
            orientation (bool): Returns False if the segment starts with J or D, otherwise True.
        """
        segment = df.sort_values(by="Start coord").head(n=3)
        single_segment = segment["Segment"].mode().iloc[0]
        return False if single_segment in ("J", "D") else True

    def run_seperate_haplotype(self, df, region, sample):
        """
        Generates and saves separate plots for each haplotype.

        Args:
            df (pd.DataFrame): DataFrame containing the segment information.
            region (str): The genomic region of the plot.
            sample (str): The sample identifier.
        """
        for hap in ["hap1", "hap2"]:
            hap_df = df.query(f"Haplotype == '{hap}'")
            if not hap_df.empty:
                orientation_check = self.orientation_calculation(
                    hap_df) if self.rotate else True
                hap_df = hap_df.sort_values(
                    by="Start coord", ascending=orientation_check)
                merged = hap_df["Short name"]
                plot = self.create_plot(hap_df, self.output_dir,
                                        region, hap, sample, merged)
                self.save_plot(plot[0], self.output_dir, region, hap)

    def get_haplotype(self, df):
        """
        Retrieves the most common haplotype from the DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing the segment information.

        Returns:
            haplotype (str): The most common haplotype.
        """
        return df["Haplotype"].mode()[0]

    def run_combined_haplotype(self, new_df, region, sample):
        """
        Generates and saves a combined plot for both haplotypes.

        Args:
            new_df (pd.DataFrame): DataFrame containing the segment information for both haplotypes.
            region (str): The genomic region of the plot.
            sample (str): The sample identifier.
        """
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
                plots[0][0], plots[1][0], self.output_dir / f"{region}.png")

    def process_region(self, df, references):
        """
        Process references within a region, aligning haplotypes if conditions are met.

        Args:
            df (pd.DataFrame): DataFrame containing the segment information.
            references (list): List of references to align.

        Returns:
            df (pd.DataFrame): Filtered DataFrame containing only the needed references.
        """
        if len(references) > 2:
            return df[df['Path'].isin(references)]
        return df

    def process_region_data(self, df, region, sample):
        """
        Process data for a specific region, running either single or combined plotting styles.

        Args:
            df (pd.DataFrame): DataFrame containing the segment information.
            region (str): The genomic region to process.
            sample (str): The sample identifier.
        """
        region_df = df.query(f"Region == '{region}'")
        references = region_df['Path'].value_counts()
        needed_references = references.head(2).index.tolist()
        new_df = self.process_region(region_df, needed_references)
        filtered = new_df.copy()
        filtered['Filtered'] = filtered['Short name'].apply(
            lambda x: self.add_filter([x])[0])
        if self.style == "single":
            self.run_seperate_haplotype(filtered, region, sample)
        else:
            self.run_combined_haplotype(filtered, region, sample)


def parse_arguments():
    """
    Parses command-line arguments for the plot generation script.

    Returns:
        args (argparse.Namespace): Parsed arguments object.
    """
    cwd = Path.cwd()
    parser = argparse.ArgumentParser(
        description="Generate plots for VDJ segments.")
    parser.add_argument('-f', '--files', nargs='+', required=True,
                        help="List of file paths to process.")
    parser.add_argument('-o', '--output-dir', type=str,
                        help="Path to the output directory.", default=cwd / "VDJ_display")
    parser.add_argument('-s', '--style', choices=['single', 'combined'],
                        help="The plotting style: 'single' for separate haplotypes, 'combined' for combined haplotypes.", default='combined')
    parser.add_argument('--no-footer', action='store_false',
                        help="Flag to not add footer to plots.", dest='footer')
    parser.add_argument('-r', '--rotate', action='store_true', default=False,
                        help="Set orientation to True, defaults to False if not provided.")
    parser.add_argument('-c', '--color-theme', choices=[
        'default', 'easy_eye', 'IBM', 'Wong', 'Tol'
    ], help="Choose a color theme for the plots.", default='default')
    parser.add_argument('--re-evaluate', action='store_true', default=False,
                        help="Reevaluate the VDJ segments to better determine the order")

    return parser.parse_args()


def main():
    """
    Main function for generating plots. Parses command-line arguments, sets up the plot generator,
    and processes the data for plotting.

    Steps:
        1. Parses the arguments and retrieves the list of input files.
        2. Initializes the PlotGenerator and sets its configuration based on the arguments.
        3. Optionally re-evaluates the input files for VDJ segment order.
        4. Loads and combines DataFrames from the input files.
        5. For each unique genomic region, processes the data and generates plots.
    """
    # Define color themes
    color_themes = {
        'default': {
            "V": '#ffdfba',
            "D": '#baffc9',
            "J": '#bae1ff',
            "Known": '#adb2fb',
            "Novel": '#ffb3ba',
            "": '#ffffff'
        },
        'easy_eye': {
            'V': '#E69F00',
            'D': '#56B4E9',
            'J': '#009E73',
            'Known': '#F0E442',
            'Novel': '#0072B2',
            '': '#FFFFFF'
        },
        'IBM': {
            'V': '#785EF0',
            'D': '#DC267F',
            'J': '#FE6100',
            'Known': '#FFB000',
            'Novel': '#648FFF',
            '': '#FFFFFF'
        },
        'Wong': {
            'V': '#E69F00',
            'D': '#F0E442',
            'J': '#CC79A7',
            'Known': '#009E73',
            'Novel': '#D55E00',
            '': '#FFFFFF'
        },
        'Tol': {
            'V': '#DDCC77',
            'D': '#332288',
            'J': '#88CCEE',
            'Known': '#117733',
            'Novel': '#882255',
            '': '#FFFFFF'
        },
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
    plot_generator.legend_labels = {
        key: key for key in plot_generator.colors if key != ''}
    reevaluated_excel = Path(Path.cwd() / "reevaluated.xlsx")
    if args.re_evaluate:
        file_log.info(
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
