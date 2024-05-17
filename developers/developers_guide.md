![Logo](../images/BPRC_logo.png)

# Developers guide VDJ-AAAP pipeline

## Table of Contents

- [Developers guide VDJ-AAAP pipeline](#developers-guide-vdj-aaap-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Authors](#authors)
  - [Introduction](#introduction)
  - [Packages](#packages)
  - [Pipeline](#pipeline)
  - [Custom scripts and tools](#custom-scripts-and-tools)
  - [region.py](#regionpy)
  - [**IMGT Scraper**](#imgt-scraper)
    - [Usage](#usage)
  - [**Annotation**](#annotation)
    - [annotation.py](#annotationpy)
    - [mapping.py](#mappingpy)
    - [write\_annotation\_report.py](#write_annotation_reportpy)
    - [RSS.py](#rsspy)
    - [Usage](#usage-1)
  - [V(D)J display](#vdj-display)
    - [VDJ\_display.py](#vdj_displaypy)
    - [reevaluate.py (BETA)](#reevaluatepy-beta)
    - [order\_segments.py (BETA)](#order_segmentspy-beta)
    - [Usage](#usage-2)

## Authors

- [@Jesse Mittertreiner](https://github.com/AntiCakejesCult)
- [@Giang Le](https://github.com/GiangLeN)

## Introduction

For a brief overview of the main working of the VDJ-AAAP pipeline, please read the [abstract](../README.md#abstract). In this document, we will outline the main working of the different created script and pipeline. We discuss parts of the scrips and pipeline that can be changed, but not with the [config file](../README.md#configuration-settings). Furthermore, the main packages and modules used for the envs are shown.

## Packages

This pipeline and its corresponding scripts use different environments to generate the different results and intermediate results. If we look closer in to environment files, you see that packages are duplicated in different environment files. This is because you cannot inherit packages from environments, and Snakemakes **Conda** argument does not allow for multiple environment files.

This list contains all the different used packages and modules and which version are used. Changing the version of the packaged can be done in the individual envs.

| Package          | Version   |
| ------------------- | --------- |
| **beautifulsoup4**  | 4.12.3 |
| **inspector**    | 1.2    |
| **matplotlib**   | 3.8.4  |
| **python**       | 3.11   |
| **python**       | 3.6    |
| **requests**     | 2.31.0 |
| **bedtools**     | 2.31.0 |
| **biopython**    | 1.83   |
| **blast**        | 2.14.1 |
| **bowtie**       | 1.3.1  |
| **bowtie2**      | 2.5.1  |
| **busco**        | 5.5.0  |
| **hifiasm**      | 0.19.8 |
| **imagemagick**  | 7.0.11_12 |
| **jq**           | 1.5    |
| **meme**         | 4.11.2 |
| **minimap2**     | 2.26   |
| **ncbidatasetscli** | 15.25.0   |
| **openpyxl**     | 3.1.2  |
| **pandas**       | 2.2.1  |
| **pyyaml**       | 6.0.1  |
| **quast**        | 5.2.0  |
| **samtools**     | 1.6    |
| **seqkit**       | 2.8.0  |
| **seqtk**        | 1.4    |
| **unzip**        | 6.0    |

## Pipeline

For the pipeline, just like mentioned, most things can be changed with the config file. But here are some examples of things that be changed in the **snakefile** itself. The main [config file](../Snakefile4#L6) can be changed. If changed, please change it as well in the following python scripts [region](../scripts/region.py), [mapping](../scripts/mapping.py), [RSS](../scripts/RSS.py) and [write_annotation_report](../scripts/write_annotation_report.py).

You can add extra values in the config file, but make sure you add them in the [snakefile](../Snakefile4#L31) and other scripts like [this](../scripts/RSS.py#L125). The input directory can also be changed if needed, please rename every **download** to something else you like. Every command in the snakefile can also be altered to your liking. Please follow the requirements of the package. For the in-house IMGT scrape, annotation tools and V(D)J display tool, please look at the **help rule (-h)**.

From previous runs, it is known that sometimes the **kmer** and **window** needed to be added to achieve better results. These parameters can be added to the [command](../Snakefile4#L365).

## Custom scripts and tools

This pipeline uses a combination of different custom created python scripts and tools to produce results.

## region.py

The script extracts of specific regions from FASTA files using flanking genes defined in the config file. It parses the alignment data from SAM files to identify the best mapping coordinates for these regions based on the bitwise flag. The final output is a set of region FASTA files, each containing the sequence of a specified region. This process is managed within a Snakemake. If wanted, you could change the [snakemake.wildcards](../scripts/region.py#L165) to something that fetches the sample name from the file itself. Lastly, you could change [snakemake.input[0]](../scripts/region.py#L191) to the name of the SAM file, then this script could also be used without the pipeline.

## **IMGT Scraper**

This script/tool automates the retrieval of immunoglobulin (IG) and T-cell receptor (TR) V(D)J segment sequences from the IMGT database for specified species. It allows users to customize the sequence fetching process based on several parameters including species type, receptor type (IG or TR), and sequence frame, which is needed to change the types of segments to fetch. The script outputs the sequences in FASTA format per region and segment. Optionally, it can compile them into a comprehensive library and remove the excess files.

### Usage

```bash
python imgt_scrape.py -S "Homo sapiens" -T IG --output /path/to/output --create-library --cleanup
```

For more information about this tool, please see the [README.md](https://github.com/BPRC-CGR/IMGT_scrape/tree/development)

## **Annotation**

This is the main tool or set of scripts that is used to generate annotation to discover novel and known V(D)J segments, present it the input data. It takes a list of regions and a library with the help of config file it generates different [annotation results](../README.md#annotation). By default, it uses minimap2, bowtie and bowtie2 to do annotation. This is a quick and simple overview of the core analysis it does.

![flowchart](../images/annotation_flowchart.png)

### annotation.py

`annotation.py` is the main script activated when running the annotation tool. It starts by checking user input, then calls the mapping script to perform initial mapping using user-specified tools.

The script evaluates the mapped segments with BLAST to identify differences between the segments and the reference library. The BLAST command and its variables can be modified in [`construct_blast_command`](../scripts/annotation.py#construct_blast_command), which includes adjustments for smaller segments, such as the D segment of the TCR (10 to 25 bases), to facilitate easier verification. These settings can be [altered](../scripts/annotation.py#L140) to accommodate different segment lengths. BLAST is executed three times with different [cutoffs](../scripts/annotation.py#:203). Information required from the FASTA header for BLAST is [parsed](../scripts/annotation.py#L169) after mapping.

Following the BLAST analysis, the script writes an initial Excel file `blast_result.xlsx` and validates the identified segments using the RSS script.

### mapping.py

This script, `mapping.py`, is the core component for initiating the genomic mapping process. It begins by ensuring the necessary configurations are properly loaded from the config file, essential for extracting the right information. If the configuration file is absent, the script terminates, as indicated by the error handling in [`load_config`](../scripts/mapping.py#L15).

Following the initial setup, the script organizes the necessary files for genomic mapping using the `MappingFiles` class, which is designed to manage file paths tailored for the used mapping tools such as Bowtie, Bowtie2, or Minimap2. [`MappingFiles constructor`](../scripts/mapping.py#L26).

With a dynamic command construction for executing mapping processes tailored to specific accuracy requirements and conditions. For example, [`make_bowtie_command`](../scripts/mapping.py#L130) and [`make_bowtie2_command`](../scripts/mapping.py#L170), can change their command based on the parsed acc score. Changing the parameters in these function will result in other mapping results. If the

The V(D)J sequence themselves is extracted with [`get_sequence`](../scripts/mapping.py#get_sequence). They are retrieved from the FASTA file based on coordinates provided in a BED file. The function, [`get_region_and_segment`](../scripts/mapping.py#L68), categorizes data into biological regions and segments based on predefined [`criteria`](../config/config.yaml#L42), adding structure to the analysis. Important that the headers in the library are separated with a `"_"` and that the region and segment are joined, like this **`TRAV`**.

Results from the mapping are parsed and compiled. The compiled data is then transformed into a pandas DataFrame by [`make_df`](../scripts/mapping.py#L255). Then it is saved as `blast_results.xlsx`

### write_annotation_report.py

The data from the `blast_results.xlsx` is processed with `write_annotation_report.py`. It takes the results from `blast_results.xlsx`. Just as before, the script loads the config file using the [`load_config`](../scripts/write_annotation_report.py#L10). This step ensures all defined settings are loaded and can be used.

The script then constructs a record dictionary from the library FASTA file through the [`make_record_dict`](../scripts/write_annotation_report.py#L17) function. This dictionary is vital for referencing the library V(D)J sequences.

Additionally, the script enhances the dataset through the [`add_values`](../scripts/write_annotation_report.py#L77) function, which adds several new columns such as mismatch percentages and sequence lengths. This key information is needed to filter more later in the script.

Then we divide the entries in `blast_results.xlsx` in to two groups, with the [`main_df`](../scripts/write_annotation_report.py#L80) function. It is split into a dataset with only novel entries, which have at least one mismatch, and a dateset with only known entries containing entries.

An important part of the script is the [`filter_df`](../scripts/write_annotation_report.py#L107) function, which refines the dataset by selecting the most relevant reference sequences based on specific criteria:

- **Specific Part Identification**: It combines "Region" and "Segment" values from each row to form a "Specific Part," used to identify relevant sequences that match the query conditions closely.
- **Reference Selection**: This function filters rows based on the presence of this specific part in the "Reference" column, prioritizing entries where the sequence length of the reference matches that of the query sequence.
- **Optimal Reference Determination**: Among the filtered rows, it selects the one with the fewest mismatches, sorting by the reference identifier to choose the best match. If no specific matches are found, it defaults to the first row.
- **Compilation of Similar References**: After selecting the best reference, it compiles a list of all other references that were not selected, providing context for potential alternative annotations.

Finally, the script prepares and formats the refined data for exporting to an Excel files. The function [`annotation_long`](../scripts/write_annotation_report.py#L283) exports the data in not filtered way to an Excel file called `annotation_report_plus.xlsx` and [`annotation`](../scripts/write_annotation_report.py#L325) exports the data into filtered reports. One of the novel matches in `annotation_report.xlsx` and the know matches in `annotation_report_100%.xlsx`.

### RSS.py

To validate the entries in the `annotation_report.xlsx` and `annotation_report_100%.xlsx` files, we use the RSS of the found segments, we use `RSS.py`. It processes and analyze RSS from the `annotation_report_100%.xlsx`.

This time the config file is really important. It is used to determine how the RSS should be extracted and what the length is of the RSS. This is different per V(D)J segment on the different region for TCR and IG. Please read the [`Important configuration settings`](../README.md#important-configuration-settings) for a deeper understanding of how it works and how to alter it. Using the [`load_config`](../scripts/RSS.py#L22) function, these settings are loaded.

It first creates the needed directories with the [`create_directory`](../scripts/RSS.py#L30) function.

We first run [`create_all_RSS_meme_files`](../scripts/RSS.py#L594). Here we create all the needed RSS extract the needed RSS for three different ["types"](../scripts/RSS.py#L619). For reference_RSS (Known), new_RSS (Novel) and combined_RSS (Known + Novel). Based on the segment and region the right values are retrieved from the config file. Using the [`calculate_position`](../scripts/RSS.py#L77) function. The start coord and end coord are adjusted to get the RSS coordinates. With [`write_fasta_file`](../scripts/RSS.py#L37) the RSS is saved into FASTA file. Each sequence file is named and organized according to its  region and segment (**TRAV.fasta**).

Then it generates MEME motif files through the [`run_meme`](../scripts/RSS.py#L169) function. This function handles the execution of the MEME suite, a tool set for motif discovery, adjusting [command](../scripts/RSS.py#L301) parameters if the amount of RSSs in the input file is higher than one. After running the MEME suite, we generated the sequence motifs and belonging text file `meme.txt`. From each RSS sequence motif regex pattern, the heptamer and nonamer are extracted with the help of a [regex pattern](../scripts/RSS.py#L334) (`\[[^\]]*\]|.`). It matches all substrings within square brackets and any single characters outside of brackets.

With the `make_reference_rss` we create a reference dictionary of RSS motifs. It processes the `meme.txt` output files to extract RSS regular expressions, for sequence validation. This dictionary is essential for comparing newly discovered RSS sequences against a reference.

When the reference RSS regex patterns are established, we use the earlier loaded datasets and extract the RSS just as before. Now we append the RSS to the dataset per entry in the dataset through the use of [`add_base_rss_parts`](../scripts/RSS.py#L254).

The script also integrates data validation features through the [`check_ref_rss`](../scripts/RSS.py#L286) function. This function compares newly generated RSS sequences against reference motifs to validate their accuracy, storing the results within the dataset. Every RSS is validated per position. When the number of mismatches is greater than 1, we put in `False`, which means there is a likely change the RSS is not correct. Otherwise, we put it as `True`.

In the final step, data is processed data for exporting to Excel. The [`create_rss_excel_file`](../scripts/RSS.py#L346) function compiles the results. Then with [`combine_df`](../scripts/RSS.py#L491) the new RSS data is merged with the old dataset. The following columns are [merged](../scripts/RSS.py#L509). Finally, the data is exported to an Excel file (`annotation_report_plus.xlsx` & `annotation_report_plus.xlsx`).

### Usage

To use the tool, you can use the following command. More information about the tool's command and settings are available on the [`README`](https://github.com/BPRC-CGR/VDJ_segment_discover/tree/devlopment) on GitHub page for the tool itself.

```bash
python annotation.py --input /path/to/sequence_data_directory --library /path/to/reference_library.fasta --output /path/to/output_directory
```

## V(D)J display

Tool to nicely visualize V(D)J segments on their haplotypes. Giving the ability to show the two haplotypes at once or all single. The tool also has the ability to redetermine the order of the segments using **`--re-evaluate`**, this option is still in beta.

### VDJ_display.py

The `VDJ_display.py` script is essential for generating various V(D)J segment plots. It starts by processing user input with the [`parse_arguments`](../scripts/VDJ_display.py#L422) function. Next, the file paths are cleaned up by removing any extraneous separators.

The script uses the `PlotGenerator` class to handle everything from data loading to the final visualization output. Upon initialization, `PlotGenerator` sets up necessary configurations such as color schemes corresponding to different genomic segments and prepares the output directory for storing generated plots.

Extra color schemes can be added by creating an extra a dictionary in [`color_themes`](../scripts/) dictionary. Make sure you also add it in the [`choices`](../scripts/) in the `--color-theme`.

Data processing begins with the `load_dataframe` function that parses data from provided Excel files into a single dataset. This dataset is then processed by unique regions through the [`process_region_data`](../scripts/VDJ_display.py) method.

The [`process_region`](../scripts/VDJ_display.py) function checks for any entries with misclassified files, where segments appear in files where they logically shouldn't. Simultaneously, the `add_filter` function refines the segment names by removing unnecessary characters, making the data easier to manage for plotting.

Depending on user settings, the script runs either [`run_seperate_haplotype`](../scripts/VDJ_display.py) for single haplotype visualizations or [`run_combined_haplotype`](../scripts/VDJ_display.py) for combined plots. The data is sorted by "Start coord" for plot creation. For combined configurations, `align_lists_by_matching` is used to synchronize two lists of genomic segments by inserting blanks to address mismatches.

Footers are dynamically determined; for single plots, `Short names` are used directly, while for combined plots, the [`generate_footer_lists`](../scripts/VDJ_display.py) method creates two lists for each haplotype based on start coordinates—one where segment values match and another where they differ. This also uses the `Short names`.

The `calculate_plot_size` function adjusts plot dimensions based on the number of segments. It returns the **block size, block spacing, width per block, font size**. The plot structure is then created by the [`configure_plot`](../scripts/VDJ_display.py) function, which calls the functions to calculate the total width needed and sets up the plotting area with `initialize_plot`. The `add_segment_blocks` method creates two blocks using `create_dual_color_features`, differentiating 'V', 'D' and 'J' for the top block and 'Novel' and 'Known' for the bottom block.

Legends are added through `add_legends`, which creates a legend for sample and haplotypes and other for the colored blocks. Additionally, an optional footer can be added to provide further details or context about the data being visualized with [`add_footer`](../scripts/).

The `save_plot` method ensures each plot is saved in the designated directory. For combined plots, `combine_and_save_plots_vertically` merges the two plots into a single image by converting them to NumPy arrays, aligning them vertically, and then saving the unified plot.

### reevaluate.py (BETA)

### order_segments.py (BETA)

### Usage

To use this tool, you can use the following command. For more information, please read the [README](https://github.com/BPRC-CGR/VDJ_display/tree/development) on its GitHub page.

```bash
python VDJ_display.py --files path/to/file1.xlsx path/to/file2.xlsx --output_dir path/to/output --style combined
```
