
![Logo](images/BPRC_logo.png)

# VDJ Analize, Assamble and Annotate Pipeline (VDJ-AAAP)

## Table of content

- [VDJ Analize, Assamble and Annotate Pipeline (VDJ-AAAP)](#vdj-analize-assamble-and-annotate-pipeline-vdj-aaap)
  - [Table of content](#table-of-content)
  - [Authors](#authors)
  - [Abstract](#abstract)
  - [Installation](#installation)
    - [Base setup](#base-setup)
    - [Environment Setup](#environment-setup)
    - [Getting Started](#getting-started)
  - [Overview tool:](#overview-tool)
  - [Pipeline:](#pipeline)
    - [Detailed Flag Descriptions](#detailed-flag-descriptions)
    - [Required:](#required)
    - [Optional:](#optional)
    - [Important Notes](#important-notes)
  - [Annotation](#annotation)
    - [Required:](#required-1)
    - [Data source (mutually exclusive):](#data-source-mutually-exclusive)
    - [Assembly specific flags:](#assembly-specific-flags)
    - [Optional flags:](#optional-flags)
    - [Important Notes](#important-notes-1)
  - [Configuration settings](#configuration-settings)
    - [Basic configuration options](#basic-configuration-options)
    - [Important configuration settings](#important-configuration-settings)
  - [Pipeline output](#pipeline-output)
    - [Annotation results](#annotation-results)
  - [Plots](#plots)
    - [Single plots](#single-plots)
    - [Interactive plot](#interactive-plot)
      - [Activation](#activation)
      - [Options](#options)
      - [Interactive plot demo](#interactive-plot-demo)
  - [Acknowledgements](#acknowledgements)

## Authors

- [@Jesse mittertreiner](https://github.com/AntiCakejesCult)
- [@Giang Le](https://github.com/GiangLeN)

## Abstract

The VDJ-AAAP pipeline offers a robust framework for analyzing, assembling, and annotating long sequence reads from Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT). Designed to uncover both novel and known VDJ segments within T-cell receptors (TCR) or B-cell receptors/immunoglobulins (IG), this part is still in beta. This versatile tool supports analysis across various species, given the availability of a reference genome and assembly report on NCBI.

Addressing the challenge of assembling the repetitive VDJ regions, the pipeline selectively processes reads exceeding 5 Kbs, minimizing erroneous mappings and noise. It generates a refined reference genome if none is specified. If a reference genome is present it can alse be selected. From the reference genome we use only known chromosomes and unplaced chromosomal fragments to prevent assembly inaccuracies. Utilizing minimap2, reads are mapped against this filtered reference genome. From the mapping files we identify the reads corresponding the TCR/IG chromosomes and convert it to fastq. 

The assembly phase leverages a hybrid approach meaning it utilizes the strengths of both PacBio and ONT reads, combining PacBio's precision with ONT's extensive read lengths for superior assembly outcomes. The files are loaded into the assembly program seperatly.

Specific genomic regions are then isolated using predefined flanking genes, these can be specified by the user or a default can be used (Default is for primates only because it is based on humans). 

The pipeline extends its functionality to annotation, which can also be ran on its own as part of this tool. Comparing the assembled region sequences against a comprehensive library of validated VDJ segments sourced from the IMGT database or a pre-existing library in **`library/library.fasta`**.
 
Outputs include detailed Excel reports on identified sequences, supplemented with an interactive plotting feature for enhanced data visualization and analysis, but also static plots showing the different haplotypes of the regions showcasing the order of the segements. Lastly it generates a easy to navigate HTML report containing al the generated results.

## Installation

### Base setup

To use VDJ-AAAP, you need to have Conda ([Conda installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)) and Python3 ([Python installation](https://realpython.com/installing-python/)) installed.

### Environment Setup

VDJ-AAAP requires a specific environment to run . You can install all necessary dependencies by setting up a Conda environment using the provided pipeline.yaml file.

  1. Open your terminal.
  2. Run the following command to create the Conda environment:

``` bash
conda env create -f pipeline.yaml --name pipeline
```
If this does not work use mamba instead.
``` bash
mamba env create -f pipeline.yaml --name pipeline
```
**--name pipeline**, can be replaced with the name of your choosing.


This command installs the following critical packages within the environment:

 1. mamba: A fast, flexible package manager that extends conda.
 2 . snakemake: A workflow management system that helps to create and manage bioinformatics pipelines.

### Getting Started

To run the VDJ-AAAP pipeline, follow these steps:

 1. Download the pipeline source code from the GitHub repository:
 [TCR_macaque at BPRC-CGR](https://github.com/BPRC-CGR/TCR_macaque). Use the “Code” button and then select "Download ZIP".
 1. Extract the downloaded TCR_macaque-main.zip file to your desired location.
 2. Open your terminal and navigate to the extracted TCR_macaque-main directory. You can do this with the cd command followed by the path to the directory.
 3. Activate the Conda environment you created earlier by running:

``` bash
conda activate pipeline
```

## Overview tool:
Within this tool, you have the option to run the whole pipeline, which generates the VDJ analysis from the input ONT and PacBio fastq data. You only need to specify the reference genome either as a NCBI genome code or as a fasta file; the type of receptor `TR` or `IG` and the species name of the animal you want to analyse. To see the options that are available enter:
```bash
python scripts/pipeline_shell.py -h
```
The following options will be shown [pipeline](#pipeline) and [annotation](#annotation). You can read more information about these tools by pressing on the highlighted names.
```txt
Tool for sequencing data processing and VDJ annotation

positional arguments:
  {pipeline,annotation}
                        Available commands
    pipeline            Run the pipeline for sequencing data processing.
    annotation          Run the annotation tool for VDJ segment analysis.

options:
  -h, --help            show this help message and exit
```

## Pipeline:
You can execute the pipeline with the following command:

``` bash
python scripts/pipeline_shell.py pipeline -ont <nanopore_data.fastq.gz> -pb <pacbio_data.fastq.gz> -ref <reference_genome> -r <receptor_type> -s <species_name> -t <threads> --default
```

### Detailed Flag Descriptions

### Required:
- `-ont <nanopore_data.fastq.gz>`, `--nanopore <nanopore_data.fastq.gz>`: This flag specifies the path to the Oxford Nanopore Technologies (ONT) reads file. The file must be in `.fastq.gz` format.
  
- `-pb <pacbio_data.fastq.gz>`, `--pacbio <pacbio_data.fastq.gz>`: This flag specifies the path to the Pacific Biosciences (PacBio) reads file. The file must be in `.fastq.gz` format.
  
- `-ref <reference_genome>`, `--reference <reference_genome>`: This optional flag specifies the path to the reference genome file if it is already available. The file can be a `.fasta` or `.fna` file or a valid accession code.
  
- `-r <receptor_type>`, `--receptor-type <receptor_type>`: This required flag specifies the type of receptor to analyze. The options are `TR` for T-cell receptor or `IG` for Immunoglobulin.
  
- `-s <species_name>`, `--species <species_name>`: This required flag specifies the scientific name of the species, e.g., "Homo sapiens".
### Optional:
- `-f <flanking_genes>`, `--flanking-genes <flanking_genes>`: This flag specifies a comma-separated list of flanking genes, e.g., `MGAM2,EPHB6`. The genes should be added as pairs.
  
- `-c <chromosomes>`, `--chromosomes <chromosomes>`: This flag specifies a list of chromosomes where TR or IG is located. All values must be integers between 1-22, or 'X', 'Y'.
  
- `-t <threads>`, `--threads <threads>`: This optional flag specifies the number of processing threads to use for the analysis. If not specified, it defaults to 8.
  
- `--default`: This optional flag uses default settings for the analysis. It cannot be used in conjunction with `-f/--flanking-genes` or `-c/--chromosomes`.

### Important Notes

- Ensure that the paths to the input files are correct and accessible.
- If using the `--default` flag, do not specify `-f/--flanking-genes` or `-c/--chromosomes` as they are mutually exclusive with `--default`.

## Annotation
Within this tool there is also a option to run the annotation tool as a stand-alone program. It is the same program that is used within the pipeline. You can execute it by running this command.

```bash
python scripts/pipeline_shell.py annotation -a <assembly_directory> -l <library.fasta> -r <receptor_type> -s <species_name> -f <flanking_genes> -t <threads>

```
### Required:
- `-l <library.fasta>`, `--library <library.fasta>`: This flag specifies the path to the library FASTA file. This is required. (default: None)

- `-r <receptor_type>`, `--receptor-type <receptor_type>`: This flag specifies the type of receptor to analyze. The options are `TR` for T-cell receptor or `IG` for Immunoglobulin. This is required. (default: None)

### Data source (mutually exclusive):
Select the data source: regions or assembly.

- `-i <input_directory>`, `--input <input_directory>`: This flag specifies the directory containing the extracted sequence regions in FASTA format, where VDJ segments can be found. Cannot be used with `-f/--flanking-genes` or `-s/--species`. (default: None)

- `-a <assembly_directory>`, `--assembly <assembly_directory>`: This flag specifies the path to the directory containing the assembly FASTA files. This is required if not using `-i/--input`. Must be used with `-f/--flanking-genes` and `-s/--species`. (default: None)

### Assembly specific flags:
These flags are required if `-a/--assembly` is chosen:

- `-f <flanking_genes>`, `--flanking-genes <flanking_genes>`: This flag specifies a comma-separated list of flanking genes, e.g., `MGAM2,EPHB6`. Add them as pairs. Required with `-a/--assembly`. (default: None)

- `-s <species_name>`, `--species <species_name>`: This flag specifies the scientific name of the species, e.g., "Homo sapiens". Required with `-a/--assembly`. (default: None)

### Optional flags:
- `-o <output_directory>`, `--output <output_directory>`: This flag specifies the output directory for the results. If not specified, it defaults to `annotation`. (default: annotation)

- `-m <mapping_tool>`, `--mapping-tool <mapping_tool>`: This flag specifies the mapping tool(s) to use. Choose from: `minimap2`, `bowtie`, `bowtie2`. Defaults to all. (default: ['minimap2', 'bowtie', 'bowtie2'])

- `-t <threads>`, `--threads <threads>`: This flag specifies the number of processing threads to use for the analysis. If not specified, it defaults to 8. (default: 8)

### Important Notes

- Ensure that the paths to the input files are correct and accessible.
- If using the `-i/--input` flag, do not specify `-f/--flanking-genes` or `-s/--species` as they are only need when using `-a/--assembly`.

## Configuration settings

The **config.yaml** file, located within the **config** directory, serves as a overview of the configuration settings that the pipeline uses. This config file is generate automatically and contains various parameters to tailor the analysis based on your given input. The config file is also different when running the pipeline or only the annotation program. For the annotation only the **SPECIES**, **FLANKING_GENES**, and the **important** settings are generated.

Although it is automatically generated, it is recommended to see what it contains.

### Basic configuration options

Below is an overview of the file's structure and the basic options you can configure:

``` yaml
ALL_CHROMOSOMES:
  - 1
  - "X"
  - "Y"

ASSEMBLY_CHROMOSOMES:
  - 3
  - 7

HAPLOTYPES:
  - 1
  - 2

SPECIES:
  name:
 "macaca mulatta"
  genome:
 GCF_003339765.1
  cell:
 TR

FLANKING_GENES:
- SALL2
- DAD1
- MGAM2
- EPHB6
- EPDR1
- VPS41
```
- **DATA**: The used input data for both ONT and PacBio samples. This is separated in the original files and the moved files, these files are always moved in a temporary downloads folder. 
- **All_CHROMOSOMES**: This list contains all chromosomes found within the reference genome. For those uncertain of the specific chromosomes included in their reference genome, it's advisable to consult the [NCBI Genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/).      
- **ASSEMBLY_CHROMOSOMES**: List contains the TCR/IG chromosomes for the analysis.
- **HAPLOTYPES**: List with each haplotype needed the analysis.
- **SPECIES: name**: The chosen species for the analysis.
- **SPECIES: genome**: Specified genome that is being used for the run. It can be found on the [NCBI Genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/). Or already present.
- **SPECIES: cell**: The receptor type that is being analyzed. This value to either `TR` (T-cell receptor) or `IG` (Immunoglobulin), depending on your study's focus.
- **FLANKING_GENES**: A list of flanking genes that are needed to identify the contig containing the different TCR/IG regions. By default these are determined based on the given receptor type.

### Important configuration settings

To enable accurate validation of the VDJ gene segments, specifying the RSS layout is crucial. This ensures proper validation can be performed. These settings are also chosen based on the receptor type.

```yaml
RSS_LAYOUT:
  TRAV:
 "23":
   "+": end_plus
   "-": start_minus
  TRAJ:
 "12":
   "+": start_minus
   "-": end_plus

RSS_LENGTH:
  "12": 28
  "23": 39

RSS_MERS:
  "12":
 - 9
 - 7
  "23":
 - 7
 - 9

```

- **RSS_LAYOUT**: Dictates the identification approach for Recombination Signal Sequences (RSS) in relation to the VDJ gene segments. This configuration is crucial for determining the precise location of RSS for accurate gene segment annotation. Within this setting, you specify:
  - The segment type (e.g., `TRAV`, `TRAJ`) to configure.
  - The RSS type (`12`, `23`) indicating the spacer length in base pairs.
  - The orientation (`+`, `-`) for identifying the RSS direction in relation to the gene segment.
  - Methods of extraction (`start_minus`, `end_plus`) which define how the RSS is located and annotated based on its positional context to the VDJ segment.
   	- `start_minus`: Extracts the RSS from the start (left side) of the VDJ segment, utilizing the segment's start coordinate and subtracting the RSS length to locate the RSS.
   	- `end_plus`: Extracts the RSS from the end (right side) of the VDJ segment, using the segment's end coordinate and adding the RSS length to pinpoint the RSS.

- **RSS_LENGTH**: Specifies the length of each RSS type, classified by the spacer length (`12`, `23`). This length is essential for correctly identifying and annotating the RSS within the genomic sequence. The numbers (e.g., `28`, `39`) represent the total length in base pairs for each RSS type.

- **RSS_MERS**: Defines the positions of key components within the RSS - specifically, the heptamer and nonamer elements, denoted as `7` and `9`, respectively. This configuration allows for detailed specification of each RSS's structural components, critical for the annotation process. The lists under each RSS type (`12`, `23`) enumerate the preferred positions of these elements, facilitating precise identification and analysis of RSS structures within the genomic data.

## Pipeline output

The pipeline creates a lot of important files locateded in different directories. The next code sample shows the tree with all the directories that are created when running pipeline.

```txt
.
├── alignments
├── annotation
├── assembly
├── benchmarks
├── BUSCO
├── busco_downloads
├── chromosomes
├── config
├── converted
├── final
├── flank_genes
├── inspector
├── library
├── logs
├── mapped_genes
├── mapping
├── merged
├── QC
├── quast
├── reference
├── region
├── RSS
├── source/html
└── split_files
```

### Annotation results

One of the most important parts of the pipeline is the finding of novel VDJ gene segments. The result of the findings are located in the folder called **annotation**. In this directory, are the following excel files located.

```txt
annotation/
├── annotation_report_100%.xlsx
├── annotation_report_100%_plus.xlsx
├── annotation_report_long.xlsx
├── annotation_report_plus.xlsx
├── annotation_report.xlsx
├── blast_results.xlsx
└── report.xlsx
```

- **report.xlsx**: In the report file are all the initial mapping results, to get a initial understanding of the amount of VDJ gene segments that are identified. This includes non-novel and novel segments.
- **blast_resutls.xlsx**: In the blast result are all the reevaluated segments. This includes the deviation between found segment and the most similar.
- **annotation_report.xlsx**, **annotation_report_100%.xlsx**, **annotation_report_long.xlsx**: In the intial annotation report are all the novel segments that are retained after the filtering of the segments. The 100% version of the annotation report is almost the same as the original report, but this includes only the known segments. Lastly, the long format is the uncondensed version, where the similar sequences are not combined in one row.
- **annotation_report_plus.xlsx** and **annotation_report_100%_plus.xlsx**: Lastly, the annotation report plus reports contain validation columns based on the RSS types of the segments. This indicates the found RSS heptamer and nonamer for a given segment and the RSS heptamer and nonamer that were used for comparison. The final report contains the flowing columns.

| Column                        | Explanation                                                                                                                                                                                                                                   |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Reference**                 | Name of the closest reference for the new segment.                                                                                                                                                                                            |
| **Old name-like**             | New segment's name, derived from the closest reference, appended with "like" to indicate similarity.                                                                                                                                          |
| **Mismatches & % Mismatches** | Number of mismatches with the reference and their percentage relative to the total length of the reference.                                                                                                                                   |
| **Start and End coord**       | Coordinates of the segment within the region of interest.                                                                                                                                                                                     |
| **Function**                  | Segment function: functional (F), open reading frame (ORF), or pseudogene (P) if an early stop codon is detected.                                                                                                                             |
| **Similar references**        | Other references for a segment with the same start and end coordinates; the best match is selected based on mutation count and reference name.                                                                                                |
| **Path**                      | Path to the FASTA file of the region containing the segment.                                                                                                                                                                                  |
| **Strand**                    | Orientation of the segment: 5' to 3' (`+`) or 3' to 5' (`-`).                                                                                                                                                                                 |
| **Region and Segments**       | Type of region and segment identified.                                                                                                                                                                                                        |
| **Haplotype**                 | Haplotype (1 or 2) on which the segment is found.                                                                                                                                                                                             |
| **Sample**                    | Name of the sample providing the genetic data.                                                                                                                                                                                                |
| **Short name**                | Only the part of the segment that includes the region, segment and variant.                                                                                                                                                                   |
| **RSS**                       | Each RSS spacer type includes six columns, with three dedicated to both the heptamer and nonamer segments. These columns represent the segment sequence, a reference sequence, and a boolean indicating if the segment matches the reference. |

## Plots

This pipeline also creates individual plots and a interactive plot to showcase the results.

### Single plots
The single plots are created with a custom V(D)J display tool using the **annotation_report_plus.xlsx** and **annotation_report_100%_plus.xlsx** excel files. These single plots are created for every region and haplotype. There is a choice the show both haplotypes of a region at once or as single plot per haplotype. For more information, please read the [README](https://github.com/BPRC-CGR/VDJ_display/tree/development) on its GitHub page.
![Logo](images/single.png)

### Interactive plot

The interactive plot is automatically generated based on the results in the **annotation_report_plus.xlsx**.

#### Activation

Launch the interactive visualization by executing the following steps:

1. **Open Terminal**: Navigate to the `TCR_macaque-main` directory using your command line interface.

2. **Start Bokeh Server**: Enter the command below and press `Enter`:

 ```bash
 bokeh serve scripts/visualisation.ipynb
 ```

 Upon execution, you should see messages similar to these:

 ```
 2024-03-08 10:08:34,162 Starting Bokeh server version 3.3.0 (running on Tornado 6.3.3).
 2024-03-08 10:08:34,163 User authentication hooks NOT provided (default user enabled).
 2024-03-08 10:08:34,166 Bokeh app running at: http://localhost:5006/visualisation.
 2024-03-08 10:08:34,166 Starting Bokeh server with process id: 3759457.
 ```

3. **Access the Visualization**: Click on [http://localhost:5006/visualisation](http://localhost:5006/visualisation) or copy and paste this URL into your web browser.

4. **Deactivate**: To deactivate the application, press `⌃C` (Control + C) on MacOS or `Ctrl + C` on Windows in the command line interface.

#### Options

The interactive plot includes several controls to customize the display:

- **Region-Haplotype**: This dropdown menu lists regions containing different segments for various haplotypes. When choosing `All`, all regions will be shown.
- **Function**: A green toggle button following the dropdown. It allows for switching between displaying only segments classified as functional (F/ORF) and showing all segments (F/ORF & P).
- **New**: Another green toggle button enabling the option to display only novel segments or both novel and non-novel segments.

#### Interactive plot demo

![Demo](images/visual.gif)

## Acknowledgements

I would like to thank [Jesse Bruijnesteijn](https://github.com/JesseBNL) and [Susan Ott](https://github.com/SusanOtt) for their contributions and insights, which have significantly enhanced the pipeline. Your expertise and suggestions have been really helpful in improving the pipelines functionality and effectiveness. Thank you both for your dedication and support.
