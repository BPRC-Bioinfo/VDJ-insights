# TCR Macaque Snakemake Pipeline README

## Overview

This Snakemake pipeline is designed for analysing TCR region of the Rhesus Macaque from publicly available data. It includes a set of rules for downloading SRA data, performing initial QC, removing adaptors, removing duplicates and filtering reads based on length.

## Dependencies

To run this pipeline, make sure you have the following software and environments:
### Already installed
1. Conda (for managing environments)
### Need to be install
1. Snakemake
2. Singularity (for containerization)
3. Questionairy

### Installation Steps

To install the required packages, you can use a Conda environment defined in **pipeline.yaml**. Run the following command to install the packages:

    conda env create -f pipeline.yaml


## Pipeline Rules

**SRA_download**: Downloads data from the SRA database. The rule switches between "wget" and a custom downloader based on the number of files.

**seqkit**: Generates basic statistics for the raw reads.

**pbAdaptFilt**: Filters adaptors specific to PacBio reads.

**hifiAdaptFilt**: Filters adaptors specific to PacBio hifi reads.

**removeDuplicateReads**: Removes duplicate reads from the dataset.

**filteredReads**: Filters reads longer than 5000 bases.

**seqkitFiltered**: Generates basic statistics for filtered reads.

## Python scripts

1. fetchall.py
2. pipeline.py

## fetchall.py

The **fetchall.py** script is responsible for automating the downloading and SRA data. It parses command-line arguments to specify the downloading method, run type, input data, and output location. The argements it uses are:

* **-t, --type**: Specifies the method to download sequence files. Options include "sra" for SRA Toolkit, "kingfisher" for Kingfisher utility, and "wget" for wget command.
* **-i, --input**: In pipeline mode the input needs to be a URL link from the ENA/SRA database. In manual mode it needs to be a file containing URLs from the ENA/SRA database.
* **-o, --output**: In pipeline mode the output is a direct output path. In manual it's a output directory where the downloaded and processed files will be saved.
* **-r, --run-type**: Choose a run type for usage. Options include pipeline and manual.

## Usage for manual mode

    python script/fetchall.py -t wget -r manual -i input/input.txt -o downloads

## Usage for pipeline mode with snakemake

    python scripts/fetchall.py -t {download_type} -r pipeline -i {file} -o {output}

## pipeline.py

The **pipeline.py** script offers an interactive interface for downloading sequence files for the main pipeline. The script employs the **questionary** library to create user-friendly prompts, allowing the user to select options easily.
