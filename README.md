# TCR Macaque Snakemake Pipeline README
## Overview

This Snakemake pipeline is designed for analysing TCR region of the Rhesus Macaque from publicly available data. It includes a set of rules for downloading SRA data, performing initial QC, removing adaptors, removing duplicates and filtering reads based on length.

## Dependencies

To run this pipeline, make sure you have the following software and environments:

1. Snakemake
2. Conda (for managing environments)
3. Singularity (for containerization)

# Pipeline
`Place holder for pipeline image`

# Python scripts
1. fetchall.py
2. pipeline.py

# fetchall.py
The **fetchall.py** script is responsible for automating the downloading and initial preprocessing of sequence data. It parses command-line arguments to specify the downloading method, input data, and output location. The argements it uses are:
* **-t, --type**: Specifies the method to download sequence files. Options include "sra" for SRA Toolkit, "kingfisher" for Kingfisher utility, and "wget" for wget command.

* **-i, --input**: Specifies the input link for the sequence files to be downloaded. For "wget", this needs to be a URL link from the ENA database.

* **-o, --output**: Specifies the output directory where the downloaded and processed files will be saved.

`python fetchall.py -i ... -t ... -o ...`

# pipeline.py
The **pipeline.py** script offers an interactive interface for downloading sequence files for the main pipeline. The script employs the **questionary** library to create user-friendly prompts, allowing the user to select options easily.
