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
    - [Usage](#usage-1)

## Authors

- [@Jesse mittertreiner](https://github.com/AntiCakejesCult)
- [@Giang Le](https://github.com/GiangLeN)

## Introduction

For a brief overview of the main working of the VDJ-AAAP pipeline, please read the [abstract](../README.md#abstract). In this document we will outline the main working of the different created script and pipeline. We discuss parts of the scrips and pipeline that can be changed but not with the [config file](../README.md#configuration-settings). Also the main packages and modules used for the envs are shown with the main package.

## Packages

This pipeline and it coresponding scripts uses different environments to generate the different results and itermidate results. If look closer in to environment files you see that packages are duplicated in different environment files. The reason for this is that you can not inherit packages from enviroments and snakemakes **conda** argument does not allow for multiple environment files.

I will list the different used packages and modules and which version i used. Changing the version of the packaged can be done in the individual envs.

| Package             | Version   |
| ------------------- | --------- |
| **beautifulsoup4**  | 4.12.3    |
| **inspector**       | 1.2       |
| **matplotlib**      | 3.8.4     |
| **python**          | 3.11      |
| **python**          | 3.6       |
| **requests**        | 2.31.0    |
| **bedtools**        | 2.31.0    |
| **biopython**       | 1.83      |
| **blast**           | 2.14.1    |
| **bowtie**          | 1.3.1     |
| **bowtie2**         | 2.5.1     |
| **busco**           | 5.5.0     |
| **hifiasm**         | 0.19.8    |
| **imagemagick**     | 7.0.11_12 |
| **jq**              | 1.5       |
| **meme**            | 4.11.2    |
| **minimap2**        | 2.26      |
| **ncbidatasetscli** | 15.25.0   |
| **openpyxl**        | 3.1.2     |
| **pandas**          | 2.2.1     |
| **pyyaml**          | 6.0.1     |
| **quast**           | 5.2.0     |
| **samtools**        | 1.6       |
| **seqkit**          | 2.8.0     |
| **seqtk**           | 1.4       |
| **unzip**           | 6.0       |

## Pipeline
For the pipeline, just like mentioned, most things can be change with the config file. But here are some examples of things that be changed in the **snakefile** itself. The main [config file](../Snakefile4#L6) can be changed. If changed please change it aswell in the following python scripts [region](../scripts/region.py), [mapping](../scripts/mapping.py), [RSS](../scripts/RSS.py) and [write_annotation_report](../scripts/write_annotation_report.py). You can add extra values in the config file, but make sure you add them in the [snakefile](../Snakefile4#L31) and others [scripts](../scripts/RSS.py#L125). The input directory can also changed if needed please rename every downloads in something else you like. Every command in the snakefile can also be altered to your liking. Please follow the requirements of the package. For the inhouse IMGT scrape and annotation tools please look at the **help rule (-h)**.

## Custom scripts and tools
This pipeline uses a combination of different custom created python scripts and tools to produce results.

## region.py
The script extracts of specific regions from FASTA files using flanking genes defined in the config file. That is parses the alignment data from SAM files to identify the best mapping coordinates for these regions based on the bitwise flag. The final output is a set of region FASTA files, each containing the sequence of a specified region. This process is managed within a Snakemake. If wanted you could change the [snakemake.wildcards](../scripts/region.py#L165) to somethinge that fethces the sample name from the file itself. Lastly you could change [snakemake.input[0]](../scripts/region.py#L191) to the name of the SAM file, then this script could also be used without the pipeline.

## **IMGT Scraper**
This script/tool automates the retrieval of immunoglobulin (IG) and T-cell receptor (TR) VDJ segment .sequences from the IMGT database for specified species. It allows users to customize the sequence fetching process based on several parameters including species type, receptor type (IG or TR), and sequence frame, which is needed to change the type segments to fetch. The script outputs the sequences in FASTA format per region and segment. Optionally it can compile them into a comprehensive library and remove the excess files.

### Usage

```bash 
python imgt_scrape.py -S "Homo sapiens" -T IG --output /path/to/output --create-library --cleanup
```
For more information about this tool, please see the [README.md](https://github.com/BPRC-CGR/IMGT_scrape/tree/development)

## **Annotation**
This is the main tool or set of scripts that is used the generate annotation to discover novel and known VDJ segments present it the input data. It takes list of regions and a library with the help with config file it generates an scala of [annotation result](../README.md#annotation). This tool can also be used with the pipeline. By default it uses minimap2, bowtie and bowtie2 to do annotation.  

### annotation.py
This is the main script that is called then running the annotation tool. It start by calling the mapping script to do the initial mapping with the different mapping tools that are specified. Than it evaluate the different found segment with the help of BLAST, to find the differences between the found segments and the used library. The command that is used to run BLAST. Than it writes the initial excel files. Lastly it validates the found segments with the help of the RSS. 

### mapping.py

### Usage 
```bash 
python vdj_annotation.py --input /path/to/sequence_data --library /path/to/reference_library.fasta --output /path/to/output_directory
```
