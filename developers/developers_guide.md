![Logo](../images/BPRC_logo.png)
# Developers guide VDJ-AAAP pipeline

## Table of Contents
- [Developers guide VDJ-AAAP pipeline](#developers-guide-vdj-aaap-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Authors](#authors)
  - [Introduction](#introduction)
  - [Packages](#packages)

## Authors

- [@Jesse mittertreiner](https://github.com/AntiCakejesCult)
- [@Giang Le](https://github.com/GiangLeN)

## Introduction

For a brief overview of the main working of the VDJ-AAAP pipeline, please read the [abstract](../README.md#abstract). In this document we will outline the main working of the different created script and pipeline. We discuss parts of the scrips and pipeline that can be changed but not with the [config file](../README.md#configuration-settings). Also the main packages and modules used for the envs are shown with the main package.

## Packages

This pipeline and it coresponding scripts uses different environments to generate the different results and itermidate results. If look closer in to environment files you see that packages are duplicated in different environment files. The reason for this is that you can not inherit packages from enviroments and snakemakes **conda** argument does not allow for multiple environment files.

I will list the different used packages and modules and which version i used. Changing the version of the packaged can be done in the individual envs.

| Package         | Version      |
|-----------------|--------------|
| **beautifulsoup4**  |  4.12.3            |
| **inspector**       | 1.2          |
| **matplotlib**      | 3.8.4             |
| **python**          | 3.11         |
| **python**          | 3.6          |
| **requests**        | 2.31.0             |
| **bedtools**        | 2.31.0             |
| **biopython**       | 1.83             |
| **blast**           | 2.14.1             |
| **bowtie**          | 1.3.1             |
| **bowtie2**         | 2.5.1             |
| **busco**           | 5.5.0             |
| **hifiasm**         | 0.19.8             |
| **imagemagick**     | 7.0.11_12             |
| **jq**              | 1.5            |
| **meme**            | 4.11.2             |
| **minimap2**        | 2.26             |
| **ncbidatasetscli** | 15.25.0             |
| **openpyxl**        | 3.1.2             |
| **pandas**          | 2.2.1             |
| **pyyaml**          | 6.0.1             |
| **quast**           | 5.2.0             |
| **samtools**        | 1.6             |
| **seqkit**          | 2.8.0             |
| **seqtk**           | 1.4             |
| **unzip**           | 6.0             |
