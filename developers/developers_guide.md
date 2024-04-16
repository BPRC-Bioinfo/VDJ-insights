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



