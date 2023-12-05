#!bin/bash
flye --nano-hq EAW_pacbio_chr7.fastq.gz EAW_nanopore_chr7.fastq.gz --iterations 0 --out-dir flye_giang --threads 24 2> flye.log
flye --pacbio-hifi EAW_pacbio_chr7.fastq.gz --resume-from polishing --out-dir flye_giang --threads 24 2> flye.log