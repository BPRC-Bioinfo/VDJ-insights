# Snakefile

import os

# List all the fastq files in the data directory.
FASTQ_FILES = [fastq_file for fastq_file in os.listdir("data/") if fastq_file.endswith(".fastq.gz")]

# Default rule to process all FASTQ files.
rule all:
    input:
        expand("data/pb/hifi/{sample}.filt.filt.fastq.gz", sample=[fastq_file.replace(".fastq.gz", "") for fastq_file in FASTQ_FILES])


# Rule to run pbadapterfilt.sh on every FASTQ file in the data directory.
rule pbAdaptFilt:
    input: 
        fastq = "data/{sample}.fastq.gz"
    output: 
        "data/pb/{sample}.filt.fastq.gz"
    params:
        data = "data",
        pb = "pb"
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.data}
        bash pbadapterfilt.sh -o {params.pb}
        """

rule hifiAdaptFilt:
    input: 
        "data/pb/{sample}.filt.fastq.gz"
    output: 
        "data/pb/hifi/{sample}.filt.filt.fastq.gz"
    params:
        data = "data/pb",
        hifi = "hifi"
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.data}
        bash pbadapterfilt.sh -o {params.hifi}
        """
