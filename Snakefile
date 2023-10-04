import os
import pandas as pd
from scripts.pipeline import get_ids

ids = get_ids("input/input.txt")

# Define the top-level rule that depends on the output from other rules
rule all:
    input:
        expand("seqkit_filtered/filtered_{accession}_stat.tsv", accession=ids.keys())

# Rule to download data from the SRA database
rule SRA_download:
    output:
        sra="downloads/sra_{accession}.fastq.gz"
    params:
        fastq_file = lambda wildcards: ids[wildcards.accession]
    conda:
        "envs/sra_download.yaml"
    shell:
        """
        python scripts/fetchall.py -t wget -i {params.fastq_file} -o {output.sra}
        """

# Rule for initial QC
rule seqkit:
    input:
        "downloads/sra_{accession}.fastq.gz"
    output:
        "seqkit/raw_read_{accession}_stat.tsv"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stats {input} -a -o {output}
        """

# Filter rule for removing adaptor for pacbio reads with seqkit rmdup.
rule pbAdaptFilt:
    input: 
        "downloads/{accession}.fastq.gz"
    output: 
        pbfilt = temp("downloads/{accession}/pb_filtered_{accession}.filt.fastq.gz")
    params:
        filt = "downloads/{accession}/{accession}.filt.fastq.gz" ,
        input_dir = "downloads",
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.input_dir}
        bash pbadapterfilt.sh -o {wildcards.accession} -p {wildcards.accession}
        cd ../
        mv {params.filt} {output.pbfilt}
        """

# Filter rule for removing adaptor for pacbio hifi reads with seqkit.
rule hifiAdaptFilt:
    input: 
        "downloads/{accession}/pb_filtered_{accession}.filt.fastq.gz"
    output: 
        filt = temp("downloads/{accession}/hifi/hifi_filtered_{accession}.fastq.gz")
    params:
        input_dir = "downloads/{accession}",
        output_dir = "hifi",
        hififilt = "downloads/{accession}/hifi/pb_filtered_{accession}.filt.filt.fastq.gz",
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.input_dir}
        bash hifiadapterfilt.sh -o {params.output_dir} -p pb_filtered_{wildcards.accession}
        cd ../../../
        mv {params.hififilt} {output.filt} 
        """

# Filter rule to remove duplicates with seqkit rmdup.
rule removeDuplicateReads:
    input:
        "downloads/{accession}/hifi/hifi_filtered_{accession}.fastq.gz"
    output:
        "downloads/{accession}/duplicate_free/no_duplicate_{accession}.fastq.gz"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit rmdup {input} -s -i -o {output} 
        """

# Filter rule to remove sequences over 5000 bps with seqkit seq.
rule filteredReads:
    input:
        "downloads/{accession}/duplicate_free/no_duplicate_{accession}.fastq.gz"
    output:
        "downloads/{accession}/cleaned/filtered_{accession}.fastq.gz"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit seq {input} -m 5000 -o {output} 
        """

# QC rule for getting statistics for filterd reads with seqkit stats.
rule seqkitFiltered:
    input:
        "downloads/{accession}/cleaned/filtered_{accession}.fastq.gz"
    output:
        "seqkit_filtered/filtered_{accession}_stat.tsv"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stats {input} -a -o {output}
        """
