import os
import pandas as pd
from scripts.pipeline import get_ids

ACCESSION, = glob_wildcards("downloads/{accs}.fastq.gz")

id_dict = get_ids()
print(id_dict)

rule all:
    input:
        expand("filtered/{accession}.filtered.fastq.gz", accession=id_dict.keys()),
#        expand("seqkit_filtered/filterd_{accession}_stat.tsv", accession=id_dict.keys())


# donwload rule to download data from SRA db with kingfisher
rule SRA_download:
    output:
        "downloads/{accession}.fastq.gz"
    params:
        lambda wildcards: id_dict[wildcards.accession]
    threads:
        8
    conda:
        "envs/kingfisher.yaml"
    shell:
        """
        kingfisher get -r {wildcards.accession} -m ena-ftp
        mv {params}.fastq.gz {output}
        """

# QC rule for getting statistics for unfilterd reads with seqkit stats.
rule seqkit:
    input:
        ancient("downloads/{accession}.fastq.gz")
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
        pbfilt=temp("downloads/{accession}/{accession}.filt.fastq.gz"),
        hififilt = temp("downloads/{accession}/hifi/{accession}.filt.filt.fastq.gz"),
        filt = "filtered/{accession}.filtered.fastq.gz"
    params:
        input_dir = "downloads",
        pb_dir = "{accession}",
        hifi_dir = "hifi"
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.input_dir}
        bash pbadapterfilt.sh -o {params.pb_dir} -p {wildcards.accession}
        cd {params.hifi_dir}
        bash pbadapterfilt.sh -o {params.hifi_dir} -p {wildcards.accession}
        cp {output.hififilt} {output.filt} 
        """
'''
# Filter rule to remove duplicates with seqkit rmdup.
rule removeDuplicateReads:
    input:
        "filtered/{accession}.filtered.fastq.gz"
    output:
        "downloads/duplicate_free/no_duplicate_{accession}.fastq.gz"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit rmdup {input} -s -i -o {output} 
        """

# Filter rule to remove sequences over 5000 bps with seqkit seq.
rule filteredReads:
    input:
        "downloads/duplicate_free/no_duplicate_{accession}.fastq.gz"
    output:
        "downloads/filtered/filtered_{accession}.fastq.gz"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit seq {input} -m 5000 -o {output} 
        """

# QC rule for getting statistics for filterd reads with seqkit stats.
rule seqkitFiltered:
    input:
        "downloads/filtered/filtered_{accession}.fastq.gz"
    output:
        "seqkit_filtered/filterd_{accession}_stat.tsv"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stats {input} -a -o {output}
        """
'''