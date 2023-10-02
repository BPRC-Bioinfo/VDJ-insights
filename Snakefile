import os
import pandas as pd
from pipeline import get_ids


id_dict = get_ids()
print(id_dict)

rule all:
    input:
        expand("seqkit_filtered/filterd_{accession}_stat.tsv", accession=id_dict.keys())


# donwload rule to download data from SRA db with kingfisher
rule SRA_download:
    output:
        "downloads/{accession}.fastq.gz"
    wildcard_constraints:
        accession="^(?!pb/).+$"
    params:
        lambda wildcards: id_dict[wildcards.accession]
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
        "downloads/{accession}.fastq.gz" 
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
        temp("downloads/pb/{accession}.filt.fastq.gz")
    wildcard_constraints:
        accession="[^/]+"
    params:
        input_dir = "downloads",
        output_dir = "pb"
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.input_dir}
        bash pbadapterfilt.sh -o {params.output_dir} -p {wildcards.accession}
        """

# Filter rule for removing adaptor for pacbio hifi reads with seqkit.
rule hifiAdaptFilt:
    input: 
        "downloads/pb/{accession}.filt.fastq.gz"
    output: 
        "downloads/pb/hifi/{accession}.filt.filt.fastq.gz"
    params:
        input_dir = "downloads/pb",
        output_dir = "hifi"
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.input_dir}
        bash pbadapterfilt.sh -o {params.output_dir} -p {wildcards.accession}
        """

# Filter rule to remove duplicates with seqkit rmdup.
rule removeDuplicateReads:
    input:
        "downloads/pb/hifi/{accession}.filt.filt.fastq.gz"
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