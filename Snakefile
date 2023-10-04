import os
import pandas as pd
from scripts.pipeline import *

current = os.getcwd()

input_file = fetchall_args_input_file()
ids = get_ids(f"input/{input_file}")

files = [f for f in os.listdir(f"{current}/downloads") if f.startswith("sra_")]
if len(ids.keys()) > len(files):
    sra_download = fetchall_args_sra_download()
else:
    sra_download = "wget"



# Define the top-level rule that depends on the output from other rules
rule all:
    input:
        expand("downloads/{accession}/alignments/{accession}.sam", accession=ids.keys()),
        expand("seqkit/raw_read_{accession}_stat.tsv", accession=ids.keys()),

# Rule to download data from the SRA database
rule SRA_download:
    output:
        sra="downloads/sra_{accession}.fastq.gz"
    params:
        fastq_file = lambda wildcards: ids[wildcards.accession],
        download_type = sra_download
    conda:
        "envs/sra_download.yaml"
    shell:
        """
        python scripts/fetchall.py -t {params.download_type} -i {params.fastq_file} -o {output.sra}
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
        "downloads/sra_{accession}.fastq.gz"
    output: 
        pbfilt = temp("downloads/{accession}/pb_filtered_{accession}.filt.fastq.gz")
    params:
        filt = "downloads/{accession}/sra_{accession}.filt.fastq.gz" ,
        input_dir = "downloads",
    singularity:
        "docker://australianbiocommons/hifiadapterfilt"
    shell:
        """
        cd {params.input_dir}
        bash pbadapterfilt.sh -o {wildcards.accession} -p sra_{wildcards.accession}
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
        cd ../../
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

rule downloadMmul10:
    output:
        "downloads/mmul10.gz/"
    shell:
        """
        wget --header="Accept: application/gzip" "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_003339765.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_003339765.1.gz" -O {output}
        """


rule minimap2:
    input:
        "downloads/{accession}/cleaned/filtered_{accession}.fastq.gz"
    output:
        "downloads/{accession}/alignments/{accession}.sam"
    params:
        mmul10 = "downloads/mmul10.gz"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax {params.mmul10} {input} > {output}
        """



# rule longqc:
#     input:
#         "downloads/{accession}/cleaned/filtered_{accession}.fastq.gz"
#     output:
#         "downloads/{accession}/longQC_results"
#     singularity:
#         "docker://cymbopogon/longqc"
#     shell:
#         """
#         python longQC.py sampleqc -x pb-sequel -o {output} -p 4 {input}
#         """
