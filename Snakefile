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
        expand("downloads/{accession}/alignments/extracted_{accession}.bam", accession=ids.keys()),
        # expand("seqkit/raw_read_{accession}_stat.tsv", accession=ids.keys()),

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
        python scripts/fetchall.py -t {params.download_type} -r pipeline -i {params.fastq_file} -o {output.sra}
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
    log:
        "logs/adaptor/log_{accession}_pb.log"
    benchmark:
        "benchmarks/adaptfilt/benchmark_{accession}_pb.txt"
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
    log:
        "logs/adaptor/log_{accession}_hifi.log"
    benchmark:
        "benchmarks/adaptfilt/benchmark_{accession}_hifi.txt"
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
    log:
        "logs/seqkit/duplicates/log_{accession}_duplicates.log"
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
    log:
        "logs/seqkit/filterd/log_{accession}_filterd.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit seq {input} -m 5000 -o {output} 
        """

# QC rule for getting statistics for filterd reads with seqkit stats.
rule seqkitFiltered:
    input:
        ancient("downloads/{accession}/cleaned/filtered_{accession}.fastq.gz")
    output:
        "seqkit_filtered/filtered_{accession}_stat.tsv"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stats {input} -a -o {output} 2> {log}
        """

rule downloadMmul10:
    output:
        ref = "downloads/mmul10.fna",
        ref_report = "downloads/reports/mmul10_assembly_report.txt",
    params:
        zipped_ref = "downloads/mmul10.fna.gz"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/339/765/GCA_003339765.3_Mmul_10/GCA_003339765.3_Mmul_10_genomic.fna.gz -O {params.zipped_ref}
        gunzip {params.zipped_ref}
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/339/765/GCA_003339765.3_Mmul_10/GCA_003339765.3_Mmul_10_assembly_report.txt -P reports -O {output.ref_report}
        """


rule minimap2:
    input:
        read = "downloads/{accession}/cleaned/filtered_{accession}.fastq.gz",
        mmul10 = "downloads/mmul10.fna",
    output:
        "downloads/{accession}/alignments/{accession}.sam"
    log:
        "logs/minimap2/log_{accession}_alignment.log"
    benchmark:
        "benchmarks/minimap2/benchmark_{accession}_alignment.txt"
    threads:
        20
    params:
        read_type = "map-hifi"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax {params.read_type} -t {threads} {input.mmul10} {input.read} > {output} 2> {log}
        """

rule extractMappedReads:
    input:
        "downloads/{accession}/alignments/{accession}.sam"
    output:
        "downloads/{accession}/alignments/extracted_{accession}.bam",
    log:
        "logs/samtools/log_{accession}_alignment.log",
    conda:
        "envs/samtools.yaml"
    threads:
        10
    shell:
        """
        samtools view -@ {threads} -b -F4 {input} > {output} 2> {log}
        """

# rule sortIndexBam:
#     input:
#         "downloads/{accession}/alignments/extracted_{accession}.bam"
#     output:

    
#     log:
#         "logs/samtools/log_chr_{accession}_alignment.log"
#     conda:
#         "envs/samtools.yaml"
#     threads:
#         10
#     #     extracted_chr = "downloads/{accession}/alignments/extracted_{chrs}_{accession}.bam",
#     #     index_bam = "downloads/{accession}/alignments/sorted_{accession}.bam.bai",
#     #     sorted_bam = "downloads/{accession}/alignments/sorted_{accession}.bam",
#     shell:
#     """
#     samtools sort -o {params.sorted_bam} {params.bam}
#     samtools index {params.sorted_bam}
#     samtool view -@ {threads} {params.bam} {wildcards.chrs}> {output.extracted_chr} 2> {log.chr_log}
#     """

# rule longqc:
#     input:
#         "downloads/{accession}/cleaned/filtered_{accession}.fastq.gz"
#     output:
#         "downloads/{accession}/longQC_results"
#     singularity:
#         "docker://quay.io/biocontainers/longqc"
#     shell:
#         """
#         LongQC sampleqc -x pb-sequel -o {output} -p 4 {input}
#         """
