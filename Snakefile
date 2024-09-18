import os
from pathlib import Path
import pandas as pd


ACCESSION, MACHINE = glob_wildcards("downloads/{accession}_{machine}.fastq.gz")
ACCESSION, MACHINE = list(set(ACCESSION)), list(set(MACHINE))

wildcard_constraints:
    accession = "|".join(ACCESSION),
    machine = "|".join(MACHINE),

# Target rule specifying the desired final output
rule all:
    input:
        expand("final/all_outputs_{accession}_{machine}_chr{assembly_chrs}_hap{haplo}.txt", accession=ACCESSION, machine=MACHINE, assembly_chrs=config["ASSEMBLY_CHROMOSOMES"], haplo=config["HAPLOTYPES"]),


# Checkpoint for dynamically splitting the FASTQ file
checkpoint split_fastq:
    input:
        fastq="downloads/{accession}_{machine}.fastq.gz"
    output:
        temp(directory("split_files/{accession}_{machine}"))
    benchmark: 
        "benchmarks/split_fastq_{accession}_{machine}.txt"
    log:
        "logs/split_fastq_{accession}_{machine}.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        # your shell commands for splitting the file
        mkdir -p {output}

        #run seqkit split2       
        seqkit split2 -s 1000000 -O {output} {input.fastq} 2> {log}
        """


def getSplitFastqFiles(wildcards):
    checkpoint_output = checkpoints.split_fastq.get(accession=wildcards.accession, machine=wildcards.machine).output[0]
    parts = glob_wildcards(os.path.join(checkpoint_output, "{accession}_{machine}.part_{i}.fastq.gz")).i
    expanded_paths = expand("QC/raw/{accession}_{machine}.part_{i}.stats",
                            accession=wildcards.accession,
                            machine=wildcards.machine,
                            i=parts)
    return expanded_paths


# Rule for initial QC
rule seqkit:
    input:
        getSplitFastqFiles
    output:
        "results_{accession}.txt",
    benchmark: 
        "benchmarks/seqkit_{accession}.txt"
    log:
        "logs/seqkit_{accession}.log"
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule rawStats:
    input:
        "downloads/{accession}_{machine}.fastq.gz"
    output:
        "QC/raw/{accession}_{machine}.stats"
    benchmark: 
        "benchmarks/rawStats_{accession}_{machine}.txt"
    log:
        "logs/rawStats_{accession}_{machine}.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit stat {input} > {output} 2> {log}
        """

## Remove duplicate reads and filter reads
rule removeDuplicateReads:
    input:
        ancient("split_files/{accession}_{machine}/{accession}_{machine}.part_{i}.fastq.gz")
    output:
        temp("filtered/no_duplicate_{accession}_{machine}.part_{i}.fastq.gz")
    benchmark: 
        "benchmarks/removeDuplicateReads_{accession}_{machine}_part_{i}.txt"
    log:
        "logs/removeDuplicateReads_{accession}_{machine}_part_{i}.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit rmdup {input} -s -i -o {output} 2> {log}
        """


rule filteredReads:
    input:
        ancient("filtered/no_duplicate_{accession}_{machine}.part_{i}.fastq.gz")
    output:
        temp("filtered/filtered_{accession}_{machine}.part_{i}.fastq.gz")
    benchmark: 
        "benchmarks/filteredReads_{accession}_{machine}_part_{i}.txt"
    log:
        "logs/filteredReads_{accession}_{machine}_part_{i}.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit seq {input} -m 5000 -o {output} 2> {log}
        """

def getProcessedStatsFiles(wildcards):
    checkpoint_output = checkpoints.split_fastq.get(accession=wildcards.accession, machine=wildcards.machine).output[0]
    parts = glob_wildcards(f"{checkpoint_output}/{{accession}}_{{machine}}.part_{{i}}.fastq.gz").i
    expanded_paths = expand("QC/filtered/filtered_{accession}_{machine}.part_{i}.stats",
                            accession=wildcards.accession,
                            machine=wildcards.machine,
                            i=parts)
    return expanded_paths


rule processedStats:
    input:
        getProcessedStatsFiles
    output:
        "processed_{accession}.txt"
    benchmark: 
        "benchmarks/processedStats_{accession}.txt"
    log:
        "logs/processedStats_{accession}.log"
    shell:
        """
        cat {input} > {output} 2> {log}
        """


## Combine splitted files
def all_fastq_files(wildcards):
    checkpoint_output = checkpoints.split_fastq.get(accession=wildcards.accession, machine=wildcards.machine).output[0]
    expanded_paths = expand("filtered/filtered_{accession}_{machine}.part_{i}.fastq.gz",
                            accession=wildcards.accession,
                            machine=wildcards.machine,
                            i=glob_wildcards(os.path.join(checkpoint_output, "{accession}.part_{i}.fastq.gz")).i)
    return expanded_paths



rule combineFastQ:
    input:
        ancient(lambda wildcards: all_fastq_files(wildcards))
    output:
        temp("combined/{accession}_{machine}.combined.fastq.gz")
    benchmark: 
        "benchmarks/combineFastQ_{accession}_{machine}.txt"
    log:
        "logs/combineFastQ_{accession}_{machine}.log"
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule seqkitFiltered:
    input:
        ancient("combined/{accession}_{machine}.combined.fastq.gz")
    output:
        "QC/filtered/filtered_{accession}_{machine}.stats"
    benchmark: 
        "benchmarks/seqkitFiltered_{accession}_{machine}.txt"
    log:
        "logs/seqkitFiltered_{accession}_{machine}.log"
    conda:
        "envs/seqkit.yaml"
    threads:
        10
    shell:
        """
        seqkit stats {input} -a -j {threads} -o {output} 2> {log}
        """

## Mapping, sorting, indexing
rule minimap2:
    input:
        read = ancient("combined/{accession}_{machine}.combined.fastq.gz"),
        reference = ancient("reference/genome/new_reference.fasta"),
    output:
        temp("alignments/{accession}_{machine}.sam")
    benchmark:
        "benchmarks/minimap2_{accession}_{machine}_alignment.txt"
    log:
        "logs/minimap2_{accession}_{machine}_alignment.log"
    threads: 
        24
    params:
        read_type_pacbio = "map-hifi",
        read_type_nanopore = "map-ont"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        if [[ "{wildcards.machine}" == "pacbio" ]]; then
            minimap2 -ax {params.read_type_pacbio} -t {threads} {input.reference} {input.read} > {output} 2> {log}
        elif [[ "{wildcards.machine}" == "nanopore" ]]; then
            minimap2 -ax {params.read_type_nanopore} -t {threads} {input.reference} {input.read} > {output} 2> {log}
        fi
        """

rule extractMappedReads:
    input:
        ancient("alignments/{accession}_{machine}.sam")
    output:
        temp("alignments/extracted_{accession}_{machine}.bam")
    benchmark: 
        "benchmarks/extractMappedReads_{accession}_{machine}.txt"
    log:
        "logs/extractMappedReads_{accession}_{machine}.log"
    conda:
        "envs/samtools.yaml"
    threads:
        24
    shell:
        """
        samtools view -@ {threads} -bh {input} > {output} 2> {log}
        """

rule sortIndexBam:
    input:
        ancient("alignments/extracted_{accession}_{machine}.bam")
    output:
        sorted_bam = temp("alignments/sorted_{accession}_{machine}.bam"),
        index_bam = temp("alignments/sorted_{accession}_{machine}.bam.bai"),
    benchmark: 
        "benchmarks/sortIndexBam_{accession}_{machine}.txt"
    log:
        "logs/sortIndexBam_{accession}_{machine}.log"
    conda:
        "envs/samtools.yaml"
    threads:
        10
    shell:
        """
        samtools sort -@ {threads} -o {output.sorted_bam} {input}
        samtools index -@ {threads} {output.sorted_bam} 2> {log}
        """


checkpoint separateChrs:
    input:
        bam = ancient("alignments/sorted_{accession}_{machine}.bam"),
        index_bam = ancient("alignments/sorted_{accession}_{machine}.bam.bai"),
        reference = ancient("reference/chromosomes/reference_chr{all_chrs}.txt")
    output:
        directory("chromosomes/{accession}_{machine}_{all_chrs}")
    benchmark: 
        "benchmarks/separateChrs_{accession}_{machine}_{all_chrs}.txt"
    log:
        "logs/separateChrs_{accession}_{machine}_{all_chrs}.log"
    conda:
        "envs/samtools.yaml"
    threads: 
        8
    shell:
        """
        mkdir -p {output}
        for i in $(cat {input.reference}); do
            # Extract individual chrs
            samtools view -b {input.bam} $i > {output}/$i".bam"
        done
        """

def combineChrFragments(wildcards):
    outdir = checkpoints.separateChrs.get(**wildcards).output[0]
    bams = glob_wildcards(os.path.join(outdir, "{bam}.bam")).bam
    return expand(os.path.join(outdir, "{bam}.bam"),
        bam = bams)

rule rejoinChrs:
    input:
        combineChrFragments
    output:
        "merged/{accession}_{machine}_chr{all_chrs}.bam"
    benchmark: 
        "benchmarks/rejoinChrs_{accession}_{machine}_{all_chrs}.txt"
    log:
        "logs/rejoinChrs_{accession}_{machine}_{all_chrs}.log"
    conda:
        "envs/samtools.yaml"
    threads: 
        8
    shell:
        """
        ulimit -Sn 4096
        samtools merge -@ {threads} {output} {input} 2> {log}
        """



rule chrFastq:
    input:
        bam = ancient("merged/{accession}_{machine}_chr{all_chrs}.bam"),
    output:
        "converted/chr{all_chrs}_{accession}_{machine}.fastq.gz",
    benchmark: 
        "benchmarks/chrFastq_{accession}_{machine}_{all_chrs}.txt"
    log:
        "logs/chrFastq_{accession}_{machine}_{all_chrs}.log"
    conda:
        "envs/samtools.yaml"
    threads: 
        8
    shell:
        """
        # Extract fastq
        samtools fastq -@ {threads} -0 {output} {input.bam} 2> {log}
        """

rule allChromosomes:
    input:
        expand("merged/{accession}_{machine}_chr{all_chrs}.bam", accession=ACCESSION, machine=MACHINE, all_chrs=config["ALL_CHROMOSOMES"])
    output:
        temp("final/allchromosomes.txt")
    benchmark:
        "benchmarks/allChromosomes.txt"
    log:
        "logs/allChromosomes.log"
    shell:
        """
        echo {input} | tr " " "\n" > {output}
        """



## Hifiasm ultra long assembly
rule hifiasmUlChromosomeAssembly:
    input:
        "final/allchromosomes.txt",
        pacbio = ancient("converted/chr{assembly_chrs}_{accession}_pacbio.fastq.gz"),
        nanopore = ancient("converted/chr{assembly_chrs}_{accession}_nanopore.fastq.gz") 
    output:
        "assembly/hifiasm/chr{assembly_chrs}/{accession}_chr{assembly_chrs}_hifiasmUL.bp.hap1.p_ctg.gfa", 
        "assembly/hifiasm/chr{assembly_chrs}/{accession}_chr{assembly_chrs}_hifiasmUL.bp.hap2.p_ctg.gfa", 
    benchmark: 
        "benchmarks/hifiasmUlChromosomeAssembly_{accession}_{assembly_chrs}.txt"
    log:
        "logs/hifiasmUlChromosomeAssembly_{accession}_{assembly_chrs}.log"
    conda:
        "envs/hifiasm.yaml"
    threads:
        24
    shell:
        """
        hifiasm -o assembly/hifiasm/chr{wildcards.assembly_chrs}/{wildcards.accession}_chr{wildcards.assembly_chrs}_hifiasmUL -t {threads} --ul {input.nanopore}  {input.pacbio} 2> {log}
        """


rule gfaToFasta:
    input:
        ancient("assembly/hifiasm/chr{assembly_chrs}/{accession}_chr{assembly_chrs}_hifiasmUL.bp.hap{haplo}.p_ctg.gfa")
    output:
        "converted/gfatofasta/chr{assembly_chrs}_{accession}_hap{haplo}.fasta"
    benchmark:
        "benchmarks/gfaToFasta_{accession}_{assembly_chrs}_{haplo}.txt"
    log:
        "logs/gfaToFasta_{accession}_{assembly_chrs}_{haplo}.log"
    conda:
        "envs/gfatools.yaml"
    shell:
        """
        gfatools gfa2fa {input} > {output} 2> {log}      
        """

rule allAssemblies:
    input:
        ancient(expand("converted/gfatofasta/chr{assembly_chrs}_{accession}_hap{haplo}.fasta", assembly_chrs=config["ASSEMBLY_CHROMOSOMES"], accession=ACCESSION, haplo=config["HAPLOTYPES"]))
    output:
        "converted/gfatofasta/check.txt"
    shell:
        """
        echo {input} | tr " " "\n" > {output} 
        """
# Quast
rule quastAssemblyStatistics:
    input:
        hap1 = ancient("converted/gfatofasta/chr{assembly_chrs}_{accession}_hap1.fasta"),
        hap2 = ancient("converted/gfatofasta/chr{assembly_chrs}_{accession}_hap2.fasta")
    output:
        directory("quast/hifiasm/chr{assembly_chrs}_{accession}")
    benchmark:
        "benchmarks/quastAssemblyStatistics_{accession}_{assembly_chrs}.txt"
    log:
        "logs/quastAssemblyStatistics_{accession}_{assembly_chrs}.log"
    conda:
        "envs/quast.yaml"
    shell:
        """
        quast.py {input.hap1} {input.hap2} -o {output} 2> {log}
        """

# BUSCO
rule BUSCO:
    input:
        ancient("converted/gfatofasta/chr{assembly_chrs}_{accession}_hap{haplo}.fasta")
    output:
        directory("BUSCO/{accession}/chr{assembly_chrs}_{accession}_hap{haplo}")
    benchmark: 
        "benchmarks/BUSCO_{accession}_{assembly_chrs}_{haplo}.txt"
    log:
        "logs/BUSCO_{accession}_{assembly_chrs}_{haplo}.log"
    params:
        dataset = config["BUSCO_DATASET"]
    conda:
        "envs/busco.yaml"
    shell:
        """
        busco -i {input} -l {params.dataset} -o {output} -m genome 2> {log}
        """

# BUSCO reference
rule BUSCOReference:
    input:
        ancient("reference/chromosomes/reference_chr{assembly_chrs}.fasta")
    output:
        directory("BUSCO/reference/chr{assembly_chrs}")
    benchmark: 
        "benchmarks/BUSCOReference_{assembly_chrs}.txt"
    log:
        "logs/BUSCOReference_{assembly_chrs}.log"
    params:
        dataset = config["BUSCO_DATASET"]
    conda:
        "envs/busco.yaml"
    shell:
        """
        busco -i {input} -l {params.dataset} -o {output} -m genome 2> {log}
        """

# Inspector
rule inspector:
    input:
        pacbio = ancient("converted/chr{assembly_chrs}_{accession}_pacbio.fastq.gz"),
        nanopore = ancient("converted/chr{assembly_chrs}_{accession}_nanopore.fastq.gz"), 
        contig = ancient("converted/gfatofasta/chr{assembly_chrs}_{accession}_hap{haplo}.fasta")
    output:
        directory("inspector/chr{assembly_chrs}_{accession}_hap{haplo}")
    benchmark: 
        "benchmarks/inspector_{accession}_{assembly_chrs}_{haplo}.txt"
    log:
        "logs/inspector_{accession}_{assembly_chrs}_{haplo}.log"
    params:
        data_type = "mixed" 
    conda:
        "envs/inspector.yaml"
    threads:
        12
    shell:
        """
        inspector.py -t {threads} -c {input.contig} -r {input.pacbio} {input.nanopore} -o {output} --datatype {params.data_type} 2> {log}
        """


rule getLibrary:
    output:
        "library/library.fasta"
    benchmark: 
        "benchmarks/getLibrary.txt"
    log:
        "logs/getLibrary.log"
    params:
        settings = config["SETTINGS"],
        species = config["SPECIES"]["name"],
        cell_type = config["SPECIES"]["cell"],
    conda:
        "envs/IMGT.yaml"
    shell:
        """
        python {params.settings}/scripts/IMGT_scrape.py -S "{params.species}" -T {params.cell_type} --create-library --cleanup --simple-headers 2> {log}
        """


rule annotation:
    input:
        ancient("library/library.fasta"),
        "converted/gfatofasta/check.txt"
    output:
        "annotation/annotation_report_novel_rss.xlsx"
    benchmark: 
        "benchmarks/annotation.txt"
    log:
        "logs/annotation.log"
    params:
        settings = config["SETTINGS"],
        species = config["SPECIES"]["name"],
        cell_type = config["SPECIES"]["cell"],
        flanking_genes = ",".join(config["FLANKING_GENES"])
    conda:
        "envs/scripts.yaml"
    threads:
        24
    shell:
        """
        python {params.settings}/scripts/annotation.py -a converted/gfatofasta -l library/library.fasta -s "{params.species}" -r {params.cell_type} -f "{params.flanking_genes}" -t {threads} 2> {log}
        """

rule VDJ_display:
    input:
        "annotation/annotation_report_novel_rss.xlsx"
    output:
        directory("VDJ_visualization"),
    benchmark: 
        "benchmarks/vdj_display.txt"
    log:
        "logs/vdj_display.log"
    params:
        settings = config["SETTINGS"]
    conda:
        "envs/display.yaml"
    shell:
        """
        python {params.settings}/scripts/VDJ_display.py -f "annotation/annotation_report_known_rss.xlsx" "annotation/annotation_report_novel_rss.xlsx" -o "VDJ_visualization" -s "combined" 2> {log}
        """

rule fetchAllInput:
    input:
        ancient("annotation/annotation_report_novel_rss.xlsx"),
        ancient("VDJ_visualization"),
        ancient("QC/raw/{accession}_{machine}.stats"),
        ancient("QC/filtered/filtered_{accession}_{machine}.stats"),
        ancient("quast/hifiasm/chr{assembly_chrs}_{accession}"),
        ancient("BUSCO/{accession}/chr{assembly_chrs}_{accession}_hap{haplo}"),
        ancient("BUSCO/reference/chr{assembly_chrs}"),
        ancient("inspector/chr{assembly_chrs}_{accession}_hap{haplo}"),
    output:
        touch("final/all_outputs_{accession}_{machine}_chr{assembly_chrs}_hap{haplo}.txt")
    benchmark: 
        "benchmarks/fetchAllInput_{accession}_{machine}_{assembly_chrs}_{haplo}.txt"
    log:
        "logs/fetchAllInput_{accession}_{machine}_{assembly_chrs}_{haplo}.log"


