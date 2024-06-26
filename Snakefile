import os
from pathlib import Path
import pandas as pd

# Setting up the global config object. 
configfile: 'config/config.yaml'

# Setup. Loading Accession, machine, flanking.
ACCESSION, MACHINE = glob_wildcards("downloads/{accession}_{machine}.fastq.gz")
ACCESSION, MACHINE = list(set(ACCESSION)), list(set(MACHINE))
FLANKING = [gene for chromosomes in config["FLANKING"].values() for region in chromosomes.values() for gene in region.values() if gene]


wildcard_constraints:
    accession = "|".join(ACCESSION),
    machine = "|".join(MACHINE),

# Target rule specifying the desired final output
rule all:
    input:
        expand("final/all_outputs_{accession}_{machine}_chr{assembly_chrs}_hap{haplo}.txt", accession=ACCESSION, machine=MACHINE, assembly_chrs=config["ASSEMBLY_CHROMOSOMES"], haplo=config["HAPLOTYPES"]),

rule downloadReference:
    output:
        temp("reference/reference.zip"),
        temp("reference/README.md"),
        ref_report = "reference/reports/assembly_report.jsonl",
        ref = "reference/genome/reference.fna",
    benchmark: 
        "benchmarks/downloadReference.txt"
    log:
        "logs/downloadReference.log"
    params:
        reference_code = config["SPECIES"]["genome"],
        prefix = "reference/ncbi_dataset/data"
    conda:
        "envs/NCBI.yaml"
    shell:
        """
        mkdir -p reference/genome
        cd reference
        datasets download genome accession {params.reference_code} --include genome --filename reference.zip
        cd ../
        datasets summary genome accession {params.reference_code} --report sequence --as-json-lines > {output.ref_report}
        unzip reference/reference.zip -d reference/
        mv {params.prefix}/{params.reference_code}/* {output.ref}
        rm -r reference/ncbi*
        """


rule splitChromosomes:
    input:
        fa = "reference/genome/reference.fna",
        report = ancient("reference/reports/assembly_report.jsonl")
    output:
        fa = "reference/chromosomes/reference_chr{all_chrs}.fasta",
        temp = "reference/chromosomes/reference_chr{all_chrs}.txt"
    benchmark: 
        "benchmarks/splitChromosomes_{all_chrs}.txt"
    log:
        "logs/splitChromosomes_{all_chrs}.log"
    conda:
        "envs/split.yaml"
    shell:
        """
        cat {input.report} | egrep '"chr_name":"{wildcards.all_chrs}"' | jq -r ".refseq_accession"  > {output.temp}
        seqkit grep -f {output.temp} {input.fa} | sed "s/>/>chr{wildcards.all_chrs}_/g" > {output.fa}
        """

# Making new reference file.
rule newReferenceFile:
    input:
        expand("reference/chromosomes/reference_chr{all_chrs}.fasta", all_chrs=config["ALL_CHROMOSOMES"])
    output:
        "reference/genome/reference_new_ref.fasta"
    benchmark: 
        "benchmarks/newReferenceFile.txt"
    log:
        "logs/newReferenceFile.log"
    shell:
        """
        cat {input} > {output}
        """

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
        seqkit split2 -s 1000000 -O {output} {input.fastq} 
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
        cat {input} > {output}
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
        seqkit stat {input} > {output}
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
        cat {input} > {output}
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
        cat {input} > {output}
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
        seqkit stats {input} -a -j {threads} -o {output}
        """

## Mapping, sorting, indexing
rule minimap2:
    input:
        read = ancient("combined/{accession}_{machine}.combined.fastq.gz"),
        reference = ancient("reference/genome/reference_new_ref.fasta"),
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
        index_bam = "alignments/sorted_{accession}_{machine}.bam.bai",
        reference = "reference/chromosomes/reference_chr{all_chrs}.txt"
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
            samtools view -b {input.bam} "chr"{wildcards.all_chrs}"_"$i > {output}/$i".bam"
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
        "{accession}_{machine}_chr{all_chrs}.bam",
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
        samtools merge -@ {threads} -o {output} {input} 
        """


rule chrFastq:
    input:
        bam = ancient("{accession}_{machine}_chr{all_chrs}.bam"),
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
        samtools fastq -@ {threads} -0 {output} {input.bam}
        """

rule allChromosomes:
    input:
        expand("{accession}_{machine}_chr{all_chrs}.bam", accession=ACCESSION, machine=MACHINE, all_chrs=config["ALL_CHROMOSOMES"])
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
        gfatools gfa2fa {input} > {output}       
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
        quast.py {input.hap1} {input.hap2} -o {output}
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
        inspector.py -t {threads} -c {input.contig} -r {input.pacbio} {input.nanopore} -o {output} --datatype {params.data_type}
        """


# Flanking genes
rule downloadFlankingGenes:
    output:
        "flanking/{fgene}/ncbi_dataset/data/gene.fna"
    benchmark: 
        "benchmarks/downloadFlankingGenes_{fgene}.txt"
    log:
        "logs/downloadFlankingGenes_{fgene}.log"
    conda:
        "envs/NCBI.yaml"
    params:
        species_name = config["SPECIES"]["name"]
    shell:
        """
        mkdir -p flanking/{wildcards.fgene}
        cd flanking/{wildcards.fgene}
        datasets download gene symbol {wildcards.fgene} --taxon "{params.species_name}" --include gene --filename {wildcards.fgene}.zip
        unzip {wildcards.fgene}.zip
        """



rule combineFlankingGenes:
    input:
        ancient(expand("flanking/{fgene}/ncbi_dataset/data/gene.fna", fgene=FLANKING))
    output:
        "flanking/all_flanking_genes.fna"
    benchmark: 
        "benchmarks/combineFlankingGenes.txt"
    log:
        "logs/combineFlankingGenes.log"
    shell:
        """
        > {output}
        for file in {input}; do
            header=$(head -1 "$file" | egrep "^>" | awk '{{print $2":" $5":"$NF}}' | tr -d '[]')
            sequence=$(cat "$file" | egrep -v "^>")
            echo -e ">$header\n$sequence" >> {output}
        done
        """


rule flankingByChromosomes:
    input:
        flanking = ancient("flanking/all_flanking_genes.fna")
    output:
        "flanking/chr{assembly_chrs}_{accession}_flanking_genes.fna"
    benchmark: 
        "benchmarks/flankingByChromosomes_{assembly_chrs}_{accession}.txt"
    log:
        "logs/flankingByChromosomes_{assembly_chrs}_{accession}.log"
    params:
        chrs_id = "flanking/chr{assembly_chrs}_{accession}_flanking_genes_id.txt" 
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        cat {input} | egrep "{wildcards.assembly_chrs}" | tr -d ">" > {params.chrs_id}
        seqtk subseq {input} {params.chrs_id} > {output}
        """


rule mapFlankingGenes:
    input:
        reads = ancient("converted/gfatofasta/chr{assembly_chrs}_{accession}_hap{haplo}.fasta"),
        flanking_genes = ancient("flanking/chr{assembly_chrs}_{accession}_flanking_genes.fna")
    output:
        sam_file = temp("flank_alignment/chr{assembly_chrs}_{accession}_hap{haplo}_aligned.sam"),
    benchmark: 
        "benchmarks/mapFlankingGenes_{accession}_{assembly_chrs}_{haplo}.txt"
    log:
        "logs/mapFlankingGenes_{accession}_{assembly_chrs}_{haplo}.log"
    conda:
        "envs/minimap2.yaml"
    threads:
        6
    shell:
        """
        minimap2 -ax asm5 --secondary=no -t {threads} {input.reads} {input.flanking_genes} > {output.sam_file}
        """


## Generate regions
rule regionFiles:
    input:
        ancient("flank_alignment/chr{assembly_chrs}_{accession}_hap{haplo}_aligned.sam")
    output:
        temp(touch("flank_alignment/chr{assembly_chrs}_{accession}_hap{haplo}_aligned.txt"))
    benchmark: 
        "benchmarks/regionFiles_{accession}_{assembly_chrs}_{haplo}.txt"
    log:
        "logs/regionFiles_{accession}_{assembly_chrs}_{haplo}.log"
    conda:
        "envs/scripts.yaml"
    script:
        "scripts/region.py" 


rule combineRegion:
    input:
        expand("flank_alignment/chr{assembly_chrs}_{accession}_hap{haplo}_aligned.txt", assembly_chrs=config["ASSEMBLY_CHROMOSOMES"], accession=ACCESSION, haplo=config["HAPLOTYPES"])
    output:
        temp("annotation/all_regions.txt")
    benchmark: 
        "benchmarks/combineRegion.txt"
    log:
        "logs/combineRegion.log"
    shell:
        """
        echo {input} | tr " " "\n" > {output}
        """

rule getLibrary:
    output:
        "library/library.fasta"
    benchmark: 
        "benchmarks/getLibrary.txt"
    log:
        "logs/getLibrary.log"
    params:
        species = config["SPECIES"]["name"],
        cell_type = config["SPECIES"]["cell"],
    conda:
        "envs/IMGT.yaml"
    shell:
        """
        python scripts/IMGT_scrape.py -S "{params.species}" -T {params.cell_type} --create-library --cleanup --simple-headers
        """


rule annotation:
    input:
        ancient("annotation/all_regions.txt"),
        ancient("library/library.fasta")
    output:
        "annotation/annotation_report_plus.xlsx"
    benchmark: 
        "benchmarks/annotation.txt"
    log:
        "logs/annotation.log"
    conda:
        "envs/scripts.yaml"
    shell:
        """
        python scripts/annotation.py -i region/ -l library/library.fasta
        """


rule fetchAllInput:
    input:
        ancient("annotation/annotation_report_plus.xlsx"),
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


