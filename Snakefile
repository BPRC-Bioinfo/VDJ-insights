# Snakefile

import os

# List all the fastq files in the data directory.
FASTQ_FILES = [fastq_file for fastq_file in os.listdir("data/") if fastq_file.endswith(".fastq")]

# Default rule to process all FASTQ files.
rule all:
    input:
        expand("results/{sample}/", sample=[fastq_file.replace(".fastq", "") for fastq_file in FASTQ_FILES])


# Rule to run pbadapterfilt.sh on every FASTQ file in the data directory.
rule HiFiAdaptFilt:
    input: 
        fastq = "data/{sample}.fastq"
    output: 
        directory("results/{sample}/")
    # checks if the right directory exist result/name of fastq file.
    shell:
        """
        cd data/
        bash pbadapterfilt.sh -p {wildcards.sample} -o ../results/{wildcards.sample}/
        if [ ! -d "../results/{wildcards.sample}/" ]; then
            echo "Error: Expected output directory was not created by the script."
            exit 1
        fi
        """

