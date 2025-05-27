
[metadata](https://github.com/BPRC-Bioinfo/VDJ-insights/blob/main/vdj_insights/metadata/metadata.xlsx)
![Logo](vdj_insights/images/BPRC_logo.png)
# VDJ-Insights

## Abstract

VDJ-Insights offers a robust framework for assembling and annotating the TCR and BCR regions using long sequence reads from Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT). We designed this tool to uncover both novel and known V(D)J gene segments within T-cell receptors (TCR) or B-cell receptors (BCR). This tool supports analysis across various species, given the availability of a reference genome on the NCBI.

---

## Installation

VDJ-Insights is currently supported only on Linux systems. The Python environment required to run the pipeline is automatically configured using PIO, eliminating the need for manual Conda environment setup. Please ensure that Python and Conda are installed on your system before running the pipeline.
To install VDJ-Insights, you can choose one of the following methods:

### Option 1: Clone the repository
1. Clone the VDJ-Insights repository:
   ```bash
   git clone https://github.com/BPRC-CGR/VDJ-insights
   ```

2. Navigate to the repository directory and run PIO to install dependencies:
   ```bash
   cd vdj_insights
   ```

3. Run the pipeline using Python's -m option:
   ```bash
   python -m vdj_insights.<pipeline|annotation|html> [arguments]
   ```
**Note:** When cloning the repository, the pipeline must always be executed using the python -m option. This ensures that Python correctly recognizes the package structure and runs the pipeline without additional installation steps.
### Option 2: Install via pip
1. Use pip to install VDJ-Insights directly:
   ```bash
   pip install vdj_insights
   ```

After installation, the environment is ready to use, and you can proceed with running the pipeline.

## Overview VDJ-Insights

Within VDJ-Insights, you have the option to assemble and analyze the V(D)J gene segments. You only need to specify the reference genome, either as a NCBI genome code or as a an input fasta file. Specify the type of receptor `TR` or `IG` and the species name of the animal you want to analyze.

## Assembly
The complete pipeline can be executed using the following command:

```bash
python vdj-insights pipeline -ont <nanopore_data.fastq.gz> -pb <pacbio_data.fastq.gz> -ref <reference_genome> -r <receptor_type> -s <species_name> -t <threads> --default
```

### **Required Arguments**
| **Argument**         | **Description**                                                                                                    | **Example**                |
|----------------------|--------------------------------------------------------------------------------------------------------------------|---------------------------|
| `-ont, --nanopore`   | Path to the Oxford Nanopore Technologies (ONT) reads file in `.fastq.gz` format.                                   | `-ont reads_nanopore.fastq.gz` |
| `-pb, --pacbio`      | Path to the Pacific Biosciences (PacBio) reads file in `.fastq.gz` format.                                         | `-pb reads_pacbio.fastq.gz`    |
| `-ref, --reference`  | Path to the reference genome file (`.fasta`, `.fna`) or valid accession code.                                      | `-ref reference.fasta`         |
| `-r, --receptor-type`| Type of receptor to analyze: `TR` (T-cell receptor) or `IG` (Immunoglobulin).                                      | `-r TR`                        |
| `-s, --species`      | Scientific name of the species (e.g., "Homo sapiens").                                                             | `-s "Homo sapiens"`            |

---

### **Optional Arguments**
| **Argument**         | **Description**                                                                                          | **Example**                |
|----------------------|----------------------------------------------------------------------------------------------------------|---------------------------|
| `-f, --flanking-genes` | Comma-separated list of flanking genes provided as key-value pairs in JSON format. If only one flanking gene is available, use `"-"` as a placeholder. | `-f '{"IGH": ["PACS2", "-"], "IGK": ["RPIA", "PAX8"], "IGL": ["GANZ", "TOP3B"]}'`           |
| `-c, --chromosomes`  | List of chromosomes where `TR` or `IG` is located. Values must be integers between 1-22 or `X`, `Y`.     | `-c 14,22`                 |
| `-t, --threads`      | Number of threads for parallel processing (default: `8`).                                                | `-t 16`                    |
| `--default`          | Uses default settings for the analysis. Cannot be used with `-f/--flanking-genes` or `-c/--chromosomes`. | `--default`                |

---

### **Important Notes**

- If using the `--default` flag, do not specify `-f/--flanking-genes` or `-c/--chromosomes`, as these are mutually exclusive with `--default`.
- The `-ref/--reference` flag is required if a reference genome file is available; otherwise, this step can be omitted if the default assembly process is to be used.

## Annotation
It is possible to analyse already assemble genomes using the FASTA file as input. You can use the following command to run the script:

```bash
python vdj-insights annotation -a <assembly_directory> | -i <region_directory> -l <library_directory/library.fasta> -r <receptor_type> -s <species_name> -f <flanking_genes> -t <threads> -m <mappingtool, mapping_tool> -M <metadata_directory> -o <output_directory> --default
```

### **Required Arguments:**
| **Argument**                      | **Description**                                                                                                                                                          | **Example**                                                                       |
|-----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| `-r, --receptor-type`             | Type of receptor to analyze. Choices are `TR` (T-cell receptor) or `IG` (immunoglobulin). **(Required when using `--default`)**                                          | `-r TR`                                                                           |
| `-i, --input` or `-a, --assembly` | Directory containing the extracted sequence regions (`--input`) **or** the assembly FASTA files (`--assembly`).                                                          | `-a /path/to/assembly` or `-i /path/to/region`                                    |
| `-l, --library`                   | Path to the library FASTA file containing V(D)J segment sequences. **(Required when using -i)**                                                                          | `-l /path/to/library.fasta`                                                       |
| `-f, --flanking-genes`            | Comma-separated list of flanking genes provided as key-value pairs in JSON format. If only one flanking gene is available, use `"-"` as a placeholder.                   | `-f '{"IGH": ["PACS2", "-"], "IGK": ["RPIA", "PAX8"], "IGL": ["GANZ", "TOP3B"]}'` |
| `-s, --species`                   | Species name, e.g., `Homo sapiens`.                                                                                                                                      | `-s "Homo sapiens"`                                                               |
| `-S, --scaffolding`               | Path to the reference genome (FASTA) containing the chromosomes of interest for the selected species. **Only supports assemblies representing a single phased contigs.** | `-S /path/to/reference.fasta`                                                     |



---

### **Optional Arguments:**
| **Argument**        | **Description**                                                                                    | **Example**              |
|---------------------|----------------------------------------------------------------------------------------------------|-------------------------|
| `-M, --metadata`    | Path to the metadata file (.xlsx) relevant to the analysis.                                        | `-M metadata.xlsx`       |
| `-o, --output`      | Output directory for the results (default: `annotation_results` in the current directory).         | `-o /path/to/output`     |
| `-m, --mapping-tool` | Available mapping tools: `minimap2`, `bowtie`, `bowtie2`. (Default: all).                          | `-m minimap2`            |
| `-t, --threads`     | Number of threads for parallel processing (default: `8`).                                          | `-t 16`                  |
| `--default`         | Use default settings (cannot be used with `--flanking-genes`).                                     | `--default`              |


### Important notes

- If using the `-i/--input` flag, do not specify `-f/--flanking-genes` as they are only needed when using `-a/--assembly`.
- If using the `-i/--input` flag, the name of the file is like this `sample-name_region.fasta` and the file should be in the directory.
- If using the `--default` flag, do not specify `-f/--flanking-genes` as they are mutually exclusive with `--default`.
- If using the `--default` flag, the annotation tool will automatically download the appropriate V(D)J segment library for the selected receptor `(-r)` type and species `(-s)`. This means you don’t need to specify the flanking genes or provide a local library file.
- If using the `--scaffolding` flag, RagTag scaffolding requires a single-haplotype assembly. If your input contains both haplotypes, split them and scaffold each one independently.

### Example
Download the assembly file (GCA_009914755.4) using the following wget command:

```bash
wget https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/genome/Homo_sapiens-GCA_009914755.4-unmasked.fa.gz
```

Run the annotation tool with the downloaded assembly:

```bash
python vdj-insights annotation -a /path/to/GCA_009914755.4-unmasked.fa.gz -r IG -s "Homo sapiens" --default
```

When using the `--default` flag, the annotation tool will automatically download the appropriate V(D)J segment library for the selected receptor `(-r)` type and species `(-s)`. This means you don’t need to specify the flanking genes or provide a local library file.

## Annotation results
The results of the V(D)J gene segment identification, generated by the annotation tool, are stored in the **`annotation`** folder. This folder contains the following Excel files: 
- `annotation_report_known.xlsx` contain information about known V(D)J gene segments, including recombination signal sequences. 
- `annotation_report_novel.xlsx` contain information about novel V(D)J gene segments, including recombination signal sequences.
- `annotation_report_all.xlsx` contain known and novel information about known V(D)J gene segments, including recombination signal sequences. 
- `blast_results.xlsx` contains BLAST search results used for validation.  
- `report.xlsx` provides a summary of the overall findings from the alignment tools.
  
Each annotation report (known | novel) contains the following columns with detailed information about the identified segments:

| **Column**                     | **Explanation**                                                                                                                                                                                                                               | **Example**                  |
|---------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------|
| **Sample**                      | The name of the sample. | `Sample_001` |
| **Haplotype**                   | The haplotype. | `1` |
| **Region**                      | The region type.  | `IGHV` |
| **Segment**                     | The gene segment. | `V` |
| **Start coord**                 | The start coordinate on the contig. | `12345`|
| **End coord**                   | The end coordinate on the contig. |`12789`|
| **Strand**                      | The orientation of the segment: `+` for 5' to 3' and `-` for 3' to 5'. | `+`|
| **Library name**                   | The name of the closest reference gene segment for the identified segment. | `IGHV3-23*01`|
| **Target name**               | The name assigned to the new segment, derived from the closest reference, with "like" appended to indicate similarity. | `IGHV3-23-like` |
| **Short name**                  | The gene name (IMGT gene name nomenclature). | `IGHV3*01` |
| **Similar references**          | Other reference gene segments with the same start and end coordinates; the best match is chosen based on the mutation count and reference name.                                                                                               | `IGHV3-33*02`                 |
| **Target sequence**           | The DNA sequence of the identified "Old name-like" segment. | `ATGGTGCAAGC...` |
| **Library sequence**               | The DNA sequence of the closest reference gene segment. | `ATGGTGCAAAC...` |
| **Mismatches**                  | The total number of mismatches between the identified segment and the reference.                                                                                                                                                               | `3`                           |
| **% Mismatches of total alignment** | The percentage of mismatches relative to the total length of the alignment between the identified segment and the reference.                                                                                                                   | `1.5%`                        |
| **% identity**                  | The percentage of identical bases between the identified segment and the reference across the entire alignment.  | `98.5%` |
| **BTOP**                        | BLAST traceback string that describes the exact location of substitutions, insertions, and deletions in the alignment.| `10A5G3T` |
| **SNPs**                        | The number of single nucleotide polymorphisms (SNPs) detected in the alignment. | `2` |
| **Insertions**                  | The number of insertions in the identified segment relative to the reference. | `1` |
| **Deletions**                   | The number of deletions in the identified segment relative to the reference. | `0` |
| **Mapping tool**                        | The name of the mapping tool used for the annotation. | `Minimap2` |
| **Function**                    | The functional classification of the segment: "F/ORF" for functional/open reading frame, "P" for potentially functional/open reading frame, or "pseudogene" if an early stop codon is detected.                                               | `F/ORF`                       |
| **Status**                      | Indicates whether the segment is classified as **Known** or **Novel**. | `Novel`|
| **Message**                     | A generated message for the segment if stop codons are detected in critical positions. | `the STOP-CODON at the 3' end of the V-REGION can be deleted by rearrangement`  |
| **Population**                  | The population group associated with the sample, if metedata is provided. | `Dutch` |                                                  

## Web interface report 
The pipeline provides an interactive web interface for visualizing and exploring the annotation results. You can generate and open the web-based Flask report using the following command:

```bash
python vdj_insights html -i /path/to/output --show
```

## Citing VDJ-Insights
If you use VDJ-Insights in your work, please cite:
<cite>

## Acknowledgements
This tool was developed by the department of Comparative genetics & Refinement of the Biomedical Primate Research Centre ([BPRC](https://www.bprc.nl/en)) in Rijswijk, the Netherlands.

- [@Jesse mittertreiner](https://github.com/AntiCakejesCult)
- [@Sayed Jamiel Mohammadi](https://github.com/sayedjm)
- [@Giang Le](https://github.com/GiangLeN)
- [@SusanOtt](https://github.com/SusanOtt)
- [@Jesse Bruijnesteijn](https://github.com/JesseBNL)

