# VDJ-Insights

## Introduction

VDJ-Insights is a robust software package for accurate annotation of the V, D, and J gene segments within immunoglobulin (IG) and T-cell receptor (TCR) genomic regions. In addition to segment annotation, it evaluates gene functionality, detects recombination signal sequences (RSS), and annotates complementary-determining regions 1 and 2 (CDR1 and CDR2). These features extend the utility of VDJ-Insights beyond gene annotation, providing a powerful framework for functional immunogenetics and enabling evolutionary and comparative analyses at individual, population, and species levels.

---

## Installation

VDJ-Insights is currently only supported on Linux systems. Before running the pipeline, please ensure that Python (version 3.7 or higher) and Conda are installed on your system. 
You can install VDJ-Insights using one of the following methods:

### Option 1: Clone the repository
1. Clone the VDJ-Insights repository:
   ```bash
   git clone https://github.com/BPRC-CGR/VDJ-insights
   ```

2. Navigate to the repository directory:
   ```bash
   cd vdj_insights
   ```

3. Run the pipeline using Python's -m option:
   ```bash
   python -m vdj_insights <annotation|html> [arguments]
   ```
**Note:** When cloning the repository, the pipeline must always be executed using the ```python -m``` option. This ensures that Python correctly recognizes the package structure and runs the pipeline without additional installation steps.

### Option 2: Install via pip
1. Use pip to install VDJ-Insights:
   ```bash
   pip install vdj_insights
   ```
2. Run the pipeline:
   ```bash
   vdj_insights <annotation|html> [arguments]
   ```

## Using VDJ-Insights
Use the following command to run the annotation script:

```bash
python vdj-insights annotation -a <assembly_directory> | -i <region_directory> -l <library_directory/library.fasta> -r <receptor_type> -s <species_name> -f <flanking_genes> -t <threads> -m <mappingtool, mapping_tool> -M <metadata_directory> -o <output_directory> --default
```

### **Required Arguments:**
| **Argument**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    | **Description**                                                                                                                                                         | **Example**                                                                    |
|-----------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
| `-r`,<br> `--receptor-type`                      | Type of receptor to analyze. Choices: `IG` (immunoglobulin) or `TR` (T-cell receptor). **Required when using `--default`.**                                            | `-r TR`                                                                        |
| `-i`,<br> `--input` <br><br> or <br><br> `-a`, `--assembly` | Directory containing either extracted sequence regions (`--input`), referring to sequences of the region of interest already isolated from a genome assembly <br><br> **or** <br><Br> genome assembly FASTA files (`--assembly`).                                                         | `-i /path/to/region` <br> `-a /path/to/assembly`                                 |
| `-l`,<br> `--library`                            | Path to the FASTA library file containing reference V(D)J segment sequences.                                                                                                     | `-l /path/to/library.fasta`                                                    |
| `-f`,<br> `--flanking-genes`                     | Comma-separated list of flanking genes provided as key-value pairs in JSON format. If only one flanking gene is present, use `"-"` as a placeholder for the missing side.              | `-f '{"IGH": ["PACS2", "-"], "IGK": ["RPIA", "PAX8"], "IGL": ["GANZ", "TOP3B"]}'` |
| `-s`,<br> `--species`                           | Scientific species name (e.g., `Homo sapiens`).                                                                                                                                   | `-s "Homo sapiens"`                                                            |

---

### **Optional Arguments:**
| **Argument**        | **Description**                                                                                    | **Example**              |
|---------------------|----------------------------------------------------------------------------------------------------|-------------------------|
| `-M, --metadata`    | Path to the metadata file (.xlsx) relevant to the analysis.                                        | `-M metadata.xlsx`       |
| `-o, --output`      | Output directory for the results (default: `annotation_results` in the current directory).         | `-o /path/to/output`     |
| `-m, --mapping-tool`| Available mapping tools: `minimap2`, `bowtie`, `bowtie2`. (Default: all).                          | `-m minimap2`            |
| `-t, --threads`     | Number of threads for parallel processing (default: `8`).                                          | `-t 16`                  |
| `--default`         | Use default settings (cannot be used with `--flanking-genes`).                                     | `--default`              |
| `-S, --scaffolding` | Path to the reference genome (FASTA). **Only supports assemblies representing a single phased contigs.** | `-S /path/to/reference.fasta`|


[Download metadata template](https://github.com/BPRC-Bioinfo/VDJ-insights/blob/main/vdj_insights/metadata/metadata.xlsx)

### Important notes

- If using the `-i/--input` flag, do not specify `-f/--flanking-genes` as they are only needed when using `-a/--assembly`.
- If using the `-i/--input` flag, the name of the file is like this `<sample-name>_<region>.fasta` and the file should be in the directory.
- If using the `--default` flag, do not specify `-f/--flanking-genes` as they are mutually exclusive with `--default`.
- If using the `--default` flag, the annotation tool will automatically download the appropriate V(D)J segment library for the selected receptor `(-r)` type and species `(-s)`. This means you don’t need to specify the flanking genes or provide a local library file.
- If using the `--scaffolding` flag, RagTag scaffolding requires a single-haplotype assembly. If your input contains both haplotypes, split them and scaffold each one independently.

### Example
Download the assembly file (GCA_009914755.4) using the following wget command:

```bash
wget https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/genome/Homo_sapiens-GCA_009914755.4-unmasked.fa.gz
```
(base) jaimy@cgr-jaimy-lnx:~/genome$ gunzip Homo_sapiens-GCA_009914755.4-unmasked.fa.gz 

gunzip Homo_sapiens-GCA_009914755.4-unmasked.fa.gz
Run the annotation tool with the downloaded assembly:

```bash
python -m vdj-insights annotation -a /path/to/GCA_009914755.4-unmasked.fa.gz -r IG -s "Homo sapiens" --default
or
vdj-insights annotation -a /path/to/GCA_009914755.4-unmasked.fa.gz -r IG -s "Homo sapiens" --default
```

When using the `--default` flag, the annotation tool will automatically download the appropriate V(D)J segment library for the selected receptor `(-r)` type and species `(-s)`. This means you don’t need to specify the flanking genes or provide a local library file.

## Annotation results
The results of the V(D)J gene segment identification, generated by the annotation tool, are stored in the **`annotation`** folder. This folder contains the following Excel files: 
- `annotation_report_known.xlsx` contain information about known V(D)J gene segments, including recombination signal sequences. 
- `annotation_report_novel.xlsx` contain information about novel V(D)J gene segments, including recombination signal sequences.
- `annotation_report_all.xlsx` contain known and novel information about known V(D)J gene segments, including recombination signal sequences. 
- `tmp/blast_results.xlsx` contains BLAST search results used for validation.  
- `tmp/report.xlsx` provides a summary of the overall findings from the alignment tools.
  
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
python -m vdj_insights.html -i /path/to/output --show
or
vdj_insights html -i /path/to/output --show
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

