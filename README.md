# Automated RNA-Seq and Differential Expression Pipeline

An end-to-end, production-ready **Nextflow DSL2** pipeline designed to automate raw RNA-Seq data processing through to downstream Differential Expression Analysis (DEG). The workflow processes paired-end reads, maps them to a reference genome, quantifies transcripts, and executes statistical analyses comparing biological conditions (e.g., late vs. early LUAD samples).

To ensure complete reproducibility, every computational step is completely containerized via Docker using community-verified BioContainers.

## 🧬 Pipeline Workflow Overview

    A[Raw FastQ Pairs] --> B(QUALITY: FastQC)
    A --> C(TRIMMING: Trimmomatic)
    D[Reference Genome] --> E(INDEXING: HISAT2)
    C --> F(MAPPING_SORT: HISAT2 & Samtools)
    E --> F
    F --> G(QUANTIFICATION: featureCounts)
    H[GTF Annotation] --> G
    G --> I(MERGING_COUNT: Python/Pandas)
    I --> J(DEG: DESeq2 in R)
    K[Sample Metadata] --> J

The pipeline executes the following sequential steps:
1. **Quality Control (`FastQC`):** Assesses the sequencing quality of raw paired-end reads.
2. **Adapter Trimming (`Trimmomatic`):** Cleans raw reads by removing adapter sequences and filtering out low-quality bases.
3. **Genome Indexing (`HISAT2`):** Generates index files for the reference genome.
4. **Read Alignment & Sorting (`HISAT2` & `Samtools`):** Maps trimmed reads to the indexed genome and produces sorted, indexed BAM files.
5. **Expression Quantification (`featureCounts`):** Quantifies gene expression counts based on the genomic annotations.
6. **Matrix Merging (`Pandas`):** Aggregates individual sample count files into a centralized gene expression matrix.
7. **Differential Gene Expression (`DESeq2`):** Processes the matrix using metadata covariates to isolate differentially expressed genes (DEGs) and outputs publication-quality visualizations.

## 🛠️ Prerequisites & Dependencies
To execute this pipeline, you only need the following two dependencies installed on your system:
* **Nextflow** (DSL2 enabled)
* **Docker Engine** 

All tool dependencies (FastQC, Trimmomatic, HISAT2, Samtools, Subread, Pandas, and DESeq2) are automatically managed, downloaded, and executed inside isolated Docker containers at runtime.

## 🚀 Quick Start

### 1. Clone the Repository
```bash
git clone https://github.com
cd YOUR_REPO_NAME
```

### 2. Configure Your Parameters
Open `main.nf` or create a custom profile to update your execution paths. The pipeline exposes the following parameters:

```nextflow
params.read    = "/path/to/raw_data/SRR*_{1,2}.fastq.gz"   // Paired-end FastQ pattern
params.ref     = "/path/to/ref_genome/hg38.fna"           // Reference genome fasta
params.adapters= "/path/to/adapters/TruSeq3-PE.fa"        // Trimmomatic adapter file
params.annot   = "/path/to/ref_genome/hg38.gtf"           // Reference GTF annotation
params.meta    = "/path/to/sample_metadata.csv"           // Experimental metadata sheet
params.output  = "/path/to/result"                        // Main output folder
```

### 3. Run the Pipeline
Execute the workflow locally using the default Docker configurations defined in your config:

```bash
nextflow run main.nf
```

To resume a failed or modified run from the last successful checkpoint, append the resume flag:
```bash
nextflow run main.nf -resume
```

## 📂 Output Directory Structure

Upon a successful run, results are neatly grouped and saved into the path specified by `params.output`:

```text
results/
├── qc_report/               # FastQC HTML reports and ZIP archives
├── trimming_output/         # Cleaned, paired FastQ files
├── indexing_output/         # Generated HISAT2 genome index files
├── mapping_output/          # Sorted BAM files and index (.bai) files
├── quantification_output/   # Individual sample count files and summaries
├── all_merged_counts/       # Raw counts matrix (raw_counts_matrix.csv)
└── DEGs_output/             # Differential expression tables and plots
    ├── DESeq2_results.csv   # Full statistical result table
    ├── significant_DEGs.csv # Screened list of significant targets
    ├── normalized_counts.csv# Normalized expression data
    ├── DESeq2_object.rds    # Downstream-ready R object
    ├── volcano_plot.png     # Volcano plot visualization
    ├── MA_plot.png          # MA plot visualization
    ├── PCA_plot.png         # Sample PCA plot
    └── heatmap_top50.png    # Heatmap of top 50 highly variable genes
```

## ⚙️ Resource Configuration

The workflow defaults to optimization profiles configured via `nextflow.config`:
* **Global Limit:** Up to `16 CPUs` allocated to the active executor profile.
* **Process Scaling:** High-compute processes like `INDEXING`, `MAPPING_SORT`, `QUANTIFICATION`, and `DEG` scale out utilizing up to `8 CPUs` each to ensure optimal throughput.

## Recommended GitHub Repository Structure:

```text
your-repo-name/
│
├── main.nf               # The workflow pipeline file you provided
├── nextflow.config       # The configuration file you provided
├── README.md             # The markdown documentation provided above
│
└── bin/                  # Nextflow automatically looks here for custom scripts
    ├── merge_counts.py   # Your custom python script using pandas
    └── DESeq2.R          # Your custom R script executing DESeq2
```
