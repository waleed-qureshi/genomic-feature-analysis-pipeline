# Genomic Feature Extraction & Analysis Pipeline

A Python-based bioinformatics pipeline for processing FASTA and GFF3 files, extracting genomic features, calculating sequence statistics, and generating publication-ready visualizations.

## Project Overview

This project demonstrates an end-to-end workflow for genomic data analysis using Python and Biopython.

The initial dataset is the **Escherichia coli K-12 MG1655 reference genome (ASM584v2)**. The project uses a matching genomic FASTA file and GFF3 annotation file from NCBI.

### Workflow

```text
Genome FASTA + GFF3 Annotation
            |
            v
      FASTA Parsing
            |
            v
 Sequence Statistics
            |
            +----> GC Content
            +----> AT Content
            +----> Sequence Length
            +----> N Content
            |
            v
       GFF3 Parsing
            |
            v
     Feature Extraction
            |
            +----> Genes
            +----> CDS
            +----> Exons
            +----> Other Features
            |
            v
    Statistical Analysis
            |
            v
      Visualization
            |
            v
 CSV Results + Figures
```

## Objectives

- Parse genomic FASTA files using Biopython.
- Parse GFF3 genome annotations.
- Extract genes, CDS, exons, and other genomic features.
- Calculate sequence-level statistics.
- Calculate genomic feature lengths and distributions.
- Generate clear, publication-quality figures.
- Export processed results as CSV files.
- Create a reproducible Python workflow suitable for research and freelance bioinformatics work.

## Dataset

### Organism

**Escherichia coli K-12 MG1655**

### Assembly

**GCF_000005845.2_ASM584v2**

### Source

NCBI Genome Assembly / RefSeq

Official assembly directory:

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/

### Input files

Download the matching files from NCBI:

- `GCF_000005845.2_ASM584v2_genomic.fna.gz`
- `GCF_000005845.2_ASM584v2_genomic.gff.gz`

After extraction, rename them for convenience:

```text
data/
├── genome.fasta
└── annotation.gff3
```

> The original NCBI filenames and accession should be retained in the project documentation so the analysis remains reproducible.

## Project Structure

```text
genomic-feature-analysis-pipeline/
│
├── data/
│   ├── genome.fasta
│   └── annotation.gff3
│
├── src/
│   ├── parse_fasta.py
│   ├── parse_gff.py
│   ├── statistics.py
│   ├── visualize.py
│   └── extract_features.py
│
├── results/
│   ├── sequence_statistics.csv
│   ├── genes.csv
│   ├── cds.csv
│   └── exons.csv
│
├── figures/
│   ├── gene_length_distribution.png
│   ├── gc_content_distribution.png
│   ├── chromosome_gene_distribution.png
│   └── feature_distribution.png
│
├── notebooks/
│   └── exploratory_analysis.ipynb
│
├── main.py
├── requirements.txt
├── .gitignore
└── README.md
```

## Technologies

- Python
- Biopython
- pandas
- NumPy
- matplotlib
- Jupyter Notebook
- Git/GitHub

## Installation

Create a virtual environment:

```bash
python -m venv bioinformatics-env
```

Activate it on Windows:

```bash
bioinformatics-env\Scripts\activate
```

Activate it on Linux/macOS:

```bash
source bioinformatics-env/bin/activate
```

Install dependencies:

```bash
pip install biopython pandas numpy matplotlib jupyter
```

Save dependencies:

```bash
pip freeze > requirements.txt
```

## Analysis Modules

### 1. FASTA Analysis

The FASTA module will calculate:

- Sequence ID
- Sequence length
- GC content
- AT content
- Number of ambiguous `N` bases

Example output:

```text
sequence_id,length,gc_percent,at_percent,n_count
...
```

### 2. GFF3 Analysis

The GFF3 parser will process annotation fields including:

- `seqid`
- `source`
- `type`
- `start`
- `end`
- `score`
- `strand`
- `phase`
- `attributes`

A derived feature-length column will be calculated as:

```text
length = end - start + 1
```

### 3. Feature Extraction

The pipeline will identify and export relevant annotation types, including:

- Genes
- CDS
- Exons
- Other annotated genomic features

The extracted tables will be saved under `results/`.

### 4. Statistical Analysis

The project will calculate statistics such as:

- Number of genomic features
- Mean feature length
- Minimum feature length
- Maximum feature length
- Feature counts by type
- Sequence length statistics
- GC-content statistics

### 5. Visualization

The pipeline will generate figures such as:

- Gene length distribution
- GC-content distribution
- Feature type distribution
- Gene distribution by sequence/replicon

Figures should be saved at high resolution for use in reports and presentations.

## Example Commands

Once the pipeline is integrated, the preferred interface will be:

```bash
python main.py --fasta data/genome.fasta --gff data/annotation.gff3
```

Expected workflow:

```text
========================================
 GENOMIC FEATURE ANALYSIS
========================================

Reading FASTA...
Reading GFF3...

Calculating sequence statistics...
Extracting genomic features...
Generating visualizations...

Results saved to: results/
Figures saved to: figures/

Analysis completed successfully.
========================================
```

## Expected Outputs

### Sequence Statistics

```text
results/sequence_statistics.csv
```

Expected columns:

```text
sequence_id
length
gc_percent
at_percent
n_count
```

### Gene Annotation

```text
results/genes.csv
```

### CDS Annotation

```text
results/cds.csv
```

### Exon Annotation

```text
results/exons.csv
```

### Figures

```text
figures/
├── gene_length_distribution.png
├── gc_content_distribution.png
├── chromosome_gene_distribution.png
└── feature_distribution.png
```

## Reproducibility

The project should make it possible for another user to reproduce the analysis by:

1. Downloading the same NCBI dataset.
2. Placing the FASTA and GFF3 files in `data/`.
3. Installing the dependencies.
4. Running `main.py`.
5. Reviewing the generated CSV files and figures.

The exact dataset accession and source should always be documented.

## Quality Checks

Before considering the pipeline complete, verify:

- FASTA files can be parsed without errors.
- GFF3 comments and metadata lines are handled correctly.
- Coordinates are interpreted correctly.
- Feature lengths are calculated correctly.
- Empty or missing attributes do not crash the pipeline.
- Output CSV files contain the expected columns.
- Figures are generated successfully.
- Results are reproducible from a clean environment.

## Future Improvements

Potential extensions include:

- FASTQ quality-control analysis.
- Sequence motif detection.
- ORF prediction.
- Protein translation.
- GC-content by genomic window.
- Codon usage analysis.
- Interactive dashboards.
- Automated HTML/PDF reports.
- Command-line configuration using `argparse`.
- Unit tests using `pytest`.
- Docker containerization.
- Workflow automation using Snakemake or Nextflow.

## Freelancing Relevance

This project is designed to demonstrate services that can be offered to researchers, students, laboratories, and biotechnology teams.

Potential services include:

- FASTA sequence analysis
- GFF3 annotation processing
- Genomic feature extraction
- Biological data cleaning
- Python bioinformatics scripting
- Genomic data visualization
- Reproducible analysis pipelines

The project should demonstrate practical ability rather than simply listing bioinformatics technologies on a CV.

## Important Scope Note

This first project focuses primarily on **genomic sequence and annotation processing**.

The next portfolio project should use a small eukaryotic dataset to demonstrate transcript-level analysis, including:

- Gene-expression analysis
- Transcript processing
- 3' UTR analysis
- PCA
- Differential expression concepts
- Heatmaps
- Volcano plots

## License

This project's source code can be released under the MIT License unless a different license is selected.

Dataset files remain subject to the terms and policies of their original data provider.
