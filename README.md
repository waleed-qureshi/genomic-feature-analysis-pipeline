# 🧬 Genomic Feature Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![Biopython](https://img.shields.io/badge/Library-Biopython-green)](https://biopython.org/)
[![Pandas](https://img.shields.io/badge/Library-Pandas-orange)](https://pandas.pydata.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end bioinformatics pipeline designed to automate the extraction, statistical analysis, and visualization of genomic features from reference genomes. This project transforms raw NCBI genomic data (FASTA and GFF3) into publication-ready datasets and figures, demonstrating a professional approach to bioinformatics software engineering.

## 🚀 Project Impact & Value

In bioinformatics, the transition from raw annotation files to a structured analysis is often a manual and error-prone process. This pipeline solves that by providing a **reproducible, scalable workflow** that ensures data integrity from parsing to visualization. 

**Key outcomes include:**
- **Automated Parsing**: Efficiently handles large-scale GFF3 and FASTA files.
- **Feature Engineering**: Extracts high-value biological features (Genes, CDS, Exons, 3' UTRs) into structured CSVs for downstream analysis.
- **Statistical Rigor**: Calculates sequence-level metrics (GC/AT content, ambiguity counts) and feature-length distributions.
- **Professional Visualization**: Generates high-resolution, publication-quality figures (DPI 300) for immediate use in research papers or reports.

---

## 🛠️ Technical Competencies Demonstrated

This project serves as a portfolio of the following professional skills:

### 🧬 Bioinformatics & Domain Expertise
- **Genome Data Processing**: Expert handling of FASTA and GFF3 standards.
- **Biological Analysis**: Implementation of GC content analysis and genomic feature length distributions.
- **Standard Compliance**: Adherence to NCBI reference genome standards.

### 💻 Software Engineering
- **Modular Architecture**: Decoupled logic across parsing, extraction, statistics, and visualization modules for maximum maintainability.
- **Robust Error Handling**: Implementation of guards against malformed GFF3 entries and empty datasets.
- **CLI Development**: Built a user-friendly command-line interface using `argparse` for flexible input/output management.

### 📊 Data Engineering & Visualization
- **Data Wrangling**: Leveraging `pandas` and `numpy` for efficient manipulation of genomic tables.
- **Data Visualization**: Advanced use of `matplotlib` to create professional histograms, pie charts, and boxplots.
- **Exploratory Analysis**: Development of an interactive Jupyter Notebook for deep-dive data inspection.

---

## 📂 Project Structure

```text
genomic-feature-analysis-pipeline/
├── data/                 # Input NCBI reference files (FASTA, GFF3)
├── src/                  # Core modular logic
│   ├── parse_fasta.py    # FASTA parsing & sequence stats
│   ├── parse_gff.py      # GFF3 parsing & attribute extraction
│   ├── extract_features.py# Feature-specific filtering logic
│   ├── statistics.py     # Summary statistics calculations
│   └── visualize.py      # Publication-quality plotting
├── results/              # Processed CSV datasets (Genes, CDS, UTRs, etc.)
├── figures/              # High-resolution PNG visualizations
├── notebooks/            # Exploratory data analysis (.ipynb)
├── main.py               # Pipeline orchestration entry point
└── requirements.txt      # Dependency management
```

---

## 🚦 Getting Started

### Installation
```bash
# Clone the repository
git clone https://github.com/yourusername/genomic-feature-analysis-pipeline.git
cd genomic-feature-analysis-pipeline

# Set up environment
python -m venv bioenv
source bioenv/bin/activate  # Linux/macOS
# bioenv\Scripts\activate    # Windows

# Install dependencies
pip install -r requirements.txt
```

### Usage
```bash
python main.py --fasta data/genome.fasta --gff data/annotation.gff3
```

---

## 📊 Data Provenance & Actual Results

To ensure scientific reproducibility and credibility, this pipeline was validated using a gold-standard public dataset from the NCBI RefSeq database.

### 🧬 Dataset Specifications
- **Organism**: *Escherichia coli* str. K-12 substr. MG1655
- **Assembly**: [ASM584v2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2/) (GCF_000005845.2)
- **Accession**: `NC_000913.3`
- **Source**: NCBI National Center for Biotechnology Information

### 📈 Pipeline Outputs (Validated Results)
The following metrics were extracted and verified from the reference genome:

| Metric | Value | File Reference |
| :--- | :--- | :--- |
| **Total Genome Length** | 4,641,652 bp | `results/sequence_statistics.csv` |
| **GC Content** | 50.79% | `results/sequence_statistics.csv` |
| **AT Content** | 49.21% | `results/sequence_statistics.csv` |
| **Ambiguity Count (N)** | 0 | `results/sequence_statistics.csv` |

**Detailed Feature Extraction:**
The pipeline successfully partitioned the genome into structured datasets, including:
- **Genes**: Comprehensive list of genomic loci (e.g., `thrL`, `thrA`, `thrB`) with precise start/end coordinates and locus tags.
- **CDS**: Protein-coding sequences mapped to their respective parent genes.
- **Exons/UTRs**: Structural boundaries for translation initiation and termination.

These results are exported as clean CSVs in the `results/` directory, ready for immediate import into downstream analysis tools like R or Python.

---

## 🎓 Professional Application
This pipeline is a blueprint for services I offer to research labs and biotech teams, including:
- Custom bioinformatics script development.
- Raw genomic data cleaning and restructuring.
- Automated figure generation for publications.
- Reproducible analysis pipeline design.

**License:** MIT  
**Dataset:** Escherichia coli K-12 MG1655 (NCBI ASM584v2)
