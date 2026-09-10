"""Visualization module for genomic feature analysis."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


def setup_style():
    """Set up consistent plot style."""
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.figsize': (10, 6),
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'font.size': 11,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
    })


def plot_gene_length_distribution(genes_df, output_path):
    """Plot histogram of gene lengths."""
    if genes_df.empty:
        print(f"Skipping {output_path}: No genes found.")
        return
    setup_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    lengths = genes_df['length'].dropna()
    ax.hist(lengths, bins=50, edgecolor='black', alpha=0.7, color='#2E86AB')
    ax.axvline(lengths.mean(), color='#E74C3C', linestyle='--', linewidth=2, label=f'Mean: {lengths.mean():.0f} bp')
    ax.axvline(lengths.median(), color='#F39C12', linestyle='--', linewidth=2, label=f'Median: {lengths.median():.0f} bp')

    ax.set_xlabel('Gene Length (bp)')
    ax.set_ylabel('Frequency')
    ax.set_title('Gene Length Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_gc_content_distribution(seq_stats_df, output_path):
    """Plot GC content distribution across sequences."""
    if seq_stats_df.empty:
        print(f"Skipping {output_path}: No sequence stats found.")
        return
    setup_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    gc_vals = seq_stats_df['gc_percent'].dropna()
    ax.hist(gc_vals, bins=30, edgecolor='black', alpha=0.7, color='#27AE60')
    ax.axvline(gc_vals.mean(), color='#E74C3C', linestyle='--', linewidth=2, label=f'Mean: {gc_vals.mean():.2f}%')
    ax.axvline(gc_vals.median(), color='#F39C12', linestyle='--', linewidth=2, label=f'Median: {gc_vals.median():.2f}%')

    ax.set_xlabel('GC Content (%)')
    ax.set_ylabel('Number of Sequences')
    ax.set_title('GC Content Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_chromosome_gene_distribution(genes_df, output_path):
    """Plot gene count per chromosome/sequence."""
    if genes_df.empty:
        print(f"Skipping {output_path}: No genes found.")
        return
    setup_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    if 'seqid' in genes_df.columns:
        chr_counts = genes_df['seqid'].value_counts().sort_index()
    else:
        chr_counts = pd.Series([len(genes_df)], index=['All'])

    bars = ax.bar(range(len(chr_counts)), chr_counts.values, color='#8E44AD', edgecolor='black', alpha=0.7)
    ax.set_xticks(range(len(chr_counts)))
    ax.set_xticklabels(chr_counts.index, rotation=45, ha='right')
    ax.set_xlabel('Chromosome / Sequence')
    ax.set_ylabel('Number of Genes')
    ax.set_title('Gene Distribution by Chromosome')
    ax.grid(True, alpha=0.3, axis='y')

    # Add value labels on bars
    for bar, val in zip(bars, chr_counts.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(chr_counts.values)*0.01,
                str(val), ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_feature_distribution(gff_df, output_path):
    """Plot feature type distribution."""
    if gff_df.empty:
        print(f"Skipping {output_path}: No features found.")
        return
    setup_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    type_counts = gff_df['type'].value_counts()
    # Limit to top 15 for readability
    if len(type_counts) > 15:
        top15 = type_counts.head(15)
        other_count = type_counts[15:].sum()
        type_counts = pd.concat([top15, pd.Series({'Other': other_count})])

    colors = plt.cm.Set3(np.linspace(0, 1, len(type_counts)))
    wedges, texts, autotexts = ax.pie(type_counts.values, labels=type_counts.index,
                                        autopct='%1.1f%%', colors=colors, startangle=90)

    ax.set_title('Genomic Feature Type Distribution')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_gene_length_boxplot(genes_df, output_path):
    """Plot boxplot of gene lengths by strand."""
    if genes_df.empty:
        print(f"Skipping {output_path}: No genes found.")
        return
    setup_style()
    fig, ax = plt.subplots(figsize=(8, 6))

    if 'strand' in genes_df.columns:
        data = [genes_df[genes_df['strand'] == s]['length'].dropna() for s in ['+', '-'] if s in genes_df['strand'].values]
        labels = [s for s in ['+', '-'] if s in genes_df['strand'].values]
        ax.boxplot(data, tick_labels=labels)
        ax.set_xlabel('Strand')
        ax.set_ylabel('Gene Length (bp)')
        ax.set_title('Gene Length Distribution by Strand')
    else:
        ax.boxplot([genes_df['length'].dropna()], tick_labels=['All'])
        ax.set_ylabel('Gene Length (bp)')
        ax.set_title('Gene Length Distribution')

    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def generate_all_figures(seq_stats_df, gff_df, genes_df, cds_df, exons_df, output_dir):
    """Generate all publication figures."""
    os.makedirs(output_dir, exist_ok=True)

    plot_gene_length_distribution(genes_df, os.path.join(output_dir, 'gene_length_distribution.png'))
    plot_gc_content_distribution(seq_stats_df, os.path.join(output_dir, 'gc_content_distribution.png'))
    plot_chromosome_gene_distribution(genes_df, os.path.join(output_dir, 'chromosome_gene_distribution.png'))
    plot_feature_distribution(gff_df, os.path.join(output_dir, 'feature_distribution.png'))
    plot_gene_length_boxplot(genes_df, os.path.join(output_dir, 'gene_length_by_strand.png'))


if __name__ == "__main__":
    import sys
    sys.path.insert(0, 'src')
    from parse_fasta import parse_fasta, calculate_sequence_stats
    from parse_gff import parse_gff
    from extract_features import extract_genes, extract_cds, extract_exons

    fasta_path = sys.argv[1] if len(sys.argv) > 1 else "data/genome.fasta"
    gff_path = sys.argv[2] if len(sys.argv) > 2 else "data/annotation.gff3"

    records = parse_fasta(fasta_path)
    seq_stats = calculate_sequence_stats(records)
    gff_df = parse_gff(gff_path)
    genes_df = extract_genes(gff_df)
    cds_df = extract_cds(gff_df)
    exons_df = extract_exons(gff_df)

    generate_all_figures(seq_stats, gff_df, genes_df, cds_df, exons_df, "figures")