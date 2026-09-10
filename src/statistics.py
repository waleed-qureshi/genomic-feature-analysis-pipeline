"""Statistical analysis for sequences and features."""

import pandas as pd
import numpy as np


def sequence_stats_summary(seq_stats_df):
    """Calculate summary statistics for sequence data."""
    if seq_stats_df.empty:
        return {}

    return {
        'num_sequences': len(seq_stats_df),
        'total_length': int(seq_stats_df['length'].sum()),
        'mean_length': round(seq_stats_df['length'].mean(), 2),
        'median_length': round(seq_stats_df['length'].median(), 2),
        'min_length': int(seq_stats_df['length'].min()),
        'max_length': int(seq_stats_df['length'].max()),
        'mean_gc_percent': round(seq_stats_df['gc_percent'].mean(), 2),
        'median_gc_percent': round(seq_stats_df['gc_percent'].median(), 2),
        'min_gc_percent': round(seq_stats_df['gc_percent'].min(), 2),
        'max_gc_percent': round(seq_stats_df['gc_percent'].max(), 2),
        'total_n_count': int(seq_stats_df['n_count'].sum()),
    }


def feature_stats_summary(features_df):
    """Calculate summary statistics for genomic features."""
    if features_df.empty:
        return {}

    stats = {
        'total_features': len(features_df),
        'feature_types': features_df['type'].value_counts().to_dict(),
        'mean_length': round(features_df['length'].mean(), 2),
        'median_length': round(features_df['length'].median(), 2),
        'min_length': int(features_df['length'].min()),
        'max_length': int(features_df['length'].max()),
        'strand_counts': features_df['strand'].value_counts().to_dict(),
    }

    # Per-type length stats
    type_stats = {}
    for ftype in features_df['type'].unique():
        subset = features_df[features_df['type'] == ftype]
        type_stats[ftype] = {
            'count': len(subset),
            'mean_length': round(subset['length'].mean(), 2),
            'median_length': round(subset['length'].median(), 2),
            'min_length': int(subset['length'].min()),
            'max_length': int(subset['length'].max()),
        }
    stats['by_type'] = type_stats

    return stats


def gene_stats(genes_df):
    """Calculate gene-specific statistics."""
    if genes_df.empty:
        return {}

    # Extract gene names if available
    gene_names = []
    if 'Name' in genes_df.columns:
        gene_names = genes_df['Name'].dropna().tolist()
    elif 'gene' in genes_df.columns:
        gene_names = genes_df['gene'].dropna().tolist()

    return {
        'num_genes': len(genes_df),
        'mean_length': round(genes_df['length'].mean(), 2),
        'median_length': round(genes_df['length'].median(), 2),
        'min_length': int(genes_df['length'].min()),
        'max_length': int(genes_df['length'].max()),
        'strand_distribution': genes_df['strand'].value_counts().to_dict(),
        'sample_genes': gene_names[:10] if gene_names else [],
    }


def cds_stats(cds_df):
    """Calculate CDS-specific statistics."""
    if cds_df.empty:
        return {}

    return {
        'num_cds': len(cds_df),
        'mean_length': round(cds_df['length'].mean(), 2),
        'median_length': round(cds_df['length'].median(), 2),
        'min_length': int(cds_df['length'].min()),
        'max_length': int(cds_df['length'].max()),
        'phase_distribution': cds_df['phase'].value_counts().to_dict() if 'phase' in cds_df.columns else {},
    }


def exon_stats(exons_df):
    """Calculate exon-specific statistics."""
    if exons_df.empty:
        return {}

    return {
        'num_exons': len(exons_df),
        'mean_length': round(exons_df['length'].mean(), 2),
        'median_length': round(exons_df['length'].median(), 2),
        'min_length': int(exons_df['length'].min()),
        'max_length': int(exons_df['length'].max()),
    }


def print_summary(stats_dict, title="Summary"):
    """Pretty print statistics dictionary."""
    print(f"\n{'='*50}")
    print(f" {title}")
    print(f"{'='*50}")
    for key, val in stats_dict.items():
        if isinstance(val, dict):
            print(f"\n{key}:")
            for k, v in val.items():
                print(f"  {k}: {v}")
        else:
            print(f"{key}: {val}")


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

    print_summary(sequence_stats_summary(seq_stats), "Sequence Statistics")
    print_summary(feature_stats_summary(gff_df), "Feature Statistics")
    print_summary(gene_stats(extract_genes(gff_df)), "Gene Statistics")
    print_summary(cds_stats(extract_cds(gff_df)), "CDS Statistics")
    print_summary(exon_stats(extract_exons(gff_df)), "Exon Statistics")