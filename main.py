#!/usr/bin/env python3
"""Genomic Feature Analysis Pipeline - Main entry point."""

import argparse
import sys
import os


from src.parse_fasta import parse_fasta, calculate_sequence_stats
from src.parse_gff import parse_gff
from src.extract_features import extract_all_features, save_features
from src.statistics import (sequence_stats_summary, feature_stats_summary,
                            gene_stats, cds_stats, exon_stats, print_summary)
from src.visualize import generate_all_figures


def run_pipeline(fasta_path, gff_path, output_dir='results', figures_dir='figures'):
    """Run the complete genomic feature analysis pipeline."""
    print("=" * 50)
    print(" GENOMIC FEATURE ANALYSIS")
    print("=" * 50)

    # Parse FASTA
    print("\nReading FASTA...")
    records = parse_fasta(fasta_path)
    print(f"  Found {len(records)} sequence(s)")

    print("Calculating sequence statistics...")
    seq_stats = calculate_sequence_stats(records)

    # Parse GFF3
    print("\nReading GFF3...")
    gff_df = parse_gff(gff_path)
    print(f"  Found {len(gff_df)} features")

    # Extract features
    print("\nExtracting genomic features...")
    features = extract_all_features(gff_df)
    print(f"  Genes: {len(features['genes'])}")
    print(f"  CDS: {len(features['cds'])}")
    print(f"  Exons: {len(features['exons'])}")
    print(f"  3' UTRs: {len(features['utr3'])}")

    # Save results
    print("\nSaving results...")
    save_features(features, output_dir)
    seq_stats.to_csv(os.path.join(output_dir, 'sequence_statistics.csv'), index=False)
    print(f"  Sequence statistics saved to {output_dir}/sequence_statistics.csv")

    # Generate statistics
    print("\nCalculating summary statistics...")
    seq_summary = sequence_stats_summary(seq_stats)
    feat_summary = feature_stats_summary(gff_df)
    gene_summary = gene_stats(features['genes'])
    cds_summary = cds_stats(features['cds'])
    exon_summary = exon_stats(features['exons'])

    print_summary(seq_summary, "Sequence Statistics")
    print_summary(feat_summary, "Feature Statistics")
    print_summary(gene_summary, "Gene Statistics")
    print_summary(cds_summary, "CDS Statistics")
    print_summary(exon_summary, "Exon Statistics")

    # Generate visualizations
    print("\nGenerating visualizations...")
    generate_all_figures(
        seq_stats, gff_df,
        features['genes'], features['cds'], features['exons'],
        figures_dir
    )

    print("\n" + "=" * 50)
    print(" Analysis completed successfully.")
    print("=" * 50)
    print(f"\nResults saved to: {output_dir}/")
    print(f"Figures saved to: {figures_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description='Genomic Feature Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
  python main.py --fasta data/genome.fasta --gff data/annotation.gff3
        '''
    )
    parser.add_argument('--fasta', required=True, help='Path to genome FASTA file')
    parser.add_argument('--gff', required=True, help='Path to GFF3 annotation file')
    parser.add_argument('--out', default='results', help='Output directory for CSV results (default: results)')
    parser.add_argument('--figures', default='figures', help='Output directory for figures (default: figures)')

    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file not found: {args.fasta}")
        sys.exit(1)
    if not os.path.exists(args.gff):
        print(f"Error: GFF3 file not found: {args.gff}")
        sys.exit(1)

    run_pipeline(args.fasta, args.gff, args.out, args.figures)


if __name__ == '__main__':
    main()