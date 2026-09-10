"""Extract specific genomic features from parsed GFF3."""

import pandas as pd


def extract_genes(gff_df):
    """Extract gene features."""
    return filter_by_type(gff_df, 'gene')


def extract_cds(gff_df):
    """Extract CDS features."""
    return filter_by_type(gff_df, 'CDS')


def extract_exons(gff_df):
    """Extract exon features."""
    return filter_by_type(gff_df, 'exon')


def extract_utr3(gff_df):
    """Extract 3' UTR features."""
    return filter_by_type(gff_df, 'three_prime_UTR')


def extract_all_features(gff_df):
    """Extract all feature types into separate DataFrames."""
    return {
        'genes': extract_genes(gff_df),
        'cds': extract_cds(gff_df),
        'exons': extract_exons(gff_df),
        'utr3': extract_utr3(gff_df),
        'all': gff_df
    }


def save_features(features_dict, output_dir):
    """Save feature DataFrames to CSV files."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Map internal keys to requested filenames
    filename_map = {
        'genes': 'genes.csv',
        'cds': 'cds.csv',
        'exons': 'exons.csv',
        'utr3': 'utr3.csv'
    }

    for key, filename in filename_map.items():
        if key in features_dict:
            df = features_dict[key]
            path = os.path.join(output_dir, filename)
            df.to_csv(path, index=False)
            print(f"Saved {len(df)} {key} to {path}")


def filter_by_type(df, feature_type):
    """Filter features by type."""
    return df[df['type'] == feature_type].copy()


if __name__ == "__main__":
    import sys
    sys.path.insert(0, 'src')
    from parse_gff import parse_gff

    gff_path = sys.argv[1] if len(sys.argv) > 1 else "data/annotation.gff3"
    df = parse_gff(gff_path)
    features = extract_all_features(df)
    save_features(features, "results")