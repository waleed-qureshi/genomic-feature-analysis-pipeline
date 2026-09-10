"""GFF3 parsing and feature extraction."""

import pandas as pd


def parse_gff(gff_path):
    """Parse GFF3 file and return DataFrame with features."""
    features = []
    with open(gff_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attributes = parts

            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val

            try:
                start_val = int(start)
                end_val = int(end)
            except ValueError:
                continue

            feature = {
                'seqid': seqid,
                'source': source,
                'type': ftype,
                'start': start_val,
                'end': end_val,
                'score': score if score != '.' else None,
                'strand': strand,
                'phase': phase if phase != '.' else None,
                'length': end_val - start_val + 1,
            }

            # Update attributes without overwriting core fields
            for key, val in attr_dict.items():
                if key not in feature:
                    feature[key] = val
            features.append(feature)

    return pd.DataFrame(features)


def get_feature_types(df):
    """Return unique feature types in the annotation."""
    return sorted(df['type'].unique())


def filter_by_type(df, feature_type):
    """Filter features by type."""
    return df[df['type'] == feature_type].copy()


if __name__ == "__main__":
    import sys
    gff_path = sys.argv[1] if len(sys.argv) > 1 else "data/annotation.gff3"
    df = parse_gff(gff_path)
    print(f"Total features: {len(df)}")
    print(f"Feature types: {get_feature_types(df)}")
    print(df.head().to_string())