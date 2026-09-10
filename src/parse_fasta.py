"""FASTA parsing and sequence statistics calculation."""

from Bio import SeqIO
import pandas as pd


def parse_fasta(fasta_path):
    """Parse FASTA file and return list of SeqRecord objects."""
    return list(SeqIO.parse(fasta_path, "fasta"))


def calculate_sequence_stats(records):
    """Calculate statistics for each sequence in records."""
    stats = []
    for record in records:
        seq = str(record.seq).upper()
        length = len(seq)
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        n_count = seq.count('N')
        total_bases = gc_count + at_count

        gc_percent = (gc_count / total_bases * 100) if total_bases > 0 else 0
        at_percent = (at_count / total_bases * 100) if total_bases > 0 else 0

        stats.append({
            'sequence_id': record.id,
            'length': length,
            'gc_percent': round(gc_percent, 2),
            'at_percent': round(at_percent, 2),
            'n_count': n_count
        })
    return pd.DataFrame(stats)


def get_sequence_dict(records):
    """Return dict mapping sequence_id to sequence string."""
    return {record.id: str(record.seq).upper() for record in records}


if __name__ == "__main__":
    import sys
    fasta_path = sys.argv[1] if len(sys.argv) > 1 else "data/genome.fasta"
    records = parse_fasta(fasta_path)
    df = calculate_sequence_stats(records)
    print(df.to_string(index=False))