import unittest
import os
import pandas as pd
from src.parse_fasta import parse_fasta, calculate_sequence_stats
from src.parse_gff import parse_gff

class TestGenomicPipeline(unittest.TestCase):
    def setUp(self):
        self.test_fasta = "tests/test.fasta"
        self.test_gff = "tests/test.gff3"

        with open(self.test_fasta, "w") as f:
            f.write(">seq1\nATGCATGC\n>seq2\nNNNNNNNN\n")

        with open(self.test_gff, "w") as f:
            f.write("##gff-version 3\n")
            f.write("seq1\tsource\tgene\t1\t10\t.\t+\t.\tID=gene1;Name=test_gene\n")
            f.write("seq1\tsource\tCDS\t2\t8\t.\t+\t.\tID=cds1;Parent=gene1\n")

    def tearDown(self):
        if os.path.exists(self.test_fasta):
            os.remove(self.test_fasta)
        if os.path.exists(self.test_gff):
            os.remove(self.test_gff)

    def test_fasta_parsing(self):
        records = parse_fasta(self.test_fasta)
        self.assertEqual(len(records), 2)
        self.assertEqual(str(records[0].seq), "ATGCATGC")

    def test_sequence_stats(self):
        records = parse_fasta(self.test_fasta)
        stats_df = calculate_sequence_stats(records)

        # seq1: ATGCATGC -> 50% GC, 50% AT
        seq1_stats = stats_df[stats_df['sequence_id'] == 'seq1'].iloc[0]
        self.assertEqual(seq1_stats['gc_percent'], 50.0)
        self.assertEqual(seq1_stats['at_percent'], 50.0)

        # seq2: NNNNNNNN -> 0% GC, 0% AT (due to total_bases = 0)
        seq2_stats = stats_df[stats_df['sequence_id'] == 'seq2'].iloc[0]
        self.assertEqual(seq2_stats['gc_percent'], 0.0)
        self.assertEqual(seq2_stats['at_percent'], 0.0)

    def test_gff_parsing(self):
        gff_df = parse_gff(self.test_gff)
        self.assertEqual(len(gff_df), 2)
        self.assertEqual(gff_df.iloc[0]['type'], 'gene')
        self.assertEqual(gff_df.iloc[1]['type'], 'CDS')

if __name__ == '__main__':
    unittest.main()
