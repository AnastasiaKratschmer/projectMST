import unittest

from sp_exact_3 import sp_exact_3, backtrack_3seq
from sp_approx import sp_approx
from helpers.utils import create_score_matrix, parse_fasta_multiple, parse_fasta, get_cost_3

#THIS IS THE CODE FOR THE METHOD SECTION OF THE REPORT

s1_short = "GTTCCGAAAGGCTAGCGCTAGGCGCC".lower()
s2_short = "ATGGATTTATCTGCTCTTCG".lower()
s3_short = "TGCATGCTGAAACTTCTCAACCA".lower()
test_short_res = 198

s1_long = "GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGATAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTCTCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTG".lower()
s2_long = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAACGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGA".lower()
s3_long = "CGCTGGTGCAACTCGAAGACCTATCTCCTTCCCGGGGGGGCTTCTCCGGCATTTAGGCCTCGGCGTTTGGAAGTACGGAGGTTTTTCTCGGAAGAAAGTTCACTGGAAGTGGAAGAAATGGATTTATCTGCTGTTCGAATTCAAGAAGTACAAAATGTCCTTCATGCTATGCAGAAAATCTTGGAGTGTCCAATCTGTTT".lower()
test_long_res = 1482

score_matrix_file = "testdata/score_matrix.txt"
score_matrix = create_score_matrix(score_matrix_file)
gap = 5

def get_cost_2(cost_matrix):
    n, m = tuple(map(lambda x: x - 1, cost_matrix.shape))
    return cost_matrix[n, m]


class TestUtils(unittest.TestCase):
    def test_score_matrix_creation(self):
        # from examples
        score_matrix = {
            'a': {'a': 0, 'c': 5, 'g': 2, 't': 5},
            'c': {'a': 5, 'c': 0, 'g': 5, 't': 2},
            'g': {'a': 2, 'c': 5, 'g': 0, 't': 5},
            't': {'a': 5, 'c': 2, 'g': 5, 't': 0}}
        created_score_matrix = create_score_matrix(score_matrix_file)
        self.assertEqual(score_matrix, created_score_matrix)


class TestAlgs(unittest.TestCase):
    def test_cost_short(self):
        res_exact_short = get_cost_3(sp_exact_3(s1_short, s2_short, s3_short, score_matrix, gap_cost = gap))
        self.assertEqual(test_short_res, res_exact_short)

        res, _ = sp_approx([s1_short, s2_short, s3_short], score_matrix, gap)
        self.assertLessEqual(res_exact_short, res)
    
    def test_cost_long(self):
        res_exact_long = get_cost_3(sp_exact_3(s1_long, s2_long, s3_long, score_matrix, gap))
        self.assertEqual(test_long_res, res_exact_long)

        res, _ = sp_approx([s1_long, s2_long, s3_long], score_matrix, gap)
        self.assertLessEqual(res_exact_long, res)

        
if __name__ == '__main__':
    unittest.main()
