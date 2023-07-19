import numpy as np
import sys
import os
from tqdm import tqdm # loading bar
from helpers.utils import linear_C, get_cost_2, get_cost_3,  get_sequence_string, parse_fasta, create_score_matrix, write_alignments_to_file


'''''
Pseudocode for recurrences fpr a nonboundary cell (i,j) 
(as written in Gusfield p. 344, section 14.6.1)

for i := 1 to n1 do
    for j := 1 to n2 do
        for k := 1 to n3 do
            begin
            if (S1(i) == S2(j)) then cij := smatch
            else cij := smis;
            if (S1(i) = S3(k)) then cik := smatch
            else cik := smis;
            if (S2(j) = S3(k)) then cjk := smatch
            else cjk := smis;

            d1 := D(i-1, j-1, k-1) + cij + cik + cjk;
            d2 := D(i-1, j-1, k) + cij + 2*sspace;
            d3 := D(i-1, j, k-1) + cik + 2*sspace;
            d4 := D(i, j-1, k-1) + cjk + 2*sspace;
            d5 := D(i-1, j, k) + 2*sspace;
            d6 := D(i, j-1, k) + 2*sspace;
            d7 := D(i, j, k-1) + 2*sspace;

            D(i, j, k) := Min ( d1, d2, d3, d4, d5, d6, d7);
            end;

'''

def sp_exact_3(seq1: str, seq2: str, seq3: str, score_matrix: dict, gap_cost: int):
    n1 = len(seq1)
    n2 = len(seq2)
    n3 = len(seq3)

    D = np.full((n1+1, n2+1, n3+1), np.infty)

    # Handle boundary cells (Gusfield section 14.6.1, p. 345)
    D[0, 0, 0] = 0
    for i in tqdm(range(1, n1+1)): 
        for j in range(1, n2+1):
            for k in range(1, n3+1):
                 D[i, j, 0] = get_cost_2(linear_C(gap_cost, score_matrix, seq1[0:i], seq2[0:j])) + (i+j)*gap_cost
                 D[i, 0, k] = get_cost_2(linear_C(gap_cost, score_matrix, seq1[0:i], seq3[0:k])) + (i+k)*gap_cost
                 D[0, j, k] = get_cost_2(linear_C(gap_cost, score_matrix, seq2[0:j], seq3[0:k])) + (j+k)*gap_cost
    
    
    # Handle the non-boundary cells (Gusfield section 14.6.1, p. 344)
    for i in range(1, n1+1):
        for j in range(1, n2+1):
            for k in range(1, n3+1):
                cij = score_matrix[seq1[i-1]][seq2[j-1]]
                cik = score_matrix[seq1[i-1]][seq3[k-1]]
                cjk = score_matrix[seq2[j-1]][seq3[k-1]]
                
                d1 = D[i-1, j-1, k-1] + cij + cik + cjk
                d2 = D[i-1, j-1, k] + cij + 2 * gap_cost
                d3 = D[i-1, j, k-1] + cik + 2 * gap_cost
                d4 = D[i, j-1, k-1] + cjk + 2 * gap_cost
                d5 = D[i-1, j, k] + 2 * gap_cost
                d6 = D[i, j-1, k] + 2 * gap_cost
                d7 = D[i, j, k-1] + 2 * gap_cost
                
                D[i, j, k] = min(d1, d2, d3, d4, d5, d6, d7)
    return D


# ------------------------------ ACTUAL PROGRAM ------------------------------
if __name__ == "__main__":
    # ------ Parameter setup
    # Read arguments ...
    if not len(sys.argv) == 6:
            sys.stderr.write("USAGE: python3 %s < three FASTA files or sequences > "
                            "< gap cost >  < score matrix >\n" % sys.argv[0])
            sys.exit(1)
    seq1_input, seq2_input, seq3_input, gap_input, score_matrix_input = sys.argv[1:]

    # Set sequences (extract from file if necessary) ...
    seq1 = get_sequence_string(seq1_input)
    seq2 = get_sequence_string(seq2_input)
    seq3 = get_sequence_string(seq3_input)

    
    # Set gap and score matrix ...
    gap = int(gap_input)

    # Set score matrix ...
    score_matrix = create_score_matrix(score_matrix_input)  


    # ----- What we see in the terminal
    print("Beep boop")
    D = sp_exact_3(seq1, seq2, seq3, score_matrix, gap)
    print(get_cost_3(D))
    print("")