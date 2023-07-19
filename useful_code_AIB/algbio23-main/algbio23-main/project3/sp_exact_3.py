import numpy as np
import sys
import os
from tqdm import tqdm  # loading bar
from helpers.utils import get_sequence_string, parse_fasta, create_score_matrix, write_alignments_to_file

def sp_exact_3(seq1: str, seq2: str, seq3: str, score_matrix: dict, gap_cost: int):
    n1 = len(seq1)
    n2 = len(seq2)
    n3 = len(seq3)

    D = np.full((n1+1, n2+1, n3+1), np.infty)
    for i in range(0, n1+1):
        for j in range(0, n2+1):
            for k in range(0, n3+1):
                d0, d1, d2, d3, d4, d5, d6, d7 = float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf')
                cij = score_matrix[seq1[i-1]][seq2[j-1]]
                cik = score_matrix[seq1[i-1]][seq3[k-1]]
                cjk = score_matrix[seq2[j-1]][seq3[k-1]]
                if i == 0 and j == 0 and k == 0:
                    d0 = 0
                if i > 0 and j > 0 and k > 0:
                    d1 = D[i - 1, j - 1, k - 1] + cij + cik + cjk #case 1
                if i > 0 and j > 0 and k >= 0:
                    d2 = D[i - 1, j - 1, k] + cij + 2 * gap_cost #case 2
                if i > 0 and j >= 0 and k>0:
                    d3 = D[i - 1, j, k - 1] + cik + 2 * gap_cost #case 3
                if i >= 0 and j > 0 and k > 0:
                    d4 = D[i, j - 1, k - 1] + cjk + 2 * gap_cost #case 4
                if i > 0 and j >= 0 and k >= 0:
                    d5 = D[i - 1, j, k] + 2 * gap_cost #case 5
                if i >= 0 and j > 0 and k >= 0:
                    d6 = D[i, j - 1, k] + 2 * gap_cost #case 6
                if i >= 0 and j >= 0 and k > 0:
                    d7 = D[i, j, k - 1] + 2 * gap_cost
                D[i, j, k] = min(d0, d1, d2, d3, d4, d5, d6, d7)
    return D

def backtrack_3seq(seq1:str,seq2:str,seq3:str, D, score_matrix:dict,gap_cost:int):
    i = len(seq1)
    j = len(seq2)
    k = len(seq3)

    alignment_of_1=[]
    alignment_of_2=[]
    alignment_of_3=[]
    
    while j>=0 or i>=0 or k>=0:

        cij = score_matrix[seq1[i-1]][seq2[j-1]]
        cik = score_matrix[seq1[i-1]][seq3[k-1]]
        cjk = score_matrix[seq2[j-1]][seq3[k-1]]

        if i>0 and j>0 and k>0 and D[i,j,k]==D[i-1,j-1,k-1]+ cij + cik + cjk: #case 1 
            alignment_of_1.append(seq1[i-1])
            alignment_of_2.append(seq2[j-1])
            alignment_of_3.append(seq3[k-1])
            i -= 1
            j -= 1
            k -= 1

        elif i>0 and j>0 and D[i,j,k]==D[i - 1, j - 1, k] + cij + 2 * gap_cost: #case 2
            alignment_of_1.append(seq1[i-1])
            alignment_of_2.append(seq2[j-1])
            alignment_of_3.append('-')
            i -= 1
            j -= 1
            
        elif i>0 and k>0 and D[i,j,k]==D[i - 1, j, k-1] + cik + 2 * gap_cost:#case 3
            alignment_of_1.append(seq1[i-1])
            alignment_of_2.append('-')
            alignment_of_3.append(seq1[k-1])
            i -= 1
            k -= 1
            
        elif j>0 and k>0 and D[i,j,k]==D[i, j - 1, k - 1] + cjk + 2 * gap_cost: #case 4
            alignment_of_1.append('-')
            alignment_of_2.append(seq2[j-1])
            alignment_of_3.append(seq3[k-1])
            j -= 1
            k -= 1

        elif i>0 and D[i,j,k]==D[i - 1, j, k] + 2 * gap_cost: #case 5
            alignment_of_1.append(seq1[i-1])
            alignment_of_2.append('-')
            alignment_of_3.append('-')
            i -= 1

        elif j>0 and D[i,j,k]==D[i, j - 1, k] + 2 * gap_cost:
            alignment_of_1.append('-')
            alignment_of_2.append(seq2[j-1])
            alignment_of_3.append('-')
            j -= 1

        elif k>0 and D[i,j,k]==D[i, j, k - 1] + 2 * gap_cost: #case 7
            alignment_of_1.append('-')
            alignment_of_2.append('-')
            alignment_of_3.append(seq3[k-1])
            k -= 1

        else:
            break
    
    alignment_of_1.reverse()
    alignment_of_2.reverse()
    alignment_of_3.reverse()
    
    return ''.join(alignment_of_1), ''.join(alignment_of_2), ''.join(alignment_of_3)

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
    print("Beep boop!\n")
    print("Computing cost of alignment...")
    D = sp_exact_3(seq1, seq2, seq3, score_matrix, gap)
    print("Backtracking...")
    alignments = backtrack_3seq(seq1,seq2,seq3, D, score_matrix, gap)
    print("Done!\n ")

    print("Cost:\n", D[-1, -1, -1])
    print("Alignments:\n", alignments[0], "\n", alignments[1], "\n", alignments[2])
    print()
