import numpy as np
import sys
import os
from tqdm import tqdm # loading bar
from helpers.utils_copy import linear_C, get_cost_2, get_sequence_string, parse_fasta_multiple, create_score_matrix, write_alignments_to_file, linear_backtrack, extend_alignment, Tree, convert_to_desired_format2,convert_format_mat_to_pseudomat,find_min_span_edges,get_visiting_order


matrix_for_MST=[]

def sp_approx(seqs: list[str], score_matrix: dict, gap_cost: int, verbose=False, return_center_string=False):
    # STEP 1: Find the center string, s1
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # loop over all distinct pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
              matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    print(matrix)
    matrix_for_MST=matrix 
    matrix_for_MST=convert_to_desired_format2(matrix_for_MST)
    min_span_edges=find_min_span_edges(matrix_for_MST)
    print(matrix_for_MST)
    print(min_span_edges)
    visiting_order=get_visiting_order(min_span_edges,"A") #A needs to be substituted at some point
    print(visiting_order) #visiting order is now letters, but we would need that as numbers/idices from the score matrix to keep track. #APPARENTLY NOT ACTUALLY IN USE

    

    # find center string/guide 
    s1_idx = np.argmin(matrix.sum(axis = 1))
    
    s1 = seqs[s1_idx]
    seqs.insert(0, seqs.pop(s1_idx)) # move guide to front of list
    if verbose: print("The center string, s1, is sequence no." + str(s1_idx+1)) # just a print statement to see which string is the center string

    # STEP 2: Construct alignment M
    M: list[list[str]] = [[letter] for letter in [*s1]]
    cost_list = []
    # print("first M = \n" + str(M))
    for i in range(1, len(seqs)):
        # if i == s1_idx: # skip the guide 
        #    continue
        cost = linear_C(gap_cost, score_matrix, s1, seqs[i])
        cost_list.append(get_cost_2(cost))
        
        # prepare A-matrix for extension
        alignment1_str, alignment2_str = linear_backtrack(s1, seqs[i], cost, score_matrix, gap_cost)
        alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
        A = [list(e) for e in zip(alignment1,alignment2)]
        
        # extend
        Mk = extend_alignment(M, A)
        M = Mk
    
    # ACTUALLY COMPUTE (approximate) COST
    total_cost = compute_cost(M, score_matrix, gap_cost)
    
    if return_center_string: return total_cost, M, s1_idx
    return total_cost, M, matrix_for_MST, visiting_order

def compute_cost(M, score_matrix, gap_cost):
    # the cost of the alignment is the sum of the cost of each column
    cost_of_columns = map(lambda m: sum_of_column(m, score_matrix, gap_cost), M)
    return sum(cost_of_columns)

def sum_of_column(col: list[str], score_matrix: dict, gap: int):
    k = len(col)
    cost = 0
    # loop through unique pairs/combinations, sum the cost alone the way
    for i in range(0,k):
        for j in range(i+1,k):
            # if -, -: add nothing to cost
            if col[i] == '-' and col[j] == '-':
                 j = j+1
            # if letter, -: add gap to cost
            # the ^ is apparently an exclusive or... a normal or should also be fine, bc we have just covered the double sitch above!
            if col[i] == '-' or col[j] == '-': #maybe change this ^to or, 'because that made it work for me
                 cost = cost + gap
            # if letter, letter: add subst (look up in score_matrix) to cost
            if col[i] != '-' and col[j] != '-':
                 cost = cost + score_matrix[col[i]][col[j]]  
    return cost


# ------------------------------ ACTUAL PROGRAM ------------------------------
if __name__ == "__main__":
    # ------ Parameter setup
    # Read arguments ...
    if not len(sys.argv) == 4:
            sys.stderr.write("USAGE: python3 %s < FASTA file of sequences > "
                            "< gap cost >  < score matrix >\n" % sys.argv[0])
            sys.exit(1)
    seqs_input, gap_input, score_matrix_input = sys.argv[1:]

    # Set sequences (extract from file if necessary) ...
    seqs, names = parse_fasta_multiple(seqs_input)
    
    # Set gap and score matrix ...
    gap = int(gap_input)

    # Set score matrix ...
    score_matrix = create_score_matrix(score_matrix_input)  

    # ----- What we see in the terminal
    print("Beep boop!\n")
    print("Computing the approximate cost of aligning the " + str(len(seqs)) + " sequences...")
    cost, M,matrix_for_MST,visiting_order = sp_approx(seqs, score_matrix, gap)
    print("vis ord:"+str(visiting_order))
    print("Done!\n")
    print("Cost: " + str(cost))
    #print(score_matrix)
    print(M)

    print()