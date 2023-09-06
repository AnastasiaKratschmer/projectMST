import Bio
import numpy as np
import sys
import os
import networkx as nx
import random as random
from tqdm import tqdm # loading bar
from utils_copy import linear_C, get_cost_2, get_sequence_string, parse_fasta_multiple, create_score_matrix, write_alignments_to_file, linear_backtrack, fill_graph,new_sp_approxi_combi
from utils_copy import convert_to_desired_format_nr_version, compute_cost, my_traversal_simply, extend_alignment_chaos, find_min_span_edges_testing, parse_fasta_multiple_remove_n
import matplotlib.pyplot as plt

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
    seqs, names = parse_fasta_multiple_remove_n(seqs_input)
    
    # Set gap and score matrix ...
    gap = int(gap_input)

    # Set score matrix ...
    score_matrix = create_score_matrix(score_matrix_input)  

    # ----- What we see in the terminal
    print("It's running!!\n")
    print("Computing the approximate cost of aligning the " + str(len(seqs)) + " sequences...")
    cost, M, matrix_for_MST, G = new_sp_approxi_combi(seqs, score_matrix, gap, verbose=True)
    print("Done!\n")
    #print("Cost: " + str(cost))
    print()
    #print(M)
    #nx.draw(G, with_labels=True)
    plt.show(block=True)