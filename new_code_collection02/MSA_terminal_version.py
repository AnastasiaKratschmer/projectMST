import Bio
import numpy as np
from Bio import SeqIO
import sys
import os
import random as random
from functions_collection_copy import *

# ------------------------------ ACTUAL PROGRAM ------------------------------
if __name__ == "__main__":
    # ------ Parameter setup
    # Read arguments ...
    if not len(sys.argv) == 6:
            sys.stderr.write("USAGE: python3 %s < FASTA file of sequences > "
                            "< gap cost >  < score matrix > < integrity check?(True or False) > < type >\n" % sys.argv[0])
            sys.exit(1)
    seqs_input, gap_input, score_matrix_input,check,type = sys.argv[1:]

    # Set sequences (extract from file if necessary) ...
    seqs, names = parse_fasta_multiple_remove_n(seqs_input)
    
    # Set gap and score matrix ...
    gap = int(gap_input)

    # Set score matrix ...
    score_matrix = create_score_matrix(score_matrix_input)

    # ----- What we see in the terminal
    print("It's running!!\n")
    print("Computing the approximate cost of aligning the " + str(len(seqs)) + " sequences...")
    if type=="k":
        matrix, MSA_list, total_cost, in_which_MSA_is_t = new_assembly_gradual_x(seqs, score_matrix, gap,check_integrity=check)
    elif type=="p":
        matrix, MSA_list, total_cost, in_which_MSA_is_t = new_assembly_Prim_x(seqs, score_matrix, gap,check_integrity=check)
    elif type=="g":
        matrix, MSA_list, total_cost, in_which_MSA_is_t = new_assembly_Gus_x(seqs, score_matrix, gap,check_integrity=check)
    print("Done!\nTotal cost was:" )
    print(total_cost)

    save_output = input("Do you want to save the output? (Y/n): ").strip().lower()

    if save_output == 'y':
        # Ask the user for a filename
        filename = input("Enter a filename: ")
   