import Bio
import numpy as np
from Bio import SeqIO
import sys
import csv
import os
import random as random
from functions_collection_copy import *

# ------------------------------ ACTUAL PROGRAM ------------------------------
if __name__ == "__main__":
    # ------ Parameter setup
    # Read arguments ...
    if not len(sys.argv) == 6:
            sys.stderr.write("USAGE: python3 %s < FASTA file of sequences > "
                            "< gap cost >  < score matrix > < integrity check?(True or False) > < type (k,p or g for Kruskal, Prim or Gusfield) >\n" % sys.argv[0])
            sys.exit(1)
    seqs_input, gap_input, score_matrix_input,check,type = sys.argv[1:]
    check=check.lower() == 'true'

    # Set sequences (extract from file if necessary) ...
    seqs, names = parse_fasta_multiple_remove_n(seqs_input)
    
    # Set gap and score matrix ...
    gap = int(gap_input)

    # Set score matrix ...
    score_matrix = create_score_matrix(score_matrix_input)

    # ----- What we see in the terminal
    print("It's running!!\n")
    print("Computing the approximate cost of aligning the " + str(len(seqs)) + " sequences...")
    print("check is:"+ str(check))
    if type=="k":
        matrix, MSA_list, total_cost, in_which_MSA_is_it = new_assembly_gradual_x(seqs, score_matrix, gap,check_integrity=check)
        names=make_names_list(in_which_MSA_is_it)
        strings=make_str_from_cols(MSA_list[0])

    elif type=="p":
        matrix, MSA_list, total_cost, in_which_MSA_is_it = new_assembly_Prim_x(seqs, score_matrix, gap,check_integrity=check)
        names=make_names_list(in_which_MSA_is_it)
        strings=make_str_from_cols(MSA_list[0])
    elif type=="g":
        matrix, MSA_list, total_cost, names = new_assembly_Gus_x(seqs, score_matrix, gap,check_integrity=check)
        strings=make_str_from_cols(MSA_list)
    else:
         print("Please choose either k (Kruskal),p (Prim) or g (Gusfield) as type")
         sys.exit()
         
    print("Done!\nTotal cost was:" )
    print(total_cost)
    print("And your alignment is:")
    print(strings)
    

    save_output = input("Do you want to save the output? (Y/n): ").strip().lower()

    if save_output == 'y':
        # Ask the user for a filename
        filename = input("Enter a filename: ")
        csv_file_path=filename+".csv"
        with open(csv_file_path, 'a', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(["type: "+type])
                csvwriter.writerow(["nr of strings: "+str(len(seqs))])
                csvwriter.writerow(["gap cost: "+str(gap)])
                csvwriter.writerow(["total cost: "+str(total_cost)])
                for i, string in enumerate(strings):
                    csvwriter.writerow([f"{names[i]} {string}"])
                csvwriter.writerow([])

   