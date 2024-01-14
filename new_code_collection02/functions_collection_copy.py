import random as random
import Bio
import numpy as np
import sys
import os
#import networkx as nx
from Bio import SeqIO
import timeit

# Basic functions

def parse_fasta(filename):
    seq = ""
    for record in SeqIO.parse(filename, "fasta"):
        seq = record.seq.lower()
    return seq

def parse_fasta_multiple(filename):
    seq = []
    name = []
    for record in SeqIO.parse(filename, "fasta"):
        seq.append(record.seq.lower())
        #print(record.name)
        name.append(record.name)
    return seq, name


def parse_fasta_multiple_remove_n(filename): #the parse fasta function that removes ns and puts a random base instead :) 
    seq = []
    name = []
    bases=['a', 'c', 't', 'g', 'A','C','G','T']
    for record in SeqIO.parse(filename, "fasta"):
        new_seq = []  # Initialize new_seq for each record
        for letter in record.seq:
            if letter not in bases:
                random_base=random.choice(['a', 'c', 't', 'g'])
                new_seq.append(random_base)
            else:
                new_seq.append(letter)
        new_seq = ''.join(new_seq)
        new_seq = new_seq.lower()
        seq.append(new_seq.lower())  # Append the modified sequence to the seq list
        name.append(record.name)
    return seq, name

def write_alignments_to_file(alignment1, alignment2, filename):
    with open(filename, 'w') as f_out:
        f_out.write('>' + 'seq1' + '\n' + alignment1 + '\n\n' +
                    '>' + 'seq2' + '\n' + alignment2)

def get_sequence_string(input: str):
    if os.path.isfile(input):
        return parse_fasta(input).lower()
    else:
        return input.lower()

def create_score_matrix(filename: str, verbose=False):
    score_matrix = {}
    with open(filename, 'r') as fh:
        alphabet_size = fh.readline()
        if verbose: print("Alphabet size: ", alphabet_size)
        alphabet = []
        lines = fh.readlines()
        for line in lines:
            alphabet.append(line[0].lower())
        if verbose: print("Alphabet: ", alphabet)
        for letter, line in zip(alphabet, lines):
            scores = line.split("\t") # N.B.: tabs are what split the characters! 
            scores_for_this_row = {}
            for l,score in zip(alphabet, scores[1:]): 
                scores_for_this_row[l] = int(score)
            score_matrix[letter] = scores_for_this_row
        if verbose: print("Score matrix: ", score_matrix)
    return score_matrix

def linear_C(gap, score_matrix, seq1, seq2, verbose=False):
    n = len(seq1)
    m = len(seq2)
    alphabet = list(score_matrix.keys())

    matrix = np.empty((n + 1, m + 1))
    # Setting the cells to not a number instead of zero in order to avoid any "mishaps"
    matrix[:] = np.nan

    # handle first column and row separately, where
    # only the neighbour cell can influence the cost
    matrix[0,0] = 0
    for i in range(1, n+1): # don't start at 0,0, it should remain 0
        matrix[i,0] = matrix[i-1,0] + gap
    
    for j in range(1, m+1): # don't start at 0,0, it should remain 0
        matrix[0,j] = matrix[0,j-1] + gap
    
    # now fill out the table from 1,1
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            if not np.isnan(matrix[i,j]): 
                # we already computed that entry
                continue

            if (not seq1[i-1] in alphabet):
                raise Exception(seq1[i-1] + " is an invalid character in sequence 1.")
            if (not seq2[j-1] in alphabet):
                raise Exception(seq2[j-1] + " is an invalid character in sequence 2.")
            
            v1 = matrix[i-1, j-1] + score_matrix[seq1[i-1]][seq2[j-1]]
            v2 = matrix[i-1, j] + gap
            v3 = matrix[i, j-1] + gap
            matrix[i,j] = min(v1, v2, v3)#, 0)
    return(matrix)

def get_cost_2(cost_matrix):
    return cost_matrix[-1,-1]

def linear_backtrack(A: str, B: str, T, score_matrix, gap: int, verbose=False):
    i = len(A) 
    j = len(B)

    alignment_of_A = []
    alignment_of_B = []
    while ( i >= 0 or j >= 0 ):
        if (i > 0) and (j > 0) and T[i,j] == (T[i-1,j-1] + score_matrix[A[i-1]][B[j-1]]):
            # “output column (A[i], B[j])”
            alignment_of_A.append(A[i-1])
            alignment_of_B.append(B[j-1])
            i -= 1
            j -= 1
        elif (i > 0) and T[i,j] == (T[i-1,j] + gap):
            #“output column(A[i], -)”
            alignment_of_A.append(A[i-1])
            alignment_of_B.append('-')
            i -= 1
        elif (j > 0) and T[i,j] == (T[i,j-1] + gap):
            #“output column(-, B[j])”
            alignment_of_A.append('-')
            alignment_of_B.append(B[j-1])
            j -= 1
        else:
            break

    # flip them because we were backtracking, so they're written backwards
    alignment_of_A.reverse()
    alignment_of_B.reverse()

    if verbose: print("Alignments: \n" + ''.join(alignment_of_A) + "\n" + ''.join(alignment_of_B))
    
    return ''.join(alignment_of_A), ''.join(alignment_of_B)

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
            elif col[i] == '-' or col[j] == '-':
                 cost = cost + gap
            # if letter, letter: add subst (look up in score_matrix) to cost
            elif col[i] != '-' and col[j] != '-':
                 cost = cost + score_matrix[col[i]][col[j]]  
    return cost

def convert_to_desired_format_nr_version(distance_matrix): #makes the distance matrix into the pseudomatrix format needed for selecting edges int the func. find_min_span_edges_testing
    # Get the number of nodes in the distance matrix
    num_nodes = len(distance_matrix)
    # Initialize an empty list to store the rows of the new format
    matrixx_rows = []

    # Process the distance_matrix to generate the rows of the new format
    for i in range(num_nodes):
        for j in range(i+1, num_nodes):  # Avoid duplicates and self-distances
            distance = int(distance_matrix[i, j])
            node1 = str(i)  # Node names as 1,2,3 ... could maybe be given better names, for instance the names from the fasta file..
            node2 = str(j)
            matrixx_rows.append(["", distance, node1, node2]) #the zero'th col is empty now but is for marking if the edge is included in the MST

    # Convert the list of rows to a np array
    matrixx_np = np.array(matrixx_rows)
    return matrixx_np

def find_min_span_edges_testing(pseudomatrix, verbose=False): #the function actually implementing Kruskal's algo
    sorted_indices = np.lexsort((pseudomatrix[:, 1].astype(int),)) #sort edges to have the shorter ones first.
    E = pseudomatrix[sorted_indices]
    if verbose: 
        print("this is E (sorted matrix without any stars yet): ")
        print(E)
    names= np.unique(pseudomatrix[:, 2:]) #extracting all node names to keep track of them in name_dict.
    name_dict= {letter:i for i, letter in enumerate(names)}
    if verbose:
        print("the names of nodes going into the first name dict are: "+str(names)+" and the name_dict is orginally "+ str(name_dict))
    E[:,0]='' #always start proces by setting zero'th col as empy. May be redundant, but just to be sure ;) 
    x=0 #to keep track of current row
    it=0 # iteration number for print statement...
    while len(set(name_dict.values()))>1: #so while all nodes are still not unified in one tree..
        it+=1 #only there for print statement below :) 
        if verbose: print("\n \n this is it "+str(it)+" of the find_min_span_edges_testing-func.") #the aforementioned print statement!
        min0,min1,min2=E[x][1],E[x][2],E[x][3] #take the next row in the pseudomatrix to get the len of the edge and the names of the two nodes it connects
        tree1_id= name_dict[min1] #getting the positions of the nodes currently in question in the current list of trees
        tree2_id= name_dict[min2]
        if tree2_id == tree1_id: #if the two nodes are already in the same tree, skip this row.
            if verbose: print("the two nodes are already in the same tree")
            x += 1
            continue
        else:
            E[x][0] = "*" #if they were not already in the same tree, then mark that the row/edge is going to be used!
            #generally making sure that the new tree assigned to the node is the one with the lowest number. Also updating name_dict to keep the name interval 'closed'.
            if tree1_id < tree2_id:
                orig_tree2_id = tree2_id
                for key, value in name_dict.items():
                    if value == orig_tree2_id:
                        name_dict[key] = tree1_id
                    elif value > orig_tree2_id:
                        name_dict[key] -= 1
            else:
                orig_tree1_id = tree1_id
                for key, value in name_dict.items():
                    if value == orig_tree1_id:
                        name_dict[key] = tree2_id
                    elif value > orig_tree1_id:
                        name_dict[key] -= 1
            x += 1
            if verbose: print("after that iteration, we end up with this dict: "+ str(name_dict)+ "and the set is: "+ str(set(name_dict.values()))+" and the len of that set is "+ str(len(set(name_dict.values()))))
    if verbose: 
        print ("Here are the edges included in the MST, marked with a star! \n")
        print(E)
            
    res_mat=E 
    return res_mat

def find_highest_second_element(d, first_element): #to find the highest col nr already in the MSA
    highest_second_element = 0
    for key, value in d.items():
        if value[0] == first_element:
            second_element = value[1]
            if second_element > highest_second_element:
                highest_second_element = second_element
    return highest_second_element

def finish_running_through_MSA1(MSA1,MSA2,j,i,l):
    how_many_cols_in_MSA2=len(MSA2[0])
    gaps_to_add=['-']*how_many_cols_in_MSA2
    new_col=list(MSA1[j])+gaps_to_add
    return new_col

def finish_running_through_MSA2(MSA1,MSA2,j,i,l):
    how_many_cols_in_MSA1=len(MSA1[0])
    gaps_to_prepend=['-']*how_many_cols_in_MSA1
    new_col=gaps_to_prepend+list(MSA2[l])
    return new_col

def equality_here_test(a,b):
    equality=False
    if a=='-' and b=='-' or a!='-' and b!='-':
        equality=True
    return equality

def new_identify_merge_case_compact(guide_PA,MSA1,MSA2,node1,node2,MSA_position_overview_dict,i,j,l):
    node1_position=MSA_position_overview_dict[node1]
    node1_col=node1_position[1]
    node2_position=MSA_position_overview_dict[node2]
    node2_col=node2_position[1]
    case=None
    move_on_in_guide_flag=False
    a=guide_PA[i][0]
    b=guide_PA[i][1]
    c= MSA1[j][node1_col]
    d=MSA2[l][node2_col]
    look_forward_flag1=False
    look_forward_flag2=False

    if equality_here_test(a,c)==True and equality_here_test(b,d)==True: case=3 #a==c and b==d
    elif equality_here_test(a,c)==False and equality_here_test(b,d)==False:
        if equality_here_test(b,c)==True:
            if guide_PA[i][0]!='-':
                k=1
                while MSA1[j+k][node1_col]=='-' and guide_PA[i][0]!='-':
                    k+=1
                if MSA1[j+k][node1_col]==guide_PA[i][0]:
                    look_forward_flag1=True
            if guide_PA[i][1]!='-':
                h=1
                while MSA2[l+h][node2_col]=='-' and guide_PA[i][1]!='-':
                    h+=1
                if MSA2[l+h][node2_col]==guide_PA[i][1]:
                    look_forward_flag2=True
        case=4
    elif equality_here_test(a,c)==True:#a==c
        if d!='-':
            case= 2 #was orginially 1 2 omg just seeing what will happen if I put it to 1. #truing to reverse to see
            move_on_in_guide_flag=True
        else: case= 1
    elif equality_here_test(a,c)==False: #a!=c
        if a=='-':
            case= 1 #was originally 1 omg!, just checking what will happen if I put it to 2 #reversing to see what will happen...
            move_on_in_guide_flag=True
        else: case=2
    else: sys.exit()
    guide_was=guide_PA[i]
    MSA2_was=MSA2[l]
    return case,move_on_in_guide_flag,guide_was,MSA2_was, look_forward_flag1, look_forward_flag2

def alt_handle_case1(MSA1,MSA2,i,j,l):
    how_many_cols_in_MSA1=len(MSA1[0])
    gaps_to_prepend=['-']*how_many_cols_in_MSA1
    new_col=gaps_to_prepend+list(MSA2[l])
    return new_col #after this, increase only l

def alt_handle_case2(MSA1,MSA2,i,j,l):
    how_many_cols_in_MSA2=len(MSA2[0])
    gaps_to_append=['-']*how_many_cols_in_MSA2
    col_now=MSA1[j]
    col_new=list(MSA1[j])+gaps_to_append
    return col_new #after this, increase only j

def alt_handle_case3(MSA1,MSA2,i,j,l):
    col_new=list(MSA1[j])+list(MSA2[l])
    return col_new #after this, increase i, j and l

def alt_handle_case4(MSA1,MSA2,i,j,l):
    #col_new=list(MSA1[j])+list(MSA2[l])
    how_many_cols_in_MSA1=len(MSA1[0])
    how_many_cols_in_MSA2=len(MSA2[0])
    #make first column!
    first_col_gaps=['-']*how_many_cols_in_MSA2
    first_col=list(MSA1[j])+first_col_gaps
    #make second colum!
    second_col_gaps=['-']*how_many_cols_in_MSA1
    second_col=second_col_gaps+list(MSA2[l])
    return first_col,second_col #increade l and j after 

def alt_alt_merge_united(guide_PA,MSA_list,in_which_MSA_is_it,node1,node2):
    #print(type(in_which_MSA_is_it))
    #print(in_which_MSA_is_it)
    #print(type(node1))
    #print(type(node2))
    MSA_of_node1=in_which_MSA_is_it[node1][0]
    MSA_of_node2=in_which_MSA_is_it[node2][0]
    spot_of_node1=in_which_MSA_is_it[node1][1]
    spot_of_node2=in_which_MSA_is_it[node2][1]
    
    #print("the node1 is right now in the alignment nr "+str(MSA_of_node1))
    #print("the node2 is right now in the alignment nr "+str(MSA_of_node2))
    MSA1=MSA_list[int(MSA_of_node1)]
    MSA2=MSA_list[int(MSA_of_node2)]

    
    j=0 
    i=0
    l=0
    MSA_new=[]
    while i<=(len(guide_PA)-1) and j<=(len(MSA1)-1) and l<=(len(MSA2)-1):
        new_col=[]
        #print("I don't have the case identified yet, but MSA1 is: " + str(MSA1[j]) + " MSA2 is: " + str(MSA2[l]) + " ,and the guide is: " + str(guide_PA[i]))
        #print("nodes pointing to cols: " + str(spot_of_node1) + " , " + str(spot_of_node2))
        #print("i,j,l: "+str(i)+','+str(j)+','+str(l))
        #print("max for those should be: "+str(len(guide_PA))+','+str(len(MSA1))+','+str(len(MSA2)))
        case,move_on_in_guide_flag,guide_was,MSA2_was,look_forward_flag1, look_forward_flag2=new_identify_merge_case_compact(guide_PA,MSA1,MSA2,node1,node2, in_which_MSA_is_it,i,j,l)
        #print("flag: "+str(move_on_in_guide_flag))
        #print("case:"+str(case))
        #print("guide_was: "+str(guide_was))
        if case==1:
            col=alt_handle_case1(MSA1,MSA2,i,j,l)
            MSA_new.append(col) 
            if move_on_in_guide_flag==True:
                #print("the damn flag was true!")
                i+=1
            l+=1
            #print("I just made the new col: "+ str(col))
            #print("MSA_new now has len: "+ str(len(MSA_new)))
            new_col.append(col)
        if case==2:
            col=alt_handle_case2(MSA1,MSA2,i,j,l)
            MSA_new.append(col)
            if move_on_in_guide_flag==True:
                #print("the flag was true dammit")
                i+=1
            j+=1
            #print("I just made the new col: "+ str(col))
            #print("MSA_new now has len: "+ str(len(MSA_new)))
            new_col.append(col)
        if case==3:
            col=alt_handle_case3(MSA1,MSA2,i,j,l)
            MSA_new.append(col)
            i+=1
            j+=1
            l+=1
            #print("I just made the new col: "+ str(col))
            #print("MSA_new now has len: "+ str(len(MSA_new)))
            new_col.append(col)
        if case==4:
            col1,col2=alt_handle_case4(MSA1,MSA2,i,j,l)
            if look_forward_flag1==True:
                MSA_new.append(col1)
                new_col.append(col1)
                j+=1
            elif look_forward_flag2==True:
                MSA_new.append(col2)
                new_col.append(col2)
                l+=1
            else:
                MSA_new.append(col1)
                MSA_new.append(col2)
                j+=1
                l+=1
                #print("I just made the new col: "+ str(col1))
                #print("I just made the new col: "+ str(col2))
                #print("MSA_new now has len: "+ str(len(MSA_new)))
                new_col.append(col1)
                new_col.append(col2)
        for object in new_col:
            #print(object)
            second_element_in_merged_pair=object[(-(len(MSA2_was))+spot_of_node2)]
            selected_stuff=[object[int(spot_of_node1)],second_element_in_merged_pair]
            #print("selected stuff from col made: ")
            #print(selected_stuff)
            if case==1 and move_on_in_guide_flag==False or case==2 and move_on_in_guide_flag==False or case==4:
                if guide_was==selected_stuff: print("ALERT!")
            if guide_was!=selected_stuff and selected_stuff!=['-', '-']: print("ALARM!")
            if case==1 or case==2 and move_on_in_guide_flag==True:
                if guide_was!=selected_stuff and selected_stuff!=['-', '-']: print("WARNING! i is gonna be increased, maybe unrightfully")
            
    if j <= (len(MSA1)-1):
        while j<= (len(MSA1)-1):
            column=finish_running_through_MSA1(MSA1,MSA2,j,i,l)
            j+=1
            #l+=1
            MSA_new.append(column)
            #print("I just made the new col: "+ str(column))
            #print("MSA_new now has len: "+ str(len(MSA_new)))
    if l<= (len(MSA2)-1):
        while l<= (len(MSA2)-1):
            column=finish_running_through_MSA2(MSA1,MSA2,j,i,l)
            l+=1
            MSA_new.append(column)
            #print("I just made the new col: "+ str(column))
            #print("MSA_new now has len: "+ str(len(MSA_new)))
    if i<=(len(guide_PA)):
        i+=1
    else:
        print('yikes, the strings have run out but the guide has not')
    print("\n\n")
    return MSA_new

def integrity_check_OBO_and_gradual(seqs, in_which_MSA_is_it, who_aligned_to_who, MSA_list,matrix, score_matrix, gap_cost):
    for i,seq in enumerate(seqs):
        col_to_extract=in_which_MSA_is_it[str(i)][1]
        j=0
        new_str_with_gaps=[]
        new_str_no_gaps=[]
        while j<=len(MSA_list[0])-1:
            found=MSA_list[0][j][col_to_extract]
            j+=1
            new_str_with_gaps.append(found)
            if found !='-':
                new_str_no_gaps.append(found)
        new_str_no_gaps=''.join(new_str_no_gaps)
        new_str_with_gaps=''.join(new_str_with_gaps)
        if new_str_no_gaps==seq:
            print("integrity check 1 passed for seq "+str(i)+": string same as original sequence")
        else:
            print("Yikes, integrity check 1 did not pas for seq "+str(i)+". constrast( new, orig): \n"+str(new_str_no_gaps)+"\n"+str(seq))
            sys.exit()
        #print("structure of who_aligned_to_who: ")
        #print(who_aligned_to_who)
    for element in who_aligned_to_who:
        seq1_nr=element[0]
        seq2_nr=element[1]
        #print("seq1_nr and seq2_nr are: "+str(seq1_nr)+" , "+str(seq2_nr))
        pos_in_MSA_seq1=in_which_MSA_is_it[seq1_nr][1]
        pos_in_MSA_seq2=in_which_MSA_is_it[seq2_nr][1]
        #print("pos_in_MSA_seq1,pos_in_MSA_seq2: "+str(pos_in_MSA_seq1)+" , "+str(pos_in_MSA_seq2))
        seq1_from_MSA=[]
        seq2_from_MSA=[]
        j=0
        while j<=len(MSA_list[0])-1:
            found=MSA_list[0][j][pos_in_MSA_seq1]
            seq1_from_MSA.append(found)
            j+=1
        j=0
        while j<=len(MSA_list[0])-1:
            found=MSA_list[0][j][pos_in_MSA_seq2]
            seq2_from_MSA.append(found)
            j+=1
        union=[]
        k=0
        len_max=max(len(seq1_from_MSA),len(seq1_from_MSA))
        while k<=(len_max-1):
            el1=seq1_from_MSA[k]
            el2=seq2_from_MSA[k]
            tuple_like_zip=[el1,el2]
            union.append(tuple_like_zip)
            k+=1
        #print("union of the two after merge looks like: "+str(union))
        cost_after_MSA=compute_cost(union,score_matrix,gap_cost)
        if cost_after_MSA==matrix[int(seq1_nr)][int(seq2_nr)]:
            print("integrity test 2 passed for: "+str(seq1_nr)+" and "+ str(seq2_nr)+": alignment cost consistent with guide alignment")
        else:
            print("Yikes, integrity check 2 did not pas for: "+str(seq1_nr)+" and "+ str(seq2_nr))
            print("Costs were before and after:"+str(matrix[0][int(seq1_nr)][int(seq2_nr)])+" and "+str(cost_after_MSA))
            cost_for_suppesed_to_have_been=linear_C(gap_cost,score_matrix,seqs[int(seq1_nr)],seqs[int(seq2_nr)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(seq1_nr)], seqs[int(seq2_nr)],cost_for_suppesed_to_have_been, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            should_have_been= [list(e) for e in zip(alignment1,alignment2)]
            print("should have been:"+ str(should_have_been))
            all_gaps_cols_removed_from_union=[sublist for sublist in union if not all(item == '-' for item in sublist)]
            print("is: "+str(all_gaps_cols_removed_from_union))
            h=0
            while h<=len(should_have_been)-1:
                if should_have_been[h]==all_gaps_cols_removed_from_union[h]:
                    h+=1
                else:
                    print("index of first error: " + str(h) + " out of approximately " + str(len(should_have_been)) + ". The cols are these (should have been, are): " + str(should_have_been[h]) + " and " + str(all_gaps_cols_removed_from_union[h]))
                    sys.exit()
                    h+=1
    #integrity check 3
    for col in MSA_list[0]:
        if all(element == '-' for element in col):
            print("integrity test 3 failed: empty columns in alignment")
            sys.exit()
    print("integrity test 3 passed: No empty columns in alignment")
    return('Pass')


def perform_updates_OBO(in_which_MSA_is_it, node1, node2,united_MSA_new, MSA_list ):
    which_spot_in_MSA_list_to_update=min(in_which_MSA_is_it[node1][0],in_which_MSA_is_it[node2][0])
    which_spot_in_MSA_list_to_remove=max(in_which_MSA_is_it[node1][0],in_which_MSA_is_it[node2][0])
    MSA_list[which_spot_in_MSA_list_to_update]=united_MSA_new
    MSA_list.pop(which_spot_in_MSA_list_to_remove)
    companions_to_update=[]
    if in_which_MSA_is_it[node1][0]<in_which_MSA_is_it[node2][0]:
        how_many_cols_already_in_MS1=find_highest_second_element(in_which_MSA_is_it,which_spot_in_MSA_list_to_update)
        for key, value in in_which_MSA_is_it.items():
            if value[0]==which_spot_in_MSA_list_to_remove:
                companions_to_update.append(key)
            if value[0]>which_spot_in_MSA_list_to_remove:
                value[0]=(value[0]-1)
        for companion in companions_to_update:
            #print(in_which_MSA_is_it[companion])
            in_which_MSA_is_it[companion][0]=which_spot_in_MSA_list_to_update
            col_of_element_in_old_MSA2=in_which_MSA_is_it[companion][1]
            in_which_MSA_is_it[companion][1]=(how_many_cols_already_in_MS1+col_of_element_in_old_MSA2+1)
    else:
        companions_to_update1=[]
        companions_to_update2=[]
        how_many_cols_already_in_MS2=find_highest_second_element(in_which_MSA_is_it,which_spot_in_MSA_list_to_remove) #get the highest nr in MSA2. need to only update dict[0] for MSA1 and only dict[1] for MSA2
        for key, value in in_which_MSA_is_it.items():
            if value[0]==which_spot_in_MSA_list_to_update: #these need their value[1] updated, because they were "pushed downwards" in their alignment
                companions_to_update1.append(key)
            if value[0]==which_spot_in_MSA_list_to_remove: #these need their value[0] updated, because they were moved to a new alignment.
                companions_to_update2.append(key)
            if value[0]>which_spot_in_MSA_list_to_remove: #these need -1 in their value[1]
                value[0]=(value[0]-1)
        for companion1 in companions_to_update1:
            in_which_MSA_is_it[companion1][1]=int(in_which_MSA_is_it[companion1][1])+int(how_many_cols_already_in_MS2)+1
        for companion2 in companions_to_update2:
            in_which_MSA_is_it[companion2][0]=which_spot_in_MSA_list_to_update
    return in_which_MSA_is_it, MSA_list

def new_assembly_OBO_x(seqs,score_matrix,gap_cost, check_integrity=False):
    # Make a matrix to hold pairwise alignment costs for all alignment combinations!
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # Loop over all pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    matrix_for_MST=matrix #copy the matrix, so that we ca keep the old matrix and make a changed version to the "pseudomatrix" version
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST) #making the "pseudomatrix"
    min_span_edges_res=find_min_span_edges_testing(matrix_for_MST)
    in_which_MSA_is_it={}
    names= np.unique(min_span_edges_res[:, 2:])
    in_which_MSA_is_it ={name: [int(name),0] for name in names}
    MSA_list=[[[char] for char in seq] for seq in seqs]
    who_aligned_to_who=[]
    k=0
    for row in min_span_edges_res:
        k+=1
        if row[0]=="*":
            node1=row[2]
            node2=row[3]
            who_aligned_to_who.append([node1,node2])

    traversal_order = []
    already_there_set = set()
    queue = []

    # Add the first element to the queue
    queue.append(who_aligned_to_who[0])

    while queue:
        element = queue.pop(0)  # Dequeue the first element from the queue
        node1, node2 = element

        if node1 in already_there_set:
            # Swap nodes if node1 is already in already_there_set
            node1, node2 = node2, node1

        if node1 not in already_there_set or node2 not in already_there_set:
            traversal_order.append([node1, node2])
            already_there_set.add(node1)
            already_there_set.add(node2)

            # Find edges connected to the last added node and enqueue them
            for edge in who_aligned_to_who:
                if edge != element and (edge[0] == node1 or edge[1] == node1 or edge[0] == node2 or edge[1] == node2):
                    queue.append(edge)

    for element in traversal_order:
            node1=element[0]
            node2=element[1]
            cost=linear_C(gap_cost,score_matrix,seqs[int(node1)],seqs[int(node2)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(node1)], seqs[int(node2)], cost, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            A = [list(e) for e in zip(alignment1,alignment2)]
            united_MSA_new=alt_alt_merge_united(A,MSA_list,in_which_MSA_is_it,node1,node2)
            in_which_MSA_is_it, MSA_list=perform_updates_OBO(in_which_MSA_is_it, node1, node2, united_MSA_new, MSA_list)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    if check_integrity==True:
       integrity_check_OBO_and_gradual(seqs,in_which_MSA_is_it,who_aligned_to_who,MSA_list, matrix,score_matrix,gap_cost)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    return(matrix,MSA_list, total_cost)


def integrity_check_Gus(seqs,M,score_matrix,gap_cost, matrix, s1_idx):
    print("Gus int check was called!!")
    #integrity check part 1, to check if each string is the same before and after, except for gaps.
    #print("Im just gonna start this int check by printint the strings/sequences dammit: "+ str(seqs))
    for i in range(len(seqs)):
        seq=seqs[i] #extract the the old nr i from seqs. now we need to recover what the new name of nr i is.
        j=0
        new_str_with_gaps=[]
        new_str_no_gaps=[]
        while j<=len(M)-1:
            found=M[j][i]
            j+=1
            new_str_with_gaps.append(found)
            if found !='-':
                new_str_no_gaps.append(found)
        new_str_no_gaps=''.join(new_str_no_gaps)
        new_str_with_gaps=''.join(new_str_with_gaps)
        #print("seq and new_str_no_gaps were: "+ str(seq)+", "+ str(new_str_no_gaps))
        if new_str_no_gaps==seq:
            print("integrity check 1 passed for seq "+str(i)+": string same as original sequence")
        else:
            print("Yikes, integrity check 1 did not pas for seq "+str(i)+". constrast( new, orig): \n"+str(new_str_no_gaps)+"\n"+str(seq))
            sys.exit()
     #part 2 lol, are the alignments preserved, expect for gaps???
    for i in range(1,len(seqs)):
        seq=seqs[i]
        seq1_nr=0
        seq2_nr=i
        seq1_from_MSA=[]
        seq2_from_MSA=[]
        j=0
        while j<=len(M)-1:
            found=M[j][0]
            seq1_from_MSA.append(found)
            j+=1
        j=0
        while j<=len(M)-1:
            found=M[j][i]
            seq2_from_MSA.append(found)
            j+=1
        union=[]
        k=0
        len_max=max(len(seq1_from_MSA),len(seq1_from_MSA))
        while k<=(len_max-1):
            el1=seq1_from_MSA[k]
            el2=seq2_from_MSA[k]
            tuple_like_zip=[el1,el2]
            union.append(tuple_like_zip)
            k+=1
        #print("union of the two after merge looks like: "+str(union))
        cost_after_MSA=compute_cost(union,score_matrix,gap_cost)

        cost_for_suppesed_to_have_been=linear_C(gap_cost,score_matrix,seqs[int(seq1_nr)],seqs[int(seq2_nr)])
        if cost_after_MSA==cost_for_suppesed_to_have_been[-1,-1]:
            print("integrity test 2 passed for: "+str(seq1_nr)+" and "+ str(seq2_nr)+": alignment cost consistent with guide alignment")
        else:
            print("Yikes, integrity check 2 did not pass for: "+str(seq1_nr)+" and "+ str(seq2_nr))
            #print("The big matrix with all alignment prices is here: ", str(matrix))
            #print("I'm just gonna print the big alignment too: "+ str(M))
            print("Costs were before and after:"+str(matrix[int(seq1_nr)][int(seq2_nr)])+" and "+str(cost_after_MSA))
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(seq1_nr)], seqs[int(seq2_nr)],cost_for_suppesed_to_have_been, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            should_have_been= [list(e) for e in zip(alignment1,alignment2)]
            print("should have been:"+ str(should_have_been))
            all_gaps_cols_removed_from_union=[sublist for sublist in union if not all(item == '-' for item in sublist)]
            print("is: "+str(all_gaps_cols_removed_from_union))
            h=0
            while h<=len(should_have_been)-1:
                if should_have_been[h]==all_gaps_cols_removed_from_union[h]:
                    h+=1
                else:
                    print("index of first error: " + str(h) + " out of approximately " + str(len(should_have_been)) + ". The cols are these (should have been, are): " + str(should_have_been[h]) + " and " + str(all_gaps_cols_removed_from_union[h]))
                    sys.exit()
                    h+=1
            sys.exit()
    #integrity check 3
    for col in M:
        if all(element == '-' for element in col):
            print("integrity test 3 failed: empty columns in alignment")
            sys.exit()
    print("integrity test 3 passed: No empty columns in alignment")
    return("Passed")

def extend_alignment(M, A):
    '''
    M = [['a','a','c'], ['-','t','t'], ['-','t','-'], ['c','c','c'], ['g','-','g'], ['t','t','a']]
    means
        a - - c g t
    M = a t t c - t 
        c t - c g a

    A = [['a','a'], ['c','c'], ['g','g'], ['-','g'], ['t','t']]
    means
    A = a c g - t
        a c g g t
    '''

    MA = []
    i = 0
    j = 0

    # print("extend M= " +str(M))
    while i < len(M) and j < len(A):
        # Case 1:
        if M[i][0] == '-' and A[j][0] == '-':
            M[i].append(A[j][1])
            MA.append(M[i])
            i = i + 1
            j = j + 1

        # Case 2:
        elif M[i][0] == '-' and A[j][0] != '-':
            M[i].append('-')
            MA.append(M[i])
            i = i + 1
        
        # Case 3:
        elif M[i][0] != '-' and A[j][0] == '-':
            c = ['-']*len(M[i])
            c.append(A[j][1])
            MA.append(c)
            j = j + 1
        
        # Case 4:
        elif M[i][0] != '-' and A[j][0] != '-':
            M[i].append(A[j][1])
            MA.append(M[i])
            i = i + 1
            j = j + 1

    if i < len(M):
        while i < len(M):
            M[i].append('-')
            MA.append(M[i])
            i = i + 1
            
    if j < len(A):
        k = len(M[0])
        while j < len(A):
            c = ['-']*(k-1)
            c.append(A[j][1])
            MA.append(c)
            j = j + 1
    return MA



def new_assembly_Gus_x(seqs, score_matrix, gap_cost, return_center_string=False, check_integrity=False):
    # STEP 1: Find the center string, s1
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # loop over all distinct pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    #making list to keep track of which string is which
    numbers = [str(i) for i in range(len(seqs))]
    # find the center string/guide 
    s1_idx = np.argmin(matrix.sum(axis=1))
    s1 = seqs[s1_idx]
    #print("umm, s1 is: ", str(s1))
    #print("and the index was: ", str(s1_idx) )
    seqs.insert(0, seqs.pop(s1_idx))  # move guide to the front of the list
    #keeping track of names
    y=[]
    y.append(numbers[s1_idx])
    numbers.remove(numbers[s1_idx])
    y.extend(numbers)
    names=y

    # STEP 2: Construct alignment M
    M: list[list[str]] = [[letter] for letter in [*s1]]
    cost_list = []
    for i in range(1, len(seqs)):
        cost = linear_C(gap_cost, score_matrix, s1, seqs[i])
        cost_list.append(get_cost_2(cost))

        # prepare A-matrix for extension
        alignment1_str, alignment2_str = linear_backtrack(s1, seqs[i], cost, score_matrix, gap_cost)
        alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
        A = [list(e) for e in zip(alignment1, alignment2)]

        # extend
        Mk = extend_alignment(M, A)
        M = Mk

    # ACTUALLY COMPUTE (approximate) COST
    total_cost = compute_cost(M, score_matrix, gap_cost)
    if check_integrity == True:
        integrity_check_Gus(seqs, M, score_matrix,gap_cost,matrix,s1_idx)
    return matrix, M, total_cost,names, s1_idx

def perform_updates_gradual(in_which_MSA_is_it, node1, node2,united_MSA_new, MSA_list):
    which_spot_in_MSA_list_to_update=min(in_which_MSA_is_it[node1][0],in_which_MSA_is_it[node2][0])
    which_spot_in_MSA_list_to_remove=max(in_which_MSA_is_it[node1][0],in_which_MSA_is_it[node2][0])
    MSA_list_copy=MSA_list
    MSA_list_copy[which_spot_in_MSA_list_to_update]=united_MSA_new
    MSA_list_copy.pop(which_spot_in_MSA_list_to_remove)
    companions_to_update=[]
    if in_which_MSA_is_it[node1][0]<in_which_MSA_is_it[node2][0]:
        #print(str(in_which_MSA_is_it[node1][0])+"is smaller than+"+str( in_which_MSA_is_it[node2][0]) +"so here I went to the original updating stategy")
        how_many_cols_already_in_MS1=find_highest_second_element(in_which_MSA_is_it,which_spot_in_MSA_list_to_update)
        for key, value in in_which_MSA_is_it.items():
            if value[0]==which_spot_in_MSA_list_to_remove:
                companions_to_update.append(key)
            if value[0]>which_spot_in_MSA_list_to_remove:
                value[0]=(value[0]-1)
        for companion in companions_to_update:
            #print(in_which_MSA_is_it[companion])
            in_which_MSA_is_it[companion][0]=which_spot_in_MSA_list_to_update
            col_of_element_in_old_MSA2=in_which_MSA_is_it[companion][1]
            in_which_MSA_is_it[companion][1]=(how_many_cols_already_in_MS1+col_of_element_in_old_MSA2+1)
    else:
        #print(str(in_which_MSA_is_it[node1][0])+"is bigger than+"+str( in_which_MSA_is_it[node2][0]) +"so here I went to the new updating stategy")
        companions_to_update1=[]
        companions_to_update2=[]
        how_many_cols_already_in_MS2=find_highest_second_element(in_which_MSA_is_it,which_spot_in_MSA_list_to_remove) #get the highest nr in MSA2. need to only update dict[0] for MSA1 and only dict[1] for MSA2
        for key, value in in_which_MSA_is_it.items():
            #how_many_cols_already_in_MS2=find_highest_second_element(in_which_MSA_is_it,which_spot_in_MSA_list_to_remove) #get the highest nr in MSA2. need to only update dict[0] for MSA1 and only dict[1] for MSA2
            if value[0]==which_spot_in_MSA_list_to_update: #these need their value[1] updated, because they were "pushed downwards" in their alignment
                companions_to_update1.append(key)
            if value[0]==which_spot_in_MSA_list_to_remove: #these need their value[0] updated, because they were moved to a new alignment.
                companions_to_update2.append(key)
            if value[0]>which_spot_in_MSA_list_to_remove: #these need -1 in their value[1]
                value[0]=(value[0]-1)
        for companion1 in companions_to_update1:
            in_which_MSA_is_it[companion1][1]=int(in_which_MSA_is_it[companion1][1])+int(how_many_cols_already_in_MS2)+1
        for companion2 in companions_to_update2:
            in_which_MSA_is_it[companion2][0]=which_spot_in_MSA_list_to_update
    return in_which_MSA_is_it, MSA_list_copy

def new_assembly_gradual_x(seqs,score_matrix,gap_cost, check_integrity=False):
    # Make a matrix to hold pairwise alignment costs for all alignment combinations!
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # Loop over all pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    matrix_for_MST=matrix #copy the matrix, so that we can keep the old matrix and make a changed version to the "pseudomatrix" version
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST) #making the "pseudomatrix"
    #print("matrix for MST: "+str(matrix_for_MST))
    min_span_edges_res=find_min_span_edges_testing(matrix_for_MST)
    in_which_MSA_is_it={}
    names= np.unique(min_span_edges_res[:, 2:])
    in_which_MSA_is_it ={name: [int(name),0] for name in names}
    MSA_list=[[[char] for char in seq] for seq in seqs]
    who_aligned_to_who=[]
    k=0
    for row in min_span_edges_res:
        k+=1
        if row[0]=="*":
            node1=row[2]
            node2=row[3]
            who_aligned_to_who.append([node1,node2])
            cost=linear_C(gap_cost,score_matrix,seqs[int(node1)],seqs[int(node2)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(node1)], seqs[int(node2)], cost, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            A = [list(e) for e in zip(alignment1,alignment2)]
            united_MSA_new=alt_alt_merge_united(A,MSA_list,in_which_MSA_is_it,node1,node2)
            in_which_MSA_is_it, MSA_list=perform_updates_gradual(in_which_MSA_is_it, node1, node2, united_MSA_new, MSA_list)
    if check_integrity==True:
        integrity_check_OBO_and_gradual(seqs, in_which_MSA_is_it, who_aligned_to_who, MSA_list,matrix,score_matrix,gap_cost)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    return(matrix,MSA_list, total_cost,in_which_MSA_is_it,who_aligned_to_who)

def find_min_span_edges_Prim(pseudomatrix, starting_node, verbose=False):
    r = str(starting_node)
    unprocessed_list = [[element] for element in pseudomatrix]

    if verbose:
        print(unprocessed_list)

    queue = []
    used = set()
    edges = []
    remove_from_unprocessed_list = []
    remove_from_queue_list = []
    run = True

    while run == True:
        for element in unprocessed_list:
            if (np.any(element[0][2] == str(r)) or np.any(element[0][3] == str(r))) and not any(
                    np.array_equal(element[0], item) for item in queue):
                if verbose:
                    print(element[0])
                queue.append(element[0])
                remove_from_unprocessed_list.append(element)

        if verbose:
            print("unprocessed_list at this time: \n ", str(unprocessed_list))
            print("this is what I wanna remove from unprocessed list:", str(remove_from_unprocessed_list))

        for y in remove_from_unprocessed_list:
            unprocessed_list = [element for element in unprocessed_list if not np.array_equal(element[0], y[0])]

        remove_from_unprocessed_list.clear()

        if verbose:
            print("queue after filling up for the iteration!: ", str(queue))

        used.add(int(r))
        u = 100000000
        edge = None

        for x in queue:
            if verbose:
                print("i just went into the queue")
                print(type(int(x[1])))
                print("x here is ", str(x), " and x[1] is", str(x[1]))

            if int(x[1]) < u:
                u = int(x[1])

                if int(x[2]) not in used:
                    r = int(x[2])
                    edge = x
                    if verbose:
                        print("I just for now assigned edge as", str(edge))
                elif int(x[3]) not in used:
                    r = int(x[3])
                    edge = x
                    if verbose:
                        print("I just for now assigned edge as", str(edge))
                else:
                    if verbose:
                        print("Yoikes, I found nothing to assign to edge because the if statements were not true.")

        used.add(r)
        edges.append(edge)

        if verbose:
            print("right now used are: ", str(used))

        for z in queue:
            if int(z[3]) in used and int(z[2]) in used:
                remove_from_queue_list.append(z)

        if verbose:
            print("queue before removing from it:", str(queue))
            print("I wanna remove: ", str(remove_from_queue_list))

        for a in remove_from_queue_list:
            queue = [element for element in queue if not np.array_equal(element, a)]

        remove_from_queue_list.clear()

        if verbose:
            print("queue after removing from queue after iteration: ", str(queue))

        if verbose:
            print("u is: ", str(u))
            print("r is: ", str(r))
            print("edges used right now are:", str(edges))

        if len(queue) < 1:
            if verbose:
                print("queue empty, you know...")
            run = False

    if verbose:
        print("end up with these edges", str(edges))
        
    return edges


def new_assembly_Prim_x(seqs,score_matrix,gap_cost, check_integrity=False):
    # Make a matrix to hold pairwise alignment costs for all alignment combinations!
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # Loop over all pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    matrix_for_MST=matrix #copy the matrix, so that we ca keep the old matrix and make a changed version to the "pseudomatrix" version
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST) #making the "pseudomatrix"
    min_span_edges_res=find_min_span_edges_Prim(matrix_for_MST, starting_node='0')
    in_which_MSA_is_it={}
    names=set()
    for element in min_span_edges_res:
        #print(element[2], element[3])
        names.add(element[2])
        names.add(element[3])
    in_which_MSA_is_it ={name: [int(name),0] for name in names}
    MSA_list=[[[char] for char in seq] for seq in seqs]
    who_aligned_to_who=[]
    for row in min_span_edges_res:
        node1=row[2]
        node2=row[3]
        who_aligned_to_who.append([node1,node2])
    for element in min_span_edges_res:
            node1=element[2]
            node2=element[3]
            cost=linear_C(gap_cost,score_matrix,seqs[int(node1)],seqs[int(node2)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(node1)], seqs[int(node2)], cost, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            A = [list(e) for e in zip(alignment1,alignment2)]
            united_MSA_new=alt_alt_merge_united(A,MSA_list,in_which_MSA_is_it,node1,node2)
            in_which_MSA_is_it, MSA_list=perform_updates_gradual(in_which_MSA_is_it, node1, node2, united_MSA_new, MSA_list)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    if check_integrity==True:
       integrity_check_OBO_and_gradual(seqs,in_which_MSA_is_it,who_aligned_to_who,MSA_list, matrix,score_matrix,gap_cost)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    return(matrix,MSA_list, total_cost,in_which_MSA_is_it,who_aligned_to_who)

def make_names_list(output_in_which_MSA:dict):
    names=[None]*len(output_in_which_MSA)
    for key, value in output_in_which_MSA.items():
        names[value[1]]=key
    return(names)

def make_str_from_cols(output_cols):
    #making each string into it's own list (instead of the lists being column)
    strings_list=[]
    for i in range(len(output_cols[0])):
        string=[]
        for col in output_cols:
            letter=col[i]
            string.append(letter)
        strings_list.append(string)
    #print(strings_list)
    #test=[char for string in strings_list for char in string]
    strings_list= [''.join(inner_list) for inner_list in strings_list]
    return strings_list

def make_related_strings(nr_of_str:int,len_of_str, degree_of_variation:float, start_string):
    string_fam_collection=[]
    first_string=[]
    if start_string==False:
        for k in range(0,int(len_of_str)):
                first_string.append(random.choice(['a','c','t','g']))
        first_string=''.join(first_string)
        string_fam_collection.append(first_string)
    else:
         first_string=start_string
         string_fam_collection.append(first_string)

    for i in range(1,nr_of_str):
        a_sequence=[]
        for element in first_string:
            if random.random() < degree_of_variation: #checking if we should change the charachter
                a_sequence.append(random.choice(['a','c','t','g']))
            else:
                 a_sequence.append(element)
        a_sequence=''.join(a_sequence)
        string_fam_collection.append(a_sequence)
    return string_fam_collection

def make_strings_in_families(nr_of_fams, nr_str_pr_fam, len_of_str, internal_var_in_fams, degree_of_var_from_first_fam):
    all_strings_coll = []
    first_family = make_related_strings(nr_str_pr_fam, len_of_str, internal_var_in_fams, start_string=False)
    all_strings_coll.append(first_family)
    
    for i in range(1, int(nr_of_fams)):
        carry_over_string = first_family[0]
        mutated_carry_over = []
        for element in carry_over_string:
            if random.random() < degree_of_var_from_first_fam:  # checking if we should change the character
                mutated_carry_over.append(random.choice(['a', 'c', 't', 'g']))
            else:
                mutated_carry_over.append(element)
        mutated_carry_over = ''.join(mutated_carry_over)  # Join the list of characters into a string
        family = make_related_strings(nr_str_pr_fam, len_of_str, internal_var_in_fams,mutated_carry_over)
        all_strings_coll.append(family)
    all_strings_coll = [item for sublist in all_strings_coll for item in sublist]
    #all_strings_coll = [''.join(sublist) for sublist in all_strings_coll]
    return all_strings_coll

import numpy as np
import heapq

def find_min_span_edges_Primzeyy(pseudomatrix, starting_node, verbose=False):
    r = str(starting_node)
    num_nodes = len(np.unique(pseudomatrix[:, 2:3]))
    edges = []
    unprocessed_set = [(int(weight), node1, node2) for _, weight, node1, node2 in pseudomatrix]
    used = set()
    used.add(r)

    while unprocessed_set:
        new_candidate_edges = [element for element in unprocessed_set if
                               (element[2] not in used and element[1] in used) or
                               (element[2] in used and element[1] not in used)]

        # Find the minimum weight edge
        element = min(new_candidate_edges)

        if element[1] not in used:
            used.add(element[1])
        if element[2] not in used:
            used.add(element[2])

        edges.append(element)

        # Update unprocessed_set
        unprocessed_set = [edge for edge in unprocessed_set if
                           edge[1] not in used or edge[2] not in used]

    edges = [np.array([''] + list(map(str, tpl)), dtype='<U21') for tpl in edges]
    return edges


def new_assembly_Primzeyy_x(seqs,score_matrix,gap_cost, check_integrity=False):
    # Make a matrix to hold pairwise alignment costs for all alignment combinations!
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # Loop over all pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    matrix_for_MST=matrix #copy the matrix, so that we ca keep the old matrix and make a changed version to the "pseudomatrix" version
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST) #making the "pseudomatrix"
    min_span_edges_res=find_min_span_edges_Primzeyy(matrix_for_MST, starting_node='0')
    print("I WENT PAST MAKING THE EDGES")
    print("just gonna print the new min span edges")
    print(min_span_edges_res)
    print(min_span_edges_res[0])
    print(min_span_edges_res[0][1])

    in_which_MSA_is_it={}
    names=set()
    print("ILL JUST PRINT MIN SPAN EGEDS")
    print(min_span_edges_res)
    for element in min_span_edges_res:
        print("gonna print an elemnet")
        print(element)
        print("gonna print element[0], 1 and 2")
        print(element[0])
        print("I said hey")
        print(element[1])
        print(element[2])
        print(element[3])
        names.add(element[2])
        names.add(element[3])
    in_which_MSA_is_it ={name: [int(name),0] for name in names}
    MSA_list=[[[char] for char in seq] for seq in seqs]
    who_aligned_to_who=[]
    for row in min_span_edges_res:
        node1=row[2]
        node2=row[3]
        who_aligned_to_who.append([node1,node2])
    for element in min_span_edges_res:
            print(element)
            node1=element[2]
            node2=element[3]
            cost=linear_C(gap_cost,score_matrix,seqs[int(node1)],seqs[int(node2)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(node1)], seqs[int(node2)], cost, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            A = [list(e) for e in zip(alignment1,alignment2)]
            united_MSA_new=alt_alt_merge_united(A,MSA_list,in_which_MSA_is_it,node1,node2)
            in_which_MSA_is_it, MSA_list=perform_updates_gradual(in_which_MSA_is_it, node1, node2, united_MSA_new, MSA_list)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    if check_integrity==True:
       integrity_check_OBO_and_gradual(seqs,in_which_MSA_is_it,who_aligned_to_who,MSA_list, matrix,score_matrix,gap_cost)
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    return(matrix,MSA_list, total_cost,in_which_MSA_is_it,who_aligned_to_who)