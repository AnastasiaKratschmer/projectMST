from Bio import SeqIO
import os
import sys
import numpy as np
import random as random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


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
        print(record.name)
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

""" def get_cost_3(cost_matrix):
     return cost_matrix[-1,-1,-1] """

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



def fill_graph(res_matrix,source_node,layout="spring"): 
    edges=[]
    edges_in_min_path=[]
    for row in res_matrix:
        if row[0]=="*": #to make sure that only the edges in min span tree are added this is added for now. Could add the others to for illustratory purposes. But easier this way..
            weight = int(row[1])  # Extract the weight 
            node1 = row[2]  # Extract the first node
            node2 = row[3]  # Extract the second node 
            edges.append((node1, node2, {'weight': weight}))  # Add the tuple to the list
        
    for row in res_matrix:
        if row[0]=="*":  # Finding only edges determined to be included in the min path.
            node1 = row[2]  # Extract the first node 
            node2 = row[3]  # Extract the second node
            edge=tuple(sorted([node1,node2])) #needs be alphabethical
            edges_in_min_path.append(edge) # Add the tuple to the list 
    G = nx.Graph()
    G.add_nodes_from(np.unique(res_matrix[:, 2:]))
    G.add_edges_from(edges)
    shortest_path=edges_in_min_path
    #the following part is to choose the layout you want, but let's just keep the standard "spring"- My experience says it makes no diff.
    if layout=="spring":
        pos = nx.spring_layout(G)
    elif layout=="planar":
        pos=nx.planar_layout(G)
    elif layout=="shell":
        pos=nx.shell_layout(G)
    elif layout=="circular":
        pos=nx.circular_layout(G)
    elif layout=="spiral":
        pos=nx.spiral_layout(G)
    elif layout=="spectral":
        pos=nx.spectral_layout(G)
    else:
        print("Error! The specified layout is not available, choose spring, planar, shell, circular, spiral or spectral. Good luck!")
    edge_colors = ['deeppink' if e in shortest_path else 'lavender' for e in G.edges()] 
    nx.draw(G, pos, with_labels=True, node_color='bisque', edge_color=edge_colors, width=2, font_size=10)
    #nx.draw_networkx_edge_labels(G,pos, rotate=False)
    # Show the plot
    plt.ion() #to make it keep running even without manually closing fig!!
    plt.show()
    return G


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

def my_traversal_simply(graph, starting_key, verbose= False): #bf inspired traversal
    neighborhood = {}
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node)) #get the neighbors in the graph
        neighborhood[node] = neighbors
        if verbose: print(f"Neighbors of node {node}: {neighbors}") #a lil' unecessary print statement
    
    alignment_pairs = {}
    queue = [starting_key]  # Initialize the queue with the starting node, coming from the 'outside' of the function
    while queue:
        current_node = queue.pop(0)  # Get the first node from the queue, call it current_node while it is removed from the queue.
        
        if current_node in neighborhood:
            for successor in neighborhood[current_node]: #take the neighbors of the current node
                if successor not in alignment_pairs: # if its not already in the alignment pairs, the neighbor is a succesor, not a predecessor and should be handled
                    alignment_pairs[successor] = current_node #set the successor to be aligned to it's own predesessor, the "current node"
                    queue.append(successor)  # Add successor to the queue for further traversal handeling
    if starting_key in alignment_pairs: #the starting key should not be aligned TO anything (though other seq(s) should probably be aligned TO IT!)
        del alignment_pairs[starting_key]
        
    if verbose: print(alignment_pairs)
    #the second part of this func. makes a dict structure to keep track of the position of each string in the merging of all the pairwise alignments in extend_alignment_chaos :) 
    index_dict = {}
    index = 0

    # Add starting key with index 0
    index_dict[starting_key] = str(index) #the starting key will be the first string but into the MSA-merge
    index += 1

    # Add values from alignment_pairs with increasing indices, because the merging order follows the alignment_pairs-dictionary
    for key, value in alignment_pairs.items(): #the function crashes is value is not there. so it's there, although unaccessed :) 
        if key not in index_dict: #filling the index_dict with the positions of each node (representing the a string ) in the alignment_pairs
            index_dict[key] = str(index)
            index += 1
    return(alignment_pairs,index_dict)


def extend_alignment_chaos(M,str1_nr,A,index_dict, verbose: False): #needs inclusion of str1_nr, to come from the outside.... str1_nr is the name "predecessor" string  
    MA = []
    i = 0
    j = 0
    col_in_M_of_parent_string=int(index_dict[str1_nr]) #getting the column from the index dict that holds the "predecessor" string, to know where it is in the MSA
    while i < len(M) and j < len(A):
        if verbose: print("i:"+str(i)+", j:"+str(j))
        if verbose: print("parent string nr: "+ str(str1_nr))
        if verbose: print("which is in col:"+str(col_in_M_of_parent_string))
        if verbose: print("I compare "+ str(M[i][col_in_M_of_parent_string])+" and "+str(A[j][0]))
        # Case 1:
        if M[i][col_in_M_of_parent_string] == '-' and A[j][0] == '-':
            if verbose: print("I go to case 1 (two gaps)")
            M[i].append(A[j][1])
            MA.append(M[i])
            if verbose: print("Now MA is this: \n "+ str(MA))
            i = i + 1
            j = j + 1

        # Case 2
        elif M[i][col_in_M_of_parent_string] == '-' and A[j][0] != '-':
            if verbose:print("I go to case 2 (gap,char)")
            M[i].append('-')
            MA.append(M[i])
            if verbose: print("Now MA is this: \n "+ str(MA))
            i = i + 1
        
        # Case 3:
        elif M[i][col_in_M_of_parent_string] != '-' and A[j][0] == '-':
            if verbose:print("I go to case 3 (char, gap)")
            c = ['-']*len(M[i])
            c.append(A[j][1])
            MA.append(c)
            if verbose: print("Now MA is this: \n "+ str(MA))
            j = j + 1
        
        # Case 4:
        elif M[i][col_in_M_of_parent_string] != '-' and A[j][0] != '-':
            if verbose: 
                print("I go to case 4 (char, char)")
            M[i].append(A[j][1])
            MA.append(M[i])
            if verbose: print("Now MA is this: \n "+ str(MA))
            i = i + 1
            j = j + 1
    #when one of the strings has ended, put in gaps, till the other ends as well. 
    if i < len(M)-1:
        while i < len(M):
            M[i].append('-')
            MA.append(M[i])
            i = i + 1

     # Old verson that I'm not sure is correct, but now testing!
     # if j < len(A)-1:
     #   k = len(M[col_in_M_of_parent_string])
     #   while j < len(A)-1:
      #      c = ['-']*(k-1)         
    while j < len(A)-1:
        c = ['-']
        c.append(A[j][1])
        MA.append(c)
        j = j + 1
    return MA

def new_sp_approxi_combi(seqs: list[str], score_matrix: dict, gap_cost: int, verbose=False, return_center_string=False,layout="spring"):
    # Make a matrix to hold pairwise alignment costs for all alignment combinations!
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # Loop over all pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2,verbose=verbose))
    if verbose:
        print("Here comes the distance matrix produced by the alignments: \n")
        print(matrix)
    matrix_for_MST=matrix #copy the matrix, so that we can keep the old matrix and make a changed version to the "pseudomatrix" version
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST) #making the "pseudomatrix"
    min_span_edges=find_min_span_edges_testing(matrix_for_MST,verbose=verbose) #Run Kruskal's algorithm on the "pseudomatrix"
    if verbose:
        print("Here comes the pseudomatrix, filled out with with the edges inclued in the MST: \n")
        print(min_span_edges)

    max_indices = np.where(matrix == np.max(matrix)) # Choosing where to start traversal. I want to start at one of the nodes that is the furthest away from any other so start from a side of graph.. hmmmm...
    max_row_index = max_indices[0][0] #just choose one of them.
    if verbose:
        print("Starting key for traversal based on max_row_idex: ")
        print(max_row_index)

    #Put the nodes and the minimum spanning edges into a graph.
    G=fill_graph(min_span_edges,str(int(max_row_index)),layout) #using the max_row_index as the starting key! (an making the graph!)
    alignment_pairs,index_dict=my_traversal_simply(G,str(int(max_row_index)),verbose=verbose) #'traverse' to get alignment_pairs (pairs of sucessors and predecessors) and their position in the MSA to come (index_dict)
    if verbose:
        print("Here come your alignment pairs and the idex dict: \n") 
        print(alignment_pairs)
        print(index_dict)

     # Constructing alignment M
    M: list[list[str]] = [[letter] for letter in [*seqs[int(max_row_index)]]] #make structure where evey column in the alignment is represented as a string in a list (in a list)
    cost_list = []
    A_dict={}
    #using the pairings of predecessors and successors in the alignment_pairs dict, align the strings.
    for key,value in alignment_pairs.items():
        if verbose: print("this is the key: "+str(key)+" and this is the value: "+str(value))
        cost = linear_C(gap_cost, score_matrix, seqs[int(value)], seqs[int(key)],verbose=verbose) #the alignment call itself :) 
        if verbose: print("\n now aligning...."+str(seqs[int(key)])+ " and "+ str(seqs[int(value)]))
        cost_list.append(get_cost_2(cost))
        
        # prepare A-matrix for extension
        alignment1_str, alignment2_str = linear_backtrack(seqs[int(value)], seqs[int(key)], cost, score_matrix, gap_cost,verbose=verbose) #backtract to get the alignments!
        str1_nr=value #the predecessor/parent string
        alignment1, alignment2 = [*alignment1_str], [*alignment2_str] #splitting up the alignments into elements to have the right format for the list of lists (M)
        
        A = [list(e) for e in zip(alignment1,alignment2)] #zipping the elements of the two aligned strings together pairwisely
        if int(key) < int(value):
            key_for_A_dict=str(key)+"_"+str(value)
        else:
            key_for_A_dict=str(value)+"_"+str(key)
        A_dict[key_for_A_dict]=A
        #pair_al=A_dict[key_for_A_dict]
        if verbose: print("A right now is: "+str(A))
        if verbose: print("M right now: "+str(M))
        # extend
        Mk = extend_alignment_chaos(M,str1_nr, A,index_dict,verbose=verbose) 
        M = Mk
    if verbose:
        print("Here is the alignment in full omg: \n")
        print(M)
    print("this is A_dict: ")
    print(A_dict)
    al_integrity_testt(A,alignment_pairs,index_dict,M,verbose=verbose)
    # ACTUALLY COMPUTE (approximate) COST
    total_cost = compute_cost(M, score_matrix, gap_cost)
    print("Total cost of MSA:"+str(total_cost))
    return total_cost, M, matrix_for_MST,G


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

def al_integrity_testt(A_dict, alignment_pairs, index_dict, MSA,verbose=False):
    integrity=True
    for key, value in alignment_pairs.items():
        position_key = index_dict[key] #find the position in the MSA
        position_value = index_dict[value]
        seq1 = []
        seq2 = []
        for element in MSA:
            seq1.append(element[int(position_key)]) #extract the 'column' for the
            seq2.append(element[int(position_value)])
        zipped = [list(e) for e in zip(seq1, seq2)] #zip the two lists to get the same format as in the pairwise alignment!
        zipped = [sublist for sublist in zipped if not all(element == '-' for element in sublist)]#removing all-gap columns

        if int(key) < int(value):
            key_for_A_dict = str(key) + "_" + str(value) #to reflect naming practice in the big function..
        else:
            key_for_A_dict = str(value) + "_" + str(key)

        # Add debugging print statement
        if verbose:
            print("key_for_A_dict:", key_for_A_dict)

        # Check if key_for_A_dict exists in A_dict
        if key_for_A_dict in A_dict:
            pair_al = A_dict[key_for_A_dict]
            for i in range(min(len(zipped), len(pair_al))): #iterate and compare columns!ch
                if zipped[i] != pair_al[i]:
                    print("A problem occurred: At position " + str(i) + ", " + str(zipped[i]) + " is not " + str(pair_al[i]))
                    integrity=False
            print("A problem occurred: key_for_A_dict not found in A_dict")
        if integrity and verbose:
            print("Alignment integrity looks fine!")