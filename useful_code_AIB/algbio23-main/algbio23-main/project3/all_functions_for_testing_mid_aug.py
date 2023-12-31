from Bio import SeqIO
import os
import sys
import numpy as np



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
                # print("Scores for this row/letter: ", scores_for_this_row)
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

    if verbose: print("Cost matrix:\n" + str(matrix))
    return(matrix)

def get_cost_2(cost_matrix):
    return cost_matrix[-1,-1]

def get_cost_3(cost_matrix):
    return cost_matrix[-1,-1,-1]

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
        if M[i][-1] == '-' and A[j][0] == '-':
            M[i].append(A[j][1])
            MA.append(M[i])
            i = i + 1
            j = j + 1

        # Case 2:
        elif M[i][-1] == '-' and A[j][0] != '-':
            M[i].append('-')
            MA.append(M[i])
            i = i + 1
        
        # Case 3:
        elif M[i][-1] != '-' and A[j][0] == '-':
            c = ['-']*len(M[i])
            c.append(A[j][1])
            MA.append(c)
            j = j + 1
        
        # Case 4:
        elif M[i][-1] != '-' and A[j][0] != '-':
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
        k = len(M[-1])
        while j < len(A):
            c = ['-']*(k-1)
            c.append(A[j][1])
            MA.append(c)
            j = j + 1
    return MA

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class Tree:
    def __init__ (self,name):
        self.name=name
        self.next=[]
        self.components=[self.name]
    def __str__(self):
        return f"{self.name} {self.next}"
    def __repr__(self):
        return f"Tree(name='{self.name}', next={self.next})"
    @classmethod
    def add(cls, tree1, tree2):
        if not tree1.next:
            tree1.next.append(tree2)
        else:
            tree1.next.append(tree2)
        tree1.components = tree1.components + tree2.components
        del tree2

def find_min_span_edges(pseudomatrix):
    sorted_indices = np.lexsort((pseudomatrix[:, 1].astype(int),))
    E = pseudomatrix[sorted_indices]
    print("this is E (sorted matrix without any stars yet): ")
    print(E)
    names= np.unique(pseudomatrix[:, 2:])
    name_dict= {letter:i for i, letter in enumerate(names)}
    E[:,0]=''
    trees=[]
    for item in names:
        trees.append(Tree(item))
        x=0
    for i,item in enumerate(E):#is this enumeration even needed????????
        while len(set(name_dict.values()))>1:
            min0,min1,min2=E[x][1],E[x][2],E[x][3] #get the next line in the matrix to get the letters/names of the trees
            tree1_id= name_dict[min1] #getting the positions of the letters in the current list of trees
            tree2_id= name_dict[min2]
            if tree2_id==tree1_id:
                x+=1
                break
            else:
                E[x][0]="*"
                tree1=trees[tree1_id] #get the first tree
                tree2=trees[tree2_id] #get he second tree
                Tree.add(tree1,tree2) #merge the two trees, meaning that our tree-list gets shorter. PROBLEM THOUGH: the last it. does not fully merge right now----
                trees.pop(tree2_id) #remove the tree that we just merged into another tree
                if tree1_id<tree2_id:
                    orig_tree2_id=tree2_id
                    name_dict[min2]=tree1_id # update the number associated with the letter/treee in the dictionary
                    for key in name_dict:
                        if name_dict[key] >= orig_tree2_id:
                            name_dict[key] -= 1
                else:
                    orig_tree1_id=tree1_id
                    name_dict[min1]=tree2_id # update the number associated with the letter/treee in the dictionary
                    for key in name_dict:
                        if name_dict[key] >= orig_tree1_id:
                            name_dict[key] -= 1
                x+=1
    res_mat=E
    return res_mat

def get_visiting_order(res_matrix,source_node,traversal="df",layout="spring"):
    edges=[]
    edges_in_min_path=[]
    for row in res_matrix:
        if row[0]=="*": #was not there before, but to make sure only the edges in min span tree are added this is added for now.
            weight = int(row[1])  # Extract the weight 
            node1 = row[2]  # Extract the first node
            node2 = row[3]  # Extract the second node 
            edges.append((node1, node2, {'weight': weight}))  # Add the tuple to the list
        
    for row in res_matrix:
        if row[0]=="*":  # Finding only edges determined to be included in the min path.
            node1 = row[2]  # Extract the first node 
            node2 = row[3]  # Extract the second node
            edge=tuple(sorted([node1,node2])) #needs be alphabethical
            edges_in_min_path.append(edge) # Add the tuple to the list #WORKS TILL HERE
    G = nx.Graph()#WORKS TILL HERE 
    G.add_nodes_from(np.unique(res_matrix[:, 2:])) #WORKS TILL HERE
    G.add_edges_from(edges) #WORKS TILL HERE 
    shortest_path=edges_in_min_path #WORKS TILL HERE
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
    #pos = nx.spring_layout(G) #WORKS TILL HERE
    edge_colors = ['deeppink' if e in shortest_path else 'lavender' for e in G.edges()] #WORKS TILL HERE
    nx.draw(G, pos, with_labels=True, node_color='bisque', edge_color=edge_colors, width=2, font_size=10) #WORKS TILL HERE
    # Show the plot
    plt.ion() #to make it keep running even without manually closing fig!!
    plt.show() #SO I GUESS: WORKS TILL HERE
    if traversal=="bf":
        order = list(nx.bfs_tree(G, source=source_node))
        print("my traversal is bf") #I guess we can use the middle string as the source node or something. or one of the ones that's the most different from the others to hopefully start "at a side".
    #and really consider using another method than depth-first! like "nx.bfs_edges" #WORKS TILL HERE
    if traversal=="df":
        order = list(nx.dfs_tree(G, source=source_node))
        print("My traversal is df")
    return order,G

def convert_format_mat_to_pseudomat(mini_mat):
    # Get the node names (excluding the first row and first column)
    node_names = mini_mat[0, 1:]
    
    # Initialize an empty list to store the rows of the new format
    matrixx_rows = []
    print(len(mini_mat))
    print(mini_mat)
    x=len(mini_mat)-1
    processed_edges=set()
    # Process the mini_mat to generate the rows of the new format
    for i, row in enumerate(mini_mat[1:, 1:]):
        for j, distance in enumerate(row):
            if i != j:
                node1 = mini_mat[i + 1, 0]
                node2 = mini_mat[0, j + 1]
                edge = (node1, node2) if node1 < node2 else (node2, node1)
                if edge not in processed_edges:
                    matrixx_rows.append(["", int(distance), node1, node2])
                    processed_edges.add(edge)
    
    # Convert the list of rows to a NumPy array
    matrixx_np = np.array(matrixx_rows)
    
    return matrixx_np

def convert_to_desired_format2(distance_matrix):
    # Get the number of nodes in the distance matrix
    num_nodes = len(distance_matrix)

    # Initialize an empty list to store the rows of the new format
    matrixx_rows = []

    # Process the distance_matrix to generate the rows of the new format
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):  # Avoid duplicates and self-distances
            distance = int(distance_matrix[i, j])
            node1 = chr(ord('A') + i)  # Node names as A, B, C, ...
            node2 = chr(ord('A') + j)
            matrixx_rows.append(["", distance, node1, node2])

    # Convert the list of rows to a NumPy array
    matrixx_np = np.array(matrixx_rows)

    return matrixx_np

def convert_to_desired_format_nr_version(distance_matrix):
    # Get the number of nodes in the distance matrix
    num_nodes = len(distance_matrix)

    # Initialize an empty list to store the rows of the new format
    matrixx_rows = []

    # Process the distance_matrix to generate the rows of the new format
    for i in range(num_nodes):
        for j in range(i+1, num_nodes):  # Avoid duplicates and self-distances
            distance = int(distance_matrix[i, j])
            node1 = str(i)  # Node names as A, B, C, ...
            node2 = str(j)
            matrixx_rows.append(["", distance, node1, node2])

    # Convert the list of rows to a NumPy array
    matrixx_np = np.array(matrixx_rows)
    return matrixx_np

def nothing(x):
    return x*33


def find_min_span_edges_testing(pseudomatrix):
    sorted_indices = np.lexsort((pseudomatrix[:, 1].astype(int),))
    E = pseudomatrix[sorted_indices]
    print("this is E (sorted matrix without any stars yet): ")
    print(E)
    names= np.unique(pseudomatrix[:, 2:])
    name_dict= {letter:i for i, letter in enumerate(names)}
    E[:,0]=''
    trees=[]
    for item in names:
        trees.append(Tree(item))
        x=0
    it=0
    for i,item in enumerate(E):#is this enumeration even needed????????
        while len(set(name_dict.values()))>1:
            it+=1
            print("\n \n this is it "+str(it))
            min0,min1,min2=E[x][1],E[x][2],E[x][3] #get the next line in the matrix to get the letters/names of the trees
            tree1_id= name_dict[min1] #getting the positions of the letters in the current list of trees
            tree2_id= name_dict[min2]
            print("before anything happens: tree1_id is "+str(tree1_id)+" and tree2_id is "+str(tree2_id))
            print(str(min0)+" "+str(min1)+" "+str(min2))
            if tree2_id == tree1_id:
                x += 1
                break
            else:
                E[x][0] = "*"
                tree1 = trees[tree1_id]
                tree2 = trees[tree2_id]
                Tree.add(tree1, tree2)
                trees.pop(tree2_id)
                print("in namedict, before merging we have the numbers: " + str(name_dict[min1]) + " and " + str(name_dict[min2]))
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
            print("in namedict, now we have the numbers: "+str(name_dict[min1])+" and "+str(name_dict[min2]))
            print("and the whole dict is:"+str(name_dict))
            
    res_mat=E
    return res_mat

def my_traversal_simply(graph, starting_key):
    neighborhood = {}
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        neighborhood[node] = neighbors
        print(f"Neighbors of node {node}: {neighbors}")
    
    alignment_pairs = {}
    queue = [starting_key]  # Initialize the queue with the starting node
    while queue:
        current_node = queue.pop(0)  # Get the first node from the queue
        
        if current_node in neighborhood:
            for successor in neighborhood[current_node]:
                if successor not in alignment_pairs:
                    alignment_pairs[successor] = current_node
                    queue.append(successor)  # Add successor to the queue for further traversal
    if starting_key in alignment_pairs:
        del alignment_pairs[starting_key]
        
    print(alignment_pairs)
    pos_dict={}
    #pos_dict_addition={}
    #pos_dict[starting_key]=0
    #pos_dict_addition={key: str(index+1) for index, key in enumerate(alignment_pairs.values())}
    #print("pos dict: "+str(pos_dict))
    #print("pos dict add: "+str(pos_dict_addition))
    #merged_dict = {**pos_dict, **pos_dict_addition}

    #return(alignment_pairs, merged_dict)
    #pos_dict = {}
    #pos_dict[starting_key] = 0
    #pos_dict.update({key: str(index) for index, key in enumerate(alignment_pairs.values())})
    #return(alignment_pairs,pos_dict)
    index_dict = {}
    index = 0

    # Add starting key with index 0
    starting_key = '0'
    index_dict[starting_key] = str(index)
    index += 1

    # Add values from alignment_pairs with increasing indices
    for key, value in alignment_pairs.items():
        if key not in index_dict:
            index_dict[key] = str(index)
            index += 1
        if value not in index_dict:
            index_dict[value] = str(index)
            index += 1
    return(alignment_pairs,index_dict)

#al,pos_dict=my_traversal_simply(graph,"0")
#print(al,pos_dict)

def extend_alignment_chaos(M,str1_nr,A,position_dictionary):#needs inclusion of str1_nr to come from the outside....
    MA = []
    i = 0
    j = 0
    col_in_M_of_parent_string=int(position_dictionary[str1_nr])
    # print("extend M= " +str(M))
    while i < len(M) and j < len(A):
        # Case 1:
        if M[i][col_in_M_of_parent_string] == '-' and A[j][0] == '-':
            M[i].append(A[j][1])
            MA.append(M[i])
            i = i + 1
            j = j + 1

        # Case 2:
        elif M[i][col_in_M_of_parent_string] == '-' and A[j][0] != '-':
            M[i].append('-')
            MA.append(M[i])
            i = i + 1
        
        # Case 3:
        elif M[i][col_in_M_of_parent_string] != '-' and A[j][0] == '-':
            c = ['-']*len(M[i])
            c.append(A[j][1])
            MA.append(c)
            j = j + 1
        
        # Case 4:
        elif M[i][col_in_M_of_parent_string] != '-' and A[j][0] != '-':
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
        k = len(M[col_in_M_of_parent_string])
        while j < len(A):
            c = ['-']*(k-1)
            c.append(A[j][1])
            MA.append(c)
            j = j + 1
    return MA

def new_sp_approxi_combi(seqs: list[str], score_matrix: dict, gap_cost: int, verbose=False, return_center_string=False,layout="spring",traversal_here="df"):
    # STEP 1: Find the center string, s1
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # loop over all distinct pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
              matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
    print(matrix)
    matrix_for_MST=matrix 
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST)
    min_span_edges=find_min_span_edges_testing(matrix_for_MST)
    print(matrix_for_MST)
    print(min_span_edges)

    max_indices = np.where(matrix == np.max(matrix))
    max_row_index = max_indices[0][0]
    print(max_row_index)

    #need to use the my_transversal_simply but already to put in the graph!
    visiting_order,G=get_visiting_order(min_span_edges,str(int(max_row_index)),traversal_here,layout) #str(int(np.argmin(matrix.sum(axis = 1)))) to get the min arg as the starting point.
    #print(visiting_order) #visiting order is now letters, but we would need that as numbers/idices from the score matrix to keep track. #APPARENTLY NOT ACTUALLY IN USE
    alignment_pairs,index_dict=my_traversal_simply(G,str(int(max_row_index)))
    print(alignment_pairs)
    print(index_dict)

    # find center string/guide 
    s1_idx = np.argmin(matrix.sum(axis = 1))
    print("nr of seqs:")
    print(len(seqs))
    #seqs_new = [None for _ in range(len(seqs))]
    #s1 = seqs[s1_idx]
    #seqs_new[0]=seqs[s1_idx]
    #for i,element in enumerate(visiting_order,start=1):
     #       seqs_new[i]=seqs[int(element)]
    seqs_new = [None for _ in range(len(seqs))]
    s1 = seqs[s1_idx] #should be "seqs[s1_idx]", but is set to 2 for testing purposes
    seqs_new[0] = seqs[int(max_row_index)] #should be "seqs[s1_idx]", but is set to 2 for testing purposes
    i=1
    print("new seqs before loop")
    print(seqs_new)
    new_s1_index=None
    for i, idx in enumerate(visiting_order):
        seqs_new[int(idx)] = seqs[i]
        if i==s1_idx:
            new_s1_index=idx
        print("idx:"+str(idx)+", i="+str(i))
    print("new s1-index:")
    print(new_s1_index)
    print("seqs new:")
    print(seqs_new)
    print("old seqs:")
    print(seqs)
    print("matrix for source node inspo:")
    print(matrix)
         
         
    #seqs.insert(0, seqs.pop(s1_idx)) # move guide to front of list
    if verbose: print("The center string, s1, is sequence no." + str(s1_idx+1)) # just a print statement to see which string is the center string

    # STEP 2: Construct alignment M
    M: list[list[str]] = [[letter] for letter in [*seqs[int(max_row_index)]]]
    print("M right now:")
    print(M)
    print("s1 right now:")
    print(s1)
    print("seqs right now")
    print(seqs)
    cost_list = []
    # print("first M = \n" + str(M))
    for key,value in alignment_pairs.items(): #was just seqs bebore
        #you need here to get the pairs from the new dict instead
        cost = linear_C(gap_cost, score_matrix, seqs[int(value)], seqs[int(key)]) #was just seqs before
        print("\n now aligning...."+str(seqs[int(value)])+ " and "+ str(seqs[int(key)]))
        cost_list.append(get_cost_2(cost))
        
        # prepare A-matrix for extension
        alignment1_str, alignment2_str = linear_backtrack(seqs[int(value)], seqs[int(key)], cost, score_matrix, gap_cost)
        str1_nr=value
        alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
        #here you need to provide the nr of the sequences to be able to index in the indexing dict in the extend_alignment_chaos func
        A = [list(e) for e in zip(alignment1,alignment2)]
        print("A right now is: "+str(A))
        print("M right now: "+str(M))
        # extend
        Mk = extend_alignment_chaos(M,str1_nr, A,index_dict) #needs the chaos version
        M = Mk
    
    # ACTUALLY COMPUTE (approximate) COST
    total_cost = compute_cost(M, score_matrix, gap_cost)

    #if return_center_string: return total_cost, M, s1_idx #IDK why this runs for no reason when it's uncommented..
    return total_cost, M, matrix_for_MST, visiting_order,G

def my_traversal_simply(graph, starting_key):
    neighborhood = {}
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        neighborhood[node] = neighbors
        print(f"Neighbors of node {node}: {neighbors}")
    
    alignment_pairs = {}
    queue = [starting_key]  # Initialize the queue with the starting node
    while queue:
        current_node = queue.pop(0)  # Get the first node from the queue
        
        if current_node in neighborhood:
            for successor in neighborhood[current_node]:
                if successor not in alignment_pairs:
                    alignment_pairs[successor] = current_node
                    queue.append(successor)  # Add successor to the queue for further traversal
    if starting_key in alignment_pairs:
        del alignment_pairs[starting_key]
        
    print(alignment_pairs)
    pos_dict={}
    #pos_dict_addition={}
    #pos_dict[starting_key]=0
    #pos_dict_addition={key: str(index+1) for index, key in enumerate(alignment_pairs.values())}
    #print("pos dict: "+str(pos_dict))
    #print("pos dict add: "+str(pos_dict_addition))
    #merged_dict = {**pos_dict, **pos_dict_addition}

    #return(alignment_pairs, merged_dict)
    #pos_dict = {}
    #pos_dict[starting_key] = 0
    #pos_dict.update({key: str(index) for index, key in enumerate(alignment_pairs.values())})
    #return(alignment_pairs,pos_dict)
    index_dict = {}
    index = 0

    # Add starting key with index 0
    #starting_key = '0' #I guess this should not be???
    index_dict[starting_key] = str(index)
    index += 1

    # Add values from alignment_pairs with increasing indices
    for key, value in alignment_pairs.items():
        if key not in index_dict:
            index_dict[key] = str(index)
            index += 1
        if value not in index_dict:
            index_dict[value] = str(index)
            index += 1
    return(alignment_pairs,index_dict)

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


print(convert_to_desired_format_nr_version(np.array([[0, 1, 2],[1, 0, 3],[2, 3, 0]])))