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
