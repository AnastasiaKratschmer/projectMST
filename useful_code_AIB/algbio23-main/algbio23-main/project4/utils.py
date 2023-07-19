
from Bio import Phylo
from io import StringIO

def print_newick_biopython(t):    
    """
    Print out the tree r in newick format using a depth-first traversal of the tree
    """
    df_print_newick_biopython(t.clade)
    print(";")

def df_print_newick_biopython(r):
    if r.is_terminal() == False:
        print("(", end="")
        for c in r.clades[:-1]:
            df_print_newick_biopython(c)
            print(",", end="")
        df_print_newick_biopython(r.clades[-1])
        print(")", end="")
    if r.name != None: print(r.name, end="")
    if r.branch_length != None: print(":", r.branch_length, sep="", end="")

'''
def parse_tree_from_file(tree_file):
    with open(tree_file, "r") as f:
        tree_string = f.read().strip()

    tree = Phylo.read(StringIO(tree_string), "newick")

    return tree
'''
def parse_tree(tree_input):
    return Phylo.read(tree_input, "newick")