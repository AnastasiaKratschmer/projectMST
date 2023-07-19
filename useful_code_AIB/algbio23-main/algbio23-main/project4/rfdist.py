import numpy as np
import sys
from Bio import Phylo
from io import StringIO
from tqdm import tqdm  # loading bar
from utils import print_newick_biopython, parse_tree

# note: i wrote these leaf functions but i guess we don't need them anyways but i'm just gonna leave them and their tests for now just in case
def get_all_leaves(tree):
    # extracts all the leaves(' names) from a Phylo tree object
    # note:  right now this function does not do anything special if there are different label names
    return set([clade.name for clade in tree.get_terminals()])

def count_shared_leaves(leaves1, leaves2):
    return len(leaves1.intersection(leaves2))

def get_splits(tree):
    all_leaves = get_all_leaves(tree)
    splits = set()
    splits_opposite = set()
    # add split for each internal node
    internal_nodes = tree.get_nonterminals()
    # we remove the first node, this is some kind of root node that phylo adds
    if len(internal_nodes)>1:
        internal_nodes.pop(0)
    for node in internal_nodes:
        # we remember opposite splits as well to check that we don't add the reverse splits :)
        split = frozenset([leaf.name for leaf in node.get_terminals()])
        split_opposite = frozenset(all_leaves - split)
        length_is_correct = (len(split) < (len(all_leaves))) and (len(split) > 1)
        # ensure that it does not add 'reverse' split
        if split not in splits and split not in splits_opposite and length_is_correct: # only add unique splits and only non-trivial splits :)
            splits.add(split)
            splits_opposite.add(split_opposite)
    return splits

def count_shared_splits(t1_splits, t2_splits):
    return len(t1_splits.intersection(t2_splits))

def rfdist_from_leaves_and_splits( t1_splits, t2_splits, num_shared_splits):
    return len(t1_splits) + len(t2_splits) - 2*num_shared_splits  

# compute RF distance from two Bio.Phylo trees
def rfdist(t1, t2):    
    # get splits and count shared splits
    t1_splits = get_splits(t1)
    t2_splits = get_splits(t2) 
    num_shared_splits = count_shared_splits(t1_splits, t2_splits)
    
    # actually compute the distance
    return rfdist_from_leaves_and_splits(t1_splits, t2_splits, num_shared_splits)


def rfdist_from_files(t1_file, t2_file, Verbose=False):
    t1 = parse_tree(t1_file)
    t2 = parse_tree(t2_file)

    if Verbose:
        print("Tree 1:")
        Phylo.draw_ascii(t1)
        print_newick_biopython(t1)
        print("Tree 2:")
        Phylo.draw_ascii(t2)
        print_newick_biopython(t2) 

    return rfdist(t1, t2)

# ------------------------------ ACTUAL PROGRAM ------------------------------
if __name__ == "__main__":
    # Read arguments ...
    if not len(sys.argv) == 3:
        sys.stderr.write("USAGE: python3 %s < two evolutionary trees in Newick format > "% sys.argv[0])
        sys.exit(1)
    t1_input, t2_input = sys.argv[1:]
    
    # ----- What we see in the terminal
    print("Beep boop, dubbedy doop!\n")
    print("Computing RF distance...")
    dist = rfdist_from_files(t1_input, t2_input) # NOTE add Verbose=True if you want the descriptive prints
    print("Done!\n ")
    print("The distance is " + str(dist) +".")
    print()