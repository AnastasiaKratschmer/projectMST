import unittest
from Bio import Phylo
from io import StringIO
from rfdist import rfdist, get_all_leaves, get_splits, rfdist_from_files
from utils import parse_tree

def make_tree_from_string(str):
    return Phylo.read(StringIO(str), "newick")

class TestRFDistance(unittest.TestCase):
    # --- LEAVES
    def test_get_all_leaves_single_leaf(self):
        # Test when there's only one leaf in the tree
        tree = make_tree_from_string("(A);")
        leaves = get_all_leaves(tree)
        self.assertEqual(leaves, set(["A"]))

    def test_get_all_leaves_multiple_leaves(self):
        # Test when there are multiple leaves in the tree
        tree = make_tree_from_string("(A,(B,C),(D,E));")
        leaves = get_all_leaves(tree)
        self.assertEqual(leaves, set(["A", "B", "C", "D", "E"]))

    def test_get_all_leaves_duplicate_leaves(self):
        # Test when there are duplicate leaf names in the tree
        tree = make_tree_from_string("(A,(B,C),(D,E),(E,D));")
        leaves = get_all_leaves(tree)
        self.assertEqual(leaves, set(["A", "B", "C", "D", "E"]))

    # --- SPLITS
    def test_get_splits_single_internal_node(self):
        # A tree with one internal node should return one split
        # All splits are trivial??? so no splits
        tree = make_tree_from_string("(B,C,A);")
        # Phylo.draw_ascii(tree)
        splits = get_splits(tree)
        self.assertEqual(len(splits), 0)

    def test_get_splits_single_internal_node_2(self):
        # A tree with one internal node should return one split
        tree = make_tree_from_string("((A,B),(C,D));")
        splits = get_splits(tree)
        self.assertEqual(len(splits), 1) 

    def test_get_splits_multiple_internal_nodes(self):
        t = make_tree_from_string("(A,B)K,(C,(D,E)H)F;")
        splits = get_splits(t)
        self.assertEqual(len(splits), 2) 
    # --- DISTANCE 
    def test_identical_trees(self):
        # Two identical trees should have distance 0
        tree = make_tree_from_string("(A,B,(C,D));") 
        dist = rfdist(tree,tree)# rfdist_from_files(tree, tree)
        self.assertEqual(dist, 0)

    def test_single_different_leaf(self):
        # Two trees with only one leaf different should have distance 2
        t1 = make_tree_from_string("(A,B,(C,D));")
        t2 = make_tree_from_string("(A,D,(C,B));")
        dist = rfdist(t1, t2)
        self.assertEqual(dist, 2)

    def test_rfdist(self):
        t1 = "./testdata/tree1.new" # NOTE if using the VS code testing these files can't be found as of right now
        t2 = "./testdata/tree2.new"
        dist = rfdist_from_files(t1, t2)
        self.assertEqual(dist, 8)

if __name__ == '__main__':
    unittest.main()