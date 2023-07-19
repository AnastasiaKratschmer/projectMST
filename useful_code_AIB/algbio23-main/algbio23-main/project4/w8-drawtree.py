# drawtree.py
#
# Use Phylo from Biopython to read in a tree in newick-format.
# See the section 'Displaying trees' at https://biopython.org/wiki/Phylo
# for more details.
#
# Christian Storm Pedersen, 15-3-2020

import sys
from Bio import Phylo

tree = Phylo.read(sys.argv[1], "newick")
Phylo.draw_ascii(tree)
