from Bio import Phylo
from io import StringIO

file = "w8-tree.new"
file3 = "three.new"
file4 = "four.new"

trees = Phylo.parse(file4, "newick")
for tree in trees:
    print(tree)

#print(tree.total_branch_length())
#print(tree.distance('Alpha', 'Delta'))

Phylo.draw(tree)
