#I wanna start understanding Kruskal's algo. First, I need to find a way to input edges. Let's start by making it as a matrix.
import numpy as np
E_matrix=np.array([[int(19),"a","b"],[int(2),"b","c"],[int(2),"d","c"],[int(4),"d","b"],[int(7),"d","a"]])
print(E_matrix)
#ok, for now I guess we can use that structure
print(E_matrix[1]) #so this is how to get a single line

#sorting the matrix based  on "edge weights"
sorted_indices = np.lexsort((E_matrix[:, 0].astype(int),))
E = E_matrix[sorted_indices]
print(E)

#Get vertices from edge input
V=[[item[1],item[2]] for item in E_matrix]
V=np.unique((sum(V,[])))
print(V)

#A structure to keep the trees in, initialize one "tree" pr element to begin with!
element_trees = {}

for element in V:
    element_trees[element] = []

print(element_trees["a"])