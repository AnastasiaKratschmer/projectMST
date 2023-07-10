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
V=list(np.unique((sum(V,[]))))
print(V)

#A structure to keep the trees in, initialize one "tree" pr element to begin with!
element_trees = {}

for element in V:
    element_trees[element] = []

print(element_trees["a"])

#this needs to keep happening until everything is in one tree: 
while len(V)>1:
#choosing the min entry, and merge the two subtrees
    min0,min1,min2=E[0][0],E[0][1],E[0][2]
    element_trees[min1]=(","+min0+":"+min1+"->"+min2+str(element_trees[min2]))
    del(element_trees[min2]) #delete tree2 after merging it into tree1
    E=np.delete(E, 0, axis=0)
    print(V)
    V.remove(min2) #remove vertix from the list of vertices
    print(V)
    print(element_trees[min1])
    print(E)

#we get to the problem: what happens if the next elemenent has already had its tree deleted but it still has potential connections in E?
