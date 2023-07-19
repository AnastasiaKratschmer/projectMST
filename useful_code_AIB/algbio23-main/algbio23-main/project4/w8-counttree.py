#!/usr/bin/env python
#
# counttree.py <n>
#
# The number T(n,k) of unrooted trees with n leaves and k internal nodes
# of arbitrary degree is given by the recursion:
#
# T(n,k) = k*T(n-1,k) + (n+k-3)*T(n-1,k-1)
#
# Either a new leaf is added attached to one of the k internal nodes of a
# tree with k internal nodes and n-1 leaves (the term k*T(n-1,k)), or a new
# node of degree is three (with the new leaf attached) is "placed" on one of
# the (n-1)+(k-1)-1=n+k-3 edges of a tree with k-1 internal nodes and n-1
# leaves (the term (n+k-3)*T(n-1,k-1))
#
# where
#
# T(n,1) = 1 and T(n,k) = 0 if k > n-2
#
# The total number of tree with n leaves is
#
# T(n,1) + T(n,2) + ... + T(n,n-2)
#
# and the total number of binary trees is
#
# T(n,n-2)
#
# Christian Storm <cstorm@daimi.au.dk>, 6-July-2004

import os
import sys
import string

# Get n form the commandline
n = int(sys.argv[1])

# Allocate a 'row'
t = [0] * (n+1)

# Fill out the 'T-table' row by row
for i in range(3, n+1):
    total = prev = t[1] = 1
    for k in range(2,i-1):
        curr = t[k]
        t[k] = k*curr + (i+k-3)*prev
        prev = curr
        total = total + t[k]

    # Print the number of binary tree and arbitrary trees with i leaves
    try:
        frc = float(total)/float(t[i-2])
    except:
        frc = total/t[i-2]
    print("%3d : \tbin = %g \n\tall = %g \n\tfrc = %f\n" % (i, t[i-2], total, frc))
