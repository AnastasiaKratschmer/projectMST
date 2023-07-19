# Thoughts :)
## `sp_approx`
Could optimise finding the column costs using bit vectors (0 as gaps, 1 as letters) and shrink no. of lookups that way.

## `sp_exact_3`


## `sp_exact_3_firstver`
Oh no, we started by implemnting the below pseudocode but then realised it was too inefficient.
Pseudocode for recurrences fpr a nonboundary cell (i,j) 
(as written in Gusfield p. 344, section 14.6.1)
```
for i := 1 to n1 do
    for j := 1 to n2 do
        for k := 1 to n3 do
            begin
            if (S1(i) == S2(j)) then cij := smatch
            else cij := smis;
            if (S1(i) = S3(k)) then cik := smatch
            else cik := smis;
            if (S2(j) = S3(k)) then cjk := smatch
            else cjk := smis;

            d1 := D(i-1, j-1, k-1) + cik + cik + cjk;
            d2 := D(i-1, j-1, k) + cij + 2*sspace;
            d3 := D(i-1, j, k-1) + cik + 2*sspace;
            d4 := D(i, j-1, k-1) + cjk + 2*sspace;
            d5 := D(i-1, j, k) + 2*sspace;
            d6 := D(i, j-1, k) + 2*sspace;
            d7 := D(i, j, k-1) + 2*sspace;

            D(i, j, k) := Min ( d1, d2, d3, d4, d5, d6, d7);
            end;
```

### On the if-cases in Gusfield-pseudocode
We simply write, e.g. `cij = score_matrix[seq1[i-1], seq2[j-1]]` instead of the if-else statement comparing the letters and scoring accordingly based on it being a match or mismatch; looking up two letters in the score matrix gives us the score of a match (`smatch`) if they're the same letters and the score of that particular mismatch (`smis`) if they are not the same. 
The `sspace` is our gap cost. 

## `sp_approx`
We use the "pseudocode" from slide 13 of SP-MSA-Approx (week 4 lecture).

### On finding $S_1$ in step 1
The slides say to set $S_1$ and call the remaining strings $S_2, S_3, ...$, which in our program would correspond to reordering the list.
Instead, we just remember the index of $S_1$. 

## Scratchpad