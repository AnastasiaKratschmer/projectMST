# Algorithms in Bioinformatics 2023

**Group**: Group 6

**Group members**: 
- Anastasia Kratschmer (201907407)
- Noa Van de Velde (202209940) 
- Stinna Danger (201907211)
- Anna Kathrine Skov (201905355)

## Project 3 - How to run
*Note that you should be in the `project3/` directory.* 

### `sp_exact_3` and `sp_approx`
The programs take the following arguments:

* For `sp_exact_3` provide sequences like so:
    * Sequence one, provided as a string or as a path to a FASTA file.
    * Sequence two, provided as a string or as a path to a FASTA file.
    * Sequence three, provided as a string or as a path to a FASTA file.
* For `sp_approx` provide sequences like so:
    * Sequences provided as a path to a FASTA file.
* A gap cost.
* A scoring matrix, provided as a path to a file containing the matrix writtin in Phylip format (note that our program uses `tab`s as separators between the numbers).

Example scoring matrix: 
```
4
A	0	5	2	5
C	5	0	5	2
G	2	5	0	5
T	5	2	5	0
```
You can run the exact algorithm with the following line in the terminal:
```
python .\sp_exact_3.py GTTCCGAAAGGCTAGCGCTAGGCGCC ATGGATTTATCTGCTCTTCG TGCATGCTGAAACTTCTCAACCA 5 .\testdata\score_matrix.txt
```
The program will compute the cost of the alignment and backtrack and output the alignment itself. The above run would look like so:
```
Beep boop!

Computing cost of alignment...
Backtracking...
Done!

Cost:
 198.0
Alignments:
 gttccgaaaggctagcgctaggc-gcc-
 a-t--g-gat-tt-at-ctgctc-ttcg
 --t--g-catgctaaaccttctcaacca
```

You can run the approximation algorithm like so:

```
python .\sp_approx.py .\testdata\brca1-testseqs.fasta 5 .\testdata\score_matrix.txt  
```

```
python .\sp_approx.py .\testdata\testdata_short_approx.fasta 5 .\testdata\score_matrix.txt
```
The program computes and outputs the approximate cost of aligning the given sequences. 

## Experiments and tests
Our tests can be run in with the `tests.py` script.
```
python .\tests.py
```
The experiments for the report can be run with the `experiments.py` script, which prompts the user on which experiment they would like to run, repeating this until the user exits. Figures are output to the `figures/` folder.
```
> python .\experiments.py
Beep boop!
Which experiment would you like to run, 1, 2, or 3?
> 3
--- Third experiment ---
100%|███████████████████████████████████████████████████████████████| 20/20 [03:52<00:00, 11.64s/it]
On average, the approximation ratio was: 
1.1789242853758934
Generated are plots in the figures-folder. :)

Experiment complete!
Which experiment would you like to run now, 1, 2, or 3? Press X to exit.
1
> 2
--- First experiment ---
The cost of alignment of the three sequences, brca1_bos_taurus, brca1_canis_lupus and brca1_gallus_gallus:
790.0

Experiment complete!
Which experiment would you like to run now, 1, 2, or 3? Press X to exit.
> x
Okay, bye!
```



The alignment of the brca1 sequences as needed by the second part of the project is run with the `brca1-alignment.py` script and the resulting alignment is output to the `brca1-alignment.fasta`.
```
> python .\brca1-alignment.py   
Beep boop!
Aligning the sequences corresponds to first 200 nucleotides in the mRNA from
the BRCA1 gene in cow, wolf, chicken, human, macaque, mouse, chimpanzee, and rat...
Done! Written to the file, brca1-alignment.fasta . :-)
Bye!
```

## Project 2 - How to run
### `global_linear` and `global_affine`
The `global_linear` and `global_affine` programs take the following arguments:
* Sequence one, provided as a string or as a path to a FASTA file.
* Sequence two, provided as a string or as a path to a FASTA file.
* For linear:
    * A gap cost, provided in the form of a number.
* For affine:
    * alpha of the cost function, provided in the form of a number.
    * beta of the cost function, provided in the form of a number.
* A scoring matrix, provided as a path to a file containing the matrix writtin in Phylip format (note that our program uses `tab`s as separators between the numbers).

Example scoring matrix: 
```
4
A	0	5	2	5
C	5	0	5	2
G	2	5	0	5
T	5	2	5	0
```
You can run case 1 from project2_examples.txt with the following line in the terminal (from our `project2/` directory):
```
python .\global_linear.py .\examples\seq1_case1.fasta .\examples\seq2_case1.fasta 5 .\examples\example_score_matrix.txt
```

The program will compute the optimal alignment cost of the given sequence under the given parameters and output the resulting matrix. Then, it will ask the user if they would like to find an optimal alignment, which can be answered with a `y` for yes or `n` for no. 
If computed, the optimal alignments will be written to the file `output/linear_alignments.fasta` or `output/affine_alignments.fasta`.

An example run of the linear program looks like so (note that this is run while navigated to the `project2/` directory):

```
> python .\global_linear.py .\examples\seq1_case1.fasta .\examples\seq2_case1.fasta 5 .\examples\example_score_matrix.txt

Hi! Welcome the the affine version of the program. :)

--- Computing cost ---
[[ 0.  5. 10. 15. 20. 25. 30. 35. 40. 45. 50. 55. 60.]
 [ 5.  0.  5. 10. 15. 20. 25. 30. 35. 40. 45. 50. 55.]
 [10.  5.  0.  5. 10. 15. 20. 25. 30. 35. 40. 45. 50.]
 [15. 10.  5.  0.  5. 10. 15. 20. 25. 30. 35. 40. 45.]
 [20. 15. 10.  5.  0.  5. 10. 15. 20. 25. 30. 35. 40.]
 [25. 20. 15. 10.  5.  5.  5. 10. 15. 20. 25. 30. 35.]
 [30. 25. 20. 15. 10.  7. 10.  5. 10. 15. 20. 25. 30.]
 [35. 30. 25. 20. 15. 10. 12. 10. 10. 15. 15. 20. 25.]
 [40. 35. 30. 25. 20. 15. 12. 15. 10. 12. 17. 20. 20.]
 [45. 40. 35. 30. 25. 20. 17. 17. 15. 12. 17. 22. 20.]
 [50. 45. 40. 35. 30. 25. 22. 19. 20. 17. 12. 17. 22.]
 [55. 50. 45. 40. 35. 30. 25. 24. 21. 20. 17. 17. 19.]
 [60. 55. 50. 45. 40. 35. 30. 25. 26. 25. 22. 17. 22.]]

Show optimal alignment? [y/n]
> y

--- Computing optimal alignment ---
Done! Written to output file.

Bye! :) 
```



### Tests and experiments
```
python .\test.py
```

```
python .\experiments.py
```

```
python .\eval.py
```

Note that you need to navigate to the `project2/` directory for it to work, since the tests and experiments use file paths relative to this folder. 

For the `experiment.py` program, the constructed plots will pop up and be saved in the `figures/` folder, e.g.:

<img src="project2/figures/linear.png" width="400" />