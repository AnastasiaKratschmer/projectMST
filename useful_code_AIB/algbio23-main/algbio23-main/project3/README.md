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