> [!NOTE]
> Go to /FINAL_code to see actually finished code :)

You can try aligning strings using either Gusfield's Sum Of Pairs Approximation Algorithm for multiple sequence alignment, or modifications of the algorithm which use Kruskal and Prim minimum spanning trees as an alignment guide tree.  

Download `FINAL_code/align_strings.py` to get going.  
If you don't have a cost matrix but still want to try aligning, download `FINAL_code/score_matrix.txt`.  
If you don't have any sequences to align but still want to try aligning, download  `FINAL_code/testseqs.fasta`.  

The program is run from the terminal (tested in Git Bash).
Example usage: 
<pre>
$ python3 align_strings.py testseqs.fasta 5 score_matrix.txt False g
</pre>
The inputs are Python program to use, alignment script to run, sequences to align, gap cost, score matrix, perform integrity checks (True/False) and lastly type of guide tree used for alignment: k: Kruskal, p: Prim, g: Gusfield.   
Once the code has run you will be asked `Do you want to save the output to a csv? (Y/n)`. Write `Y` for yes or `n` for no in the terminal. If you say yes, you'll be asked to provide a name for the output csv file.

See example output in `FINAL_code/output_sample.py`. (Notice that this file contains output from 3 seperate runs of the script!)

𝑬𝓃𝓳𝑜𝔂 :) 
