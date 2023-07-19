# Algorithms in Bioinformatics 2023

**Group**: Group 6

**Group members**: 
- Anastasia Kratschmer (201907407)
- Noa Van de Velde (202209940) 
- Stinna Danger (201907211)
- Anna Kathrine Skov (201905355)

## Project 4 - How to run
*Note that you should be in the `project4/` directory.* 

### `rfdist`
The program takes two file paths each to an evolutionary tree in Newick format.
It can be run ffrom the terminal with the following command.
```
python .\rfdist.py .\testdata\tree1.new .\testdata\tree2.new
```
The program will compute and output the RF distance of the two trees. The above run would look like so:

```
Beep boop, dubbedy doop!

Computing RF distance...
Done!

The distance is 8.
```

If needed, the script can be altered to draw the two given trees in the terminal.


### Experiments
You can run the experiments by running the experiments.py file. 
The results for the first, second and third experiment can be found in the files normaldist, permuteddist and exp3 respectively. 