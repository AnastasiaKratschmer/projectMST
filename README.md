# projectMST
Ideas:
* Try to check if the order of alignments in 1/2 approx matters for the score. Do this by sorting the sequences based on a min spanning tree you get from the dist matrix in the linear al.s. ``` ! DID NOT WORK... GOTTA THINK OF A NEW IDEA THEN!```
* Make MST, use the standard traversals df and bf to guide merge. use one of the most extreme differences in the distance matrix as a starting point for traversal order. df works better than bf for some reason, but the traversal order seems logically off. ```KINDA WORKS BUT RESULTS STILL WORSE THAN WITH JUST CENTER STRING, BUT ONLY SLIGTLY. (28/7 2023)```
* Using the idea from above: now just making my own traversal order very much like bf, just more handheld and where absolutely everything is aligned to its own closest match. ```A GOOD LOOKING TRAVERSAL ORDER IS BEING SPIT OUT, NOW I JUST NEED TO INCORPORATE IT (3/8 2023).```

```diff
- text in red
+ text in green
! text in orange
# text in gray
@@ text in purple (and bold)@@
```
