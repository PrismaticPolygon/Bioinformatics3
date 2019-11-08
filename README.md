# Functions
* Basic DP (quadratic time and space) [50]
* DP (linear space) [15]
* Heuristic (sub-quadratic) [20]


Should be called dynprog, dynproglin, and heuralign
Signature should be: (alphabet, scoring_matrix, sequence1, sequence2)

Return should be a list of the score of the bst local alignment plus two lists
of indices, one for each sequence, that realises this score. 

Marked automatically. 

For Q2: read the linear-space algorithm in the textbook, which
solves the global-alignment problem. To solve the local, modify
this to work recursively:
* It is easy to find the end-point of an optimal local alignment
* If the corresponding starting-point can be found, we're done
* If the starting-point is the the left, descend into that side;
otherwise, descend into the right. 
* In this manner some of the recursive subproblems become global alignment ones,
while others become local alignment with a given endpoint. 
