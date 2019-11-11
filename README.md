# Functions
* Basic DP (quadratic time and space) [50]
* DP (linear space) [15]
* Heuristic (sub-quadratic) [20]


Should be called dynprog, dynproglin, and heuralign
Signature should be: (alphabet, scoring_matrix, sequence1, sequence2)

Return should be a list of the score of the bst local alignment plus two lists
of indices, one for each sequence, that realises this score. 

Marked automatically. 

BLAST is mostly involved in finding ungapped, locally optimal sequence alignments

### BLAST
* True match alignments are likely to contain somewhere within them a short stretch of identities
(very high scoring matches)
* These can be used as 'seeds' from which to extend out in seach for longer alignments
* By keeping seed sequences short, the query sequence can be pre-processed and a table of all possible seeds
generated

