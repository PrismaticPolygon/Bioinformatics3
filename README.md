# Bioinformatics

Implement different algorithms for local alignment.

## Part 1

> Basic dynamic programming that runs in quadratic time and space [50 marks]

Implemented Smith-Waterman based off the Needleman-Wunsch algorithm described 
[here](http://biorecipes.com/DynProgBasic/code.html). Quadratic time and space.

> Dynamic programming that runs in linear space [15 marks]

Implemented Hirschberg, described [here](https://en.wikipedia.org/wiki/Hirschberg's_algorithm).
Quadratic time, linear space. However, this is a *global* alignment algorithm. Very few
linear local alignment algorithms exist; two options are:
* A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons (Waterman and Eggert, 1987)
 (available [here](https://www.sciencedirect.com/science/article/pii/0022283687904785)).
* Huang and Miller (1991) (available [here](https://www.sciencedirect.com/science/article/pii/019688589190017D)).

> A heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [20 marks] 

Implemented a simple variant of [FASTA](https://en.wikipedia.org/wiki/FASTA), with `ktup = 4`.

## Part 2

> Design your own substitution-cost function that operates on pairs of sequences of letters instead of on pairs of 
letters. Clearly describe it on at most one page. [15 marks]