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

> Full credit will be given to solutions that run in quadratic time, while the trivial solution, which runs in cubic time,
will count for very little (and thus a partially correct quadratic-time solution will most likely be given a higher mark
than a fully correct but cubic-time one).

Use the two rows technique to find the cell with the maximum value in the table.
Align the reverse string starting from the indices of the cell with the maximum value
Then we can find the start of the best local alignment between the two strings
Use hirschberg from the start to the end to find the traceback.


Rumour has it that it can be done in cubic time. The intuition is: find the end, find the start, and go global between 
those two points. Then we backtrack using Waterman.
Fuck it: if it works, then I'm game.

Aligning two sequences within a specified diagonal band (Chao, Pearson, Miller)


> A heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [20 marks] 

Implemented a simple variant of [FASTA](https://en.wikipedia.org/wiki/FASTA), with `ktup = 4`.

## Part 2

> Design your own substitution-cost function that operates on pairs of sequences of letters instead of on pairs of 
letters. Clearly describe it on at most one page  [15 marks], e.g.
> * give a cost of a multi-letter substitution of `ABC` by `CBB` which is different from the simple
addition of the single-letter costs, i.e. the mismatches `(A, C)`, `(C, B)` and the match `(B, B)`
> * have a fixed cost for inserting / deleting a sequence of `C` irrespective of its length, e.g. the first `C`
incurs a cost of `1\2`, the second `1\4`, and so on.

The lectures slides available [here](https://www.site.uottawa.ca/~lucia/courses/5126-10/lecturenotes/03-05SequenceSimilarity.pdf)
introduce a "general gap penalty" `g(q)` for a gap of length `q`, based on the idea
that mutations (which may cause the insertion of deletion of a substring) may be as likely
as the indel of a single base.

#### Ideas
* Mismatches of the same character are not scored poorly. 

`g(q) = h + qs`, where `h` is the penalty for initiating the gap, and `s` is a penalty
proportional to the length of the gap.



