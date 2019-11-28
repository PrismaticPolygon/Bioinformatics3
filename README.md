# Bioinformatics

Implement different algorithms for local alignment.

## Part 1

> Basic dynamic programming that runs in quadratic time and space [50 marks]

Implemented Smith-Waterman, based off the Needleman-Wunsch algorithm described 
[here](http://biorecipes.com/DynProgBasic/code.html). Quadratic time and space.

> Dynamic programming that runs in linear space [15 marks]

Implemented a two-pass local alignment algorithm, which identifies the start and end points of the best local
alignment between two strings. Hirschberg's algorithm, described [here](https://en.wikipedia.org/wiki/Hirschberg's_algorithm),
is then used for find the optimal global alignment within that section.

> Full credit will be given to solutions that run in quadratic time, while the trivial solution, which runs in cubic time,
will count for very little (and thus a partially correct quadratic-time solution will most likely be given a higher mark
than a fully correct but cubic-time one).

To the best of my understanding, my algorithm, as currently implemented. is quadratic.

> A heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST) [20 marks] 

Implemented a simple variant of [FASTA](https://en.wikipedia.org/wiki/FASTA), with `ktup = 4`.

# How does FASTA work?

Pseudo-code, or even detailed information, on how FASTA works is simply not available. I have compiled the following
from a variety of sources, listed below.

### 1. Identify hot-spots

Hot-spots are exact matches of length `ktup` between a query and a database string. This stage can be done using
a look-up table or a hash. Pre-process the database and store the location of each possible `ktup`. Move a 
`ktup`-sized window along the query sequence and record the position of matching locations in the database.

Spaces between hot-spots are allowed.

Determine all exact matches of length `k` between two sequences (hot-spots). A hot-spot is given by `(i, j)`, where 
`(i, j)` are the locations (i.e. start positions) of an exact match of length `k` in the query and database sequence 
respectively. Any such hotspot lies on the diagonal `(i - j)` of the dot-plot. The main diagonal has a number 0. 
Diagonals above the main one have positive numbers, and ones below, negative.

It's not clear if we should merge these, or score them.

Look for `hot-spots`.

### 2. Identify the 10 best diagonal runs

A `diagonal run` is a set of hot-spots that lie in a consecutive sequence on the same diagonal. They correspond to 
gapless local alignment. Each run can be indexed by the offset `y - x`, such that the main diagonal has an offset of 
`0`. Diagonals below the main one have a negative index, and those above a positive index.

Each diagonal run is assigned a score. A positive score is given to each match, and a negative score for gaps in the run,
dependent on length.

and top 10 proceed to the next step.

A score is assigned to each diagonal run. This is done by giving a positive score to each 
match, and a negative score for gaps in the run. The latter scores decrease with increasing length of gaps between 
hot-spots. The algorithm then locates the 10 best diagonal runs.

Diagonal sums are computed into a table, indexed with the offset, and initialised to `0`. Go through the `k`-words of 
`l`. Look for matches and  update the diagonal sums.
Ah. So that's just the sum of the actual hits.

Index diagonal `l` by the offset `y - x`

In the first step, the 10 best diagonal regions are found using a simple formula
based on the number of ktup matches and the distance between the matches without considering
shorter runs of idenities, conservative replacements, insertions, or deletions.

Find the best `diagonal runs`. Each `hot spot` gets a positive score. The distance between hot spots is negative
 and length dependent. Store the 10 best `diagonal runs`.

Local regions are diagonals of a certain length in a dot plot, evaluated by counting matches and penalising intervening mismatches (without using a score matrix)

### 3. Score the 10 best diagonal runs.

Each of the 10 diagonal runs with the highest score are further processed. Within each of these runs, an optimal local 
alignment is computed using the match score substitution matrix. These alignments are called initial regions. The score 
of the best sub-alignment found in this phase is `init_1`.

Each high-scoring diagonal chosen in the previous step is rescored to find subregions with identities shorter than `k`. 
Non-matching ends of the diagonal are trimmed.

Compute `init_1` and filter. Diagonal runs specify a potential alignment. Evaluate each run using a substitution
matrix. Define the best scoring run as `init_1`. Discard any much lower scoring runs.

In the second step, these 10 regions are rescored using a scoring matrix that allows
conservative replacements and runs of idenities shorter than ktup to contribute to the similarity
score. For each of these best diagonal regions, a subregion with maximal score is identified.
is identified: the initial region.

Rescan the regions taken using the scoring matrices. Trim the ends of the region to include only those contributing to the highest 
score. While rescoring, conservative replacements that contribute to the similarity score are taken.
A subregion with maximum score is identigied. Te initial scores are used to rank the library sequences.
The highest score is called `init_1`

### 4. Combine high-scoring sub-alignments



Then, combine high-scoring sub-alignments into a single larger alignment, allowing the introduction of gaps
into the alignment. The score of this alignment is `init_n`.

Then, calculate an optimal alignment of intiail regions as a combination of maximal regions with maximal score
Use the resulting score to rank the library sequences (?). Limit the degradation of selectivity by including
in this step only those initial regions with scores above a threshold.

Combine diagonal runs and compute `init_n`. Take the good alignments from the previous stage.
Allow gaps / indels. Combine into a single, better-scoring, alignment Construct a directed weighted graph,
where vertices are the runs, and edge weights represent gap penalties. Find the best path through the graph: `init_n`.

In an alignment, if several regions with scores greater than a CUTOFF value are found,
check whether the trimmed initial regions can be joined to form an approximate alignment with gaps.
Calculate a similiarity score that is the sum of the joined regions, penalising each gap 20 points.
This initial similarity score `init_n` is used to rank the library sequences. 

### 5. Use banded DP to produce an optimal local alignment along the best matched regions.

Finally, used banded SW to produce an optimal local alignment along the best matched regions.
Center of the band is determined by the region with the score `init_1`, and the band has width `8`.
The score of the resulting alignment is `opt`.

Find the best local alignment. 

The score of the single best initial region is `init_1`

Use a cut-off value that is approximately one standard-deviation above the average score
expected from unrelated sequences in the library. Use banded SW for each alignment of query sequences in a database sequence.
Use a band of 32 residues centered on the `init_1` region of step2 for calcilating.

Find the best local alignment. Use the alignments from the previous stage to define a narrow band through the search space.
Go through that band using a DP approach. The best local alignment is called `opt`.

Finally, align the highest scoring library sequences using a modification of the optimisation method
described by NeedlemanWunsch / SmithWaterm. This comparison considers all possible alignments of the queyr
and library sequence that fall within a band centred around the highest-scoring initial region.

Use a band of 32 residues centered on the `init_1` region of step2 for calcilating.


####Sources

* [Slides](http://www.inf.ed.ac.uk/teaching/courses/bio2/lectures05/lecture2.pdf)
* [Wikipedia](https://en.wikipedia.org/wiki/FASTA)
* [More slides](https://www.cs.helsinki.fi/bioinformatiikka/mbi/courses/07-08/itb/slides/itb0708_slides_83-116.pdf)
* [Original paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC280013/pdf/pnas00260-0036.pdf)
* [Even more slides](https://ksvi.mff.cuni.cz/~mraz/bioinf/BioAlg10-8.pdf)
* [A website](http://www.biology.wustl.edu/gcg/fasta.html#algorithm)
* [StackOverflow](https://stackoverflow.com/questions/8366581/fasta-algorithm-explanation)
* [Lecture notes](http://www.cs.tau.ac.il/~rshamir/algmb/98/scribe/pdf/lec03.pdf)