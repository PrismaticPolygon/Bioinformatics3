# Part 2

> Design your own substitution-cost function that operates on pairs of sequences of letters instead of on pairs of 
letters. Clearly describe it on at most one page. [15 marks]

This substitution-cost function will operate on the standard DNA alphabet, i.e `ACGT`. The following rules are not
defined with a biological purpose in mind. They are applicable to sequences of letters of any length.

## Rules

### 1. The cost of mismatch of a sequence of the same characters of length `n` is `n!`. 

For example, matching `AAA` against `CGG` incurs of a cost of `3!`. This penalty is independent of the cost of matching 
the characters in each subsequence with each other. 

Naturally, this heavily penalises such sequences, with the result that matching characters
with gaps is significantly cheaper: matching `AAA___` against `___CDD`, for instance, would incur a much lower penalty 
(as governed by the next rule).

### 2. The cost of indel of a sequence is logarithmic.

For example, matching `AGT` against `___` would incur a cost of `4 + log(3)`. More generally, the cost of indel of a
sequence `q` is `h + log(|q|)`, where `h` is a constant, chosen here to be 4. The cost of indel is therefore
significantly lower for large sequences than in a standard linear cost function.

### 3. The score of aligning a sequence is the product of the scores of aligning individual characters.

For example, matching `ATG` with `ATG` would score `(score(A, A) * score(T, T) * score(G, G)`. The score of matching
now increase exponentially, rather than linearly, and is consequently much more important.

Interestingly, this also means that aligning, for example, `AC` with `AG` would likely result in a negative score 
(i.e. a cost). A sequence containing two mismatches and one match would become positive again. 

## Application

As these rules are applicable to sequences of any length, it is not possible to present a one-size-fits-all substitution
matrix. Examples for sequence of length `2` are presented below instead, using the following 
one-character: substitution matrix:

|   | A  | C  | G  | T  | _  |
|---|----|----|----|----|----|
| A |  2 | -1 | 1  | -1 | -2 |
| C | -1 | 2  | -1 | 1  | -2 |
| G | 1  | -1 | 2  | -1 | -2 |
| T | -1 | 1  | -1 | 2  | -2 |
| _ | -2 | -2 | -2 | -2 | 0  |
|   |    |    |    |    |    |
 
 In the event that rules are combined, scores / costs are additive.

* `(AA, AA)`: rule 3. `2 * 2 = 4`
* `(G_, TA)`: rules 1 and 2. `4 + log(1) + 1 = 6` 
* `(__, __)`: rule 3. `0 * 0 = 0`
* `(TC, TT)`: rule 3. `2 * 1 = 1`
* `(CC, __)`: rule 2. `4 + log(2) = 4.301...`


##### Sources for ideas

* [Slides](https://www.site.uottawa.ca/~lucia/courses/5126-10/lecturenotes/03-05SequenceSimilarity.pdf)
* [More slides](https://www.cs.cmu.edu/~02710/Lectures/ScoringMatrices2015.pdf)