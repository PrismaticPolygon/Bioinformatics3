import numpy as np

# https://science.umd.edu/labs/delwiche/bioinform/lec/AlignHeuristic.html
# http://etutorials.org/Misc/blast/Part+III+Practice/Chapter+5.+BLAST/5.1+The+Five+BLAST+Programs/
# https://www.cs.cmu.edu/~02710/Lectures/SeqAlign2015.pdf
# file:///C:/Users/user/Downloads/els_blast_2003.pdf
# https://bio.libretexts.org/Bookshelves/Cell_and_Molecular_Biology/Book%3A_Investigations_in_Molecular_Cell_Biology_(O'Connor)/9%3A_Protein_Conservation/9.3%3A_BLAST_algorithms_are_used_to_search_databases
# http://mirror.ufs.ac.za/blast/demo/ieee_talk.pdf
# http://resources.qiagenbioinformatics.com/manuals/clcmainworkbench/800/index.php?manual=How_does_BLAST_work.html
# http://www.gersteinlab.org/courses/452/09-spring/pdf/Altschul.pdf
# https://is.muni.cz/el/1431/jaro2015/C2135/um/55046697/BLAST_algorithm.pdf

# A heuristic procedure that runs in sub-quadratic time (similar to FASTA and BLAST)
# Assessed on trade-offs between running time and quality of output.
# A typical input will come with a planted alignment, which consists of segments of matches of different lengths
# separated by random stuff

# BLAST: attempts to find a short fragment of a query sequence that aligns perfectly with a fragment
# of a subject sequence in a database,
# Runs in sub-quadratic time.
# BLAST requires a hella huge database to actually work.

def heuralign(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    BLAST(s, t)

    return 0, [], [], "", ""

# Construct neighbourhood words by iteratively replacing letters in a word with characters from the alphabet

# Word: a set-length string of adjacent letters in a sequence
# Hit: given two sequences, and a threshold T, a pair of words, one from each sequence, whose aligned score > T
# A segment pair is locally maximal if its score cannot be improved either by extending or shortening both segments

# MSP score can be calculated in O(n^2) using DP.

# In searching a DB of thousands of sequences, generally only a handful will be homoglous.
# So we're interested only in those sequence entries with MSP scores about some threshold T.

# Search a long sequence of all occurrences of certain short sequences
# w = 4; map each word to an integer between 1 and 20^4.
# So a word can be used as an index into an array of size 20^4
# Let the ith entry of such an array point to a list of all occurences in the query sequence of the ith word.


# Let a word pair be a segment of fixed length w.
# BLAST seeks only segment pairs containing a word pair with a score of at least T.
# Can quickly determine whether it contains a word of length w that can pair with the query sequence to produce
# a word pair with a score greater than or equal to threshold T.
# For each word in the query sequence, make a list of all neighbouring 'words'
# that score above some threshold. Scan the database for these words
# Perform ungapped hit extension until score < threshold


def BLAST(query, target):

    global ALPHABET

    alphabet = ALPHABET[:-1]

    w = 3               # Word length. 3 for protein sequences, 11 for nucleic acids
    T = (w * 2)         # Word score threshold. Typically 2 bits per residue

    # So we could preallocate. But do we need to?
    # No, of course not.

    words = dict()

    for i in range(len(query) - w + 1):

        query_word = query[i:i + w]

        if query_word not in words:

            words[query_word] = [i]

        else:

            words[query_word].append(i)

    # Having done that, we look for keys in common, and then use those as seeds if they exceed
    # a certain threshold.

    # Then for every word we find matches about a certain threshold?
    # So: we make a list of all the words that we have in the query string
    # This is an O(n) operation. For each, we store its index.
    # That is correct.
    # It's a dictionary.
    # Then we do the same thing with the target string.

    # The starting point is the words that two sequences have in common.

    # The starting point is WORDS that the two sequences have in common.
    # So if we can find such a hit...

    # Okay. Now we need to think. How should this work?
    # Build a look-up table in sub-quadratic time.

    # The query sequence is broken in 'words' that will act as seeds in alignments
    # A look-up table is made of all the 'words' (short subsequences_ in the query sequence
    # In many types of searches, neighbouring words are included


    # BLAST searches for matches (or synonyms) in target entries in the database

    # But that's quadratic!

    print(query, target)

    # For each word in the query query sequence, a compilation of neighbourhood words which exceed
    # the threshold of T are also generated.
    # NEIGHBOURHOOD WORD: a word obtaining a score of at least T when comparing.

    # Then we compare against all words in the database. So we're already in quadratic territory.
    # List all possible matching words, then filter by threshold.
    # A table of all possible seeds, along with their start points.
    # And then we compare against THOSE?
    # Is that really local? It seems a bit of a gimmick

    #

    # Okay, nice.
    seeds = []

    for i in range(len(query) - w + 1):

        query_word = query[i:i + w]
        hit_word = target[i:i + w]
        word_score = score(query_word, hit_word)

        print("score({}, {}) = {}".format(query_word, hit_word, word_score))

        if word_score > T:

            seeds.append((i, query_word))

    print(seeds)

    # It then scans through the database, and whenever it finds a word in this set,
    # it starts a hit extension process to extend the possible match as an ungapped alignment in both directions
    # stopping at the maximum scoring extension.

    # FASTA cares about all the common words in the database;
    # BLAST only cares about high-scoring words.
    # Compare the word in the list in step2 with ALL three letter words.

    # Ah. The words are three long.

    # For each word from the query sequence, find the list of words with high score using a substitution matrix
    # matched against what?

    # print("Word count: ", word_count)



def score(s, t):

    score = 0

    for i in range(len(s)):

        score += cost(s[i], t[i])

    return score



def cost(c1, c2):

    global ALPHABET, SCORING_MATRIX

    i, j = ALPHABET.index(c1), ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]
