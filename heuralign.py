import numpy as np

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

    return 0, [], []

def BLAST(query, target):

    w = 3               # Word length. 3 for protein sequences, 11 for nucleic acids
    T = (w * 2)    # Word score threshold. Typically 2 bits per residue

    # That's the idea behind the table

    # The starting point is WORDS that the two sequences have in common.
    # So if we can find such a hit...

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
