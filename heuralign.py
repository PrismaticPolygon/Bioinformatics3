import numpy as np
from main import traceback

# Nah. Let's just use what we have now.

def get_lookup(s, ktup):

    lookup = dict()

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]

        if word not in lookup:

            lookup[word] = [y]

        else:

            lookup[word].append(y)

    return lookup

def get_diagonals(alphabet, scoring_matrix, s, t, ktup, lookup, threshold):

    diagonals = dict()

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        if word in lookup:

            for y in lookup[word]:

                diagonal = x - y

                if diagonal not in diagonals:

                    seed = extend(alphabet, scoring_matrix, s, t, (x, y, ktup), threshold)

                    diagonals[diagonal] = [seed]

                else:

                    # diagonals[diagonal].append((x, y, ktup))

                    last_x, last_y, last_length = diagonals[diagonal][-1]  # Get the last entry on this diagonal

                    if last_x + last_length >= x and last_y + last_length >= y:  # If they overlap, update the old one

                        new_length = x - last_x + ktup

                        seed = extend(alphabet, scoring_matrix, s, t, (last_x, last_y, new_length), threshold)

                        diagonals[diagonal][-1] = seed
                    #
                    elif (x, y) not in diagonals[diagonal]:

                        seed = extend(alphabet, scoring_matrix, s, t, (x, y, ktup), threshold)

                        diagonals[diagonal].append(seed)

    return diagonals

# http://www.cs.tau.ac.il/~rshamir/algmb/98/scribe/html/lec03/node2.html

# The docs really SUCK, huh?
# There's something too.
# Combine diagonal runs.
# Allow gaps / indels.

# Construct a directed, weighted graph.
# Vertices are the runs.
# Edge weights represent gap peanlites.
# BLAST is the fastest. Apparently.
# Fuck it
#

# Find all k-length identities.
# Then find locally similar regions by selecting those dense with identities.
# Recore the regions by applying a substitution matrix.
# Take the 10 best initial regions
# Create an alignment of the trimmed initial regions using DP with a gap penalty of 20.
# Exclude regions with too low a score.
# Optimise the alignment using banded DP.

def best_diagonal(alphabet, scoring_matrix, s, t, diagonals):

    # best_x, best_y = None, None
    best_diagonal = None
    best_diagonal_score = 0

    # We should use score here. Not length.

    for diagonal, seeds in diagonals.items():

        # That's wrong.
        # Oh, I see. Fair enough!
        # The calculation is... close.
        # Wait. That's right.
        # That is an unacceptably poor result.
        # Lol.
        # FUCK.

        start_x, start_y, _ = seeds[0]
        end_x, end_y, length = seeds[-1]

        # print("seed_s", s[start_y: end_y + length])
        # print("seed_t", t[start_x: end_x + length])

        diagonal_score = score(alphabet, scoring_matrix, s[start_y: end_y + length], t[start_x: end_x + length])

        # print(diagonal, diagonal_score)

        if diagonal_score > best_diagonal_score:

            best_diagonal = (diagonal, start_x, start_y, end_x + length, end_y + length)
            best_diagonal_score = diagonal_score

    return best_diagonal


def heuralign(alphabet, scoring_matrix, s, t):

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)

    ktup = 2    # 6 for DNA sequence matching (|alphabet| = 4). 2 for protein sequence (|alphabet| = 20)
    threshold = get_threshold(alphabet, scoring_matrix, s, t)

    lookup = get_lookup(s, ktup)

    diagonals = get_diagonals(alphabet, scoring_matrix, s, t, ktup, lookup, threshold)

    # The hardest part is combining the diagonals, right?
    # How tf is he going to check anyway/

    print(diagonals)

    diagonal, start_x, start_y, end_x, end_y = best_diagonal(alphabet, scoring_matrix, s, t, diagonals)
    #
    print(diagonal, start_x, start_y, end_x, end_y)
    #
    s = s[start_y: end_y]
    t = t[start_x: end_x]

    D, _ = score_matrix(alphabet, scoring_matrix, s, t, band_width=12, threshold=0)
    #
    score, s_matches, t_matches, s_align, t_align = traceback(alphabet, scoring_matrix, s, t, D, local=True)

    return score, s_matches + start_x, t_matches + start_y, s_align, t_align
    #
    # return 0, 0, 0, 0, 0

def score(alphabet, scoring_matrix, s, t):

    score = 0

    for i in range(min(len(s), len(t))):

        score += cost(alphabet, scoring_matrix, s[i], t[i])

    return score

def score_seed(alphabet, scoring_matrix, s, t, seed):

    x, y, length = seed
    score = 0

    for i in range(length):

        score += cost(alphabet, scoring_matrix, s[y + i], t[x + i])

    return score

def evaluate(alphabet, scoring_matrix, s, t, diagonal):

    if diagonal < 0:    #

        x, y = 0, abs(diagonal)
        limit = len(s) - y

    else:

        x, y = diagonal, 0
        limit = len(t) - x

    score = 0

    for i in range(limit):

        score += cost(alphabet, scoring_matrix, s[y + i], t[x + i])

    return score






    # If we have -2, for example.

    # -2 = x - y
    # As we know one of them must be 0.




    # x, y, length = seed

    # So we keep going back until one of them is 0.
    # Is that a hard or easy problem? I can't even tell.

def extend(alphabet, scoring_matrix, s, t, seed, threshold):

    x, y, length = seed

    seed_score = score_seed(alphabet, scoring_matrix, s, t, seed)

    # So. For each seed, we evaluate the score of that diagonal.

    up, down = 1, 1

    while up or down:

        if x == 0 or y == 0:

            up = 0

        else:

            new_score = seed_score + cost(alphabet, scoring_matrix, t[x - 1], s[y - 1])

            if new_score > threshold:

                x -= 1
                y -= 1
                length += 1

                seed_score = new_score

            else:

                up = 0

        if x + length == len(t) or y + length == len(s):

            down = 0

        else:

            new_score = seed_score + cost(alphabet, scoring_matrix, t[x + 1], s[y + 1])

            if new_score > threshold:

                length += 1
                seed_score = new_score

            else:

                down = 0

    # print("Done", x, y, length)
    # print("Score", seed_score)

    return x, y, length

def get_threshold(alphabet, scoring_matrix, s, t):

    char_freq = {i: (s + t).count(i) / len(s + t) for i in alphabet}

    pairs = [(i, j) for i in alphabet for j in alphabet]

    pair_weights = [char_freq[x] + char_freq[y] for (x, y) in pairs]

    return np.average([cost(alphabet, scoring_matrix, x, y) for (x, y) in pairs], weights=pair_weights)

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def rev(l):

    return l[::-1]

# It is local, and it is not sublinear.

def score_matrix(alphabet, scoring_matrix, s, t, band_width=None, threshold=0):

    print(s, t)

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    D = np.zeros((size_y, size_x), dtype="int16")

    max_i = None  # The index of the cell containing the maximum score
    max_score = None  # The maximum score observed

    for y in range(size_y):  # y

        min_x = 0 if band_width is None else y - band_width  # If banded, the minimum x that we explore
        max_x = size_x if band_width is None else y + band_width + 1  # If banded, the maximum x that we explore.

        for x in range(min_x, max_x):

            if 0 < x < size_x and 0 < y:

                D[y, x] = max(
                    D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),
                    D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"),
                    D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_"),
                    threshold
                )

                if max_i is None or D[y, x] > max_score:

                    max_i = y, x
                    max_score = D[y, x]


    return D, max_i


# alphabet = "ABCD_"
# scoring_matrix = [
#             [ 1,-5,-5,-5,-1],
#             [-5, 1,-5,-5,-1],
#             [-5,-5, 5,-5,-4],
#             [-5,-5,-5, 6,-4],
#             [-1,-1,-4,-4,-9]]
# # s = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
# # t = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
#
# s = "ABCDABC"
# t = "ABCAABC"

# heuralign(alphabet, scoring_matrix, s, t)

#
# alphabet = "ABC_"
# scoring_matrix = np.array([
#             [1, -1, -2, -1],
#             [-1, 2, -4, -1],
#             [-2, -4, 3, -2],
#             [-1, -1, -2, 0]
# ])
# s = "AABBAACA"
# t = "CBACCCBA"
#
# heuralign(alphabet, scoring_matrix, s, t)