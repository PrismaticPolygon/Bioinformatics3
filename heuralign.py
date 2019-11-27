import numpy as np

def score_matrix(alphabet, scoring_matrix, s, t, band_width=None, local=True, sublinear=False, threshold=0):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    shape = (2, size_x) if sublinear else (size_y, size_x)  # The shape of the score matrix. If sublinear, we only use 2 rows.

    D = np.zeros(shape, dtype="int16")

    max_i = None    # The index of the cell containing the maximum score
    max_score = None    # The maximum score observed

    for y in range(size_y):  # y

        min_x = 0 if band_width is None else y - band_width # If banded, the minimum x that we explore
        max_x = size_x if band_width is None else y + band_width + 1  # If banded, the maximum x that we explore.

        for x in range(min_x, max_x):

            if not local and y == 0 and x > 0:    # First row: the cost of matching t with all gaps. If local, should be 0

                D[0, x] = D[0, x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_")

            elif not local and x == 0 and y > 0:  # First column: the cost of matching s with all gaps. If local, should be 0

                D_y = 1 if sublinear else y

                D[D_y, 0] = D[D_y - 1, 0] + cost(alphabet, scoring_matrix, s[y - 1], "_")

            elif 0 < x < size_x and 0 < y:

                D_y = 1 if sublinear else y

                D[D_y, x] = max(
                    D[D_y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),  # The cost of matching two characters
                    D[D_y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"), # The cost of matching a gap in s with a character in t
                    D[D_y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_") # The cost of matching a gap in t with a character in s
                )

                D[D_y, x] = 0 if (local and D[D_y, x] < threshold) else D[D_y, x]

                if max_i is None or D[D_y, x] > max_score:

                    max_i = y, x
                    max_score = D[D_y, x]

        if y > 0 and sublinear:

            D[0] = D[1].copy()

    return (D[1] if sublinear else D), max_i

def build_lookup(s, ktup):

    lookup = dict()

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]

        if word not in lookup:

            lookup[word] = [y]

        else:

            lookup[word].append(y)

    return lookup

def build_lookup_II(t, ktup, lookup):

    lookup_II = dict()

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        if word in lookup:

            if word not in lookup_II:

                lookup_II[word] = [x]

            else:

                lookup_II[word].append(x)

    return lookup_II

def build_diagonals(alphabet, scoring_matrix, s, t, ktup, lookup, threshold):

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

def best_diagonal(s, t, diagonals):

    best_x, best_y = None, None
    best_diagonal_length = 0

    for diagonal, seeds in diagonals.items():

        diagonal_length = 0
        diagonal_x, diagonal_y = len(s), len(t)

        for seed in seeds:

            x, y, length = seed

            diagonal_length += length

            if x < diagonal_x:

                diagonal_x = x

            if y < diagonal_y:

                diagonal_y = y

        # print("{}: {} ({}, {})".format(diagonal, diagonal_length, diagonal_x, diagonal_y), end="\n")

        if diagonal_length > best_diagonal_length:

            best_x = diagonal_x
            best_y = diagonal_y
            best_diagonal_length = diagonal_length

    return best_x, best_y

# Find the 10 best diagonal runs.
# Give each hot spot a socre, and give the space between between consecut

def heuralign(alphabet, scoring_matrix, s, t):

    scoring_matrix = np.array(scoring_matrix)
    alphabet += "_"
    ktup = 2
    threshold = get_threshold(alphabet, scoring_matrix, s, t)

    lookup = build_lookup(s, ktup)

    print(lookup)

    lookup_ii = build_lookup_II(t, ktup, lookup)

    print(lookup_ii)

    diagonals = build_diagonals(alphabet, scoring_matrix, s, t, ktup, lookup, threshold)

    # If I did it Peter's way...

    print(diagonals)

    best_x, best_y = best_diagonal(s, t, diagonals)

    s = s[best_x:]
    t = t[best_y:]

    D, _ = score_matrix(alphabet, scoring_matrix, s, t, band_width=10, threshold=threshold)

    score, s_matches, t_matches, s_align, t_align = traceback(alphabet, scoring_matrix, s, t, D)

    return score, s_matches + best_x, t_matches + best_y, s_align, t_align

def traceback(alphabet, scoring_matrix, s, t, D, local=True):

    score = np.amax(D) if local else D[-1][-1]
    y, x = np.unravel_index(D.argmax(), D.shape) if local else (len(s), len(t))

    s_align, t_align = "", ""
    s_matches, t_matches = [], []

    while y != 0 or x != 0:

        current = D[y][x]

        if local and current == 0:  # The end of the best local alignment

            break

        elif current == D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]): # Match (D)

            x -= 1
            y -= 1

            s_align = s[y] + s_align
            t_align = t[x] + t_align

            s_matches.append(y)
            t_matches.append(x)

        elif current == D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"):  # Matching a gap in s with a character in t (L)

            x -= 1
            s_align = "_" + s_align
            t_align = t[x] + t_align

        elif current == D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_"):  # Matching a gap in t with a character in s (U)

            y -= 1
            s_align = s[y] + s_align
            t_align = "_" + t_align

        else:

            raise ValueError("Something's fucked!")

    return score, np.array(rev(s_matches)), np.array(rev(t_matches)), s_align, t_align

def extend(alphabet, scoring_matrix, s, t, seed, threshold):

    x, y, length = seed
    up, down = 1, 1

    while up or down:

        if x == 0 or y == 0 or cost(alphabet, scoring_matrix, t[x - 1], s[y - 1]) < threshold:

            up = 0

        else:

            x -= 1
            y -= 1
            length += 1

        if x + length == len(t) or y + length == len(s) or cost(alphabet, scoring_matrix, t[x + 1], s[y + 1]) < threshold:

            down = 0

        else:

            length += 1

    return x, y, length

def get_threshold(alphabet, scoring_matrix, s, t):

    char_freq = {i: (s + t).count(i) / len(s + t) for i in alphabet}

    pairs = [(i, j) for i in alphabet for j in alphabet]

    pair_weights = [char_freq[x] + char_freq[y] for (x, y) in pairs]

    return np.average([cost(alphabet, scoring_matrix, x, y) for (x, y) in pairs], weights=pair_weights)

def align_score(alphabet, scoring_matrix, s_align, t_align):

    assert(len(s_align) == len(t_align))

    score = 0

    for i in range(len(s_align)):

        score += cost(alphabet, scoring_matrix, s_align[i], t_align[i])

    return score

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def rev(l):

    return l[::-1]

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