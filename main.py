import numpy as np
from itertools import product

# http://biorecipes.com/DynProgBasic/code.html
def dynprog(alphabet, scoring_matrix, s, t):

    alphabet = alphabet if alphabet[-1] == "_" else alphabet + "_"
    scoring_matrix = np.array(scoring_matrix)

    D, max_i = score_matrix(alphabet, scoring_matrix, s, t)

    return traceback(alphabet, scoring_matrix, s, t, D)

# Smith Waterman is used for local
# Needleman-Wunsch is used for global

def dynproglin(alphabet, scoring_matrix, s, t):

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)

    # Get max index from forwards pass
    # Get min index from backwards pass.
    # He literally just reverses, it, mind you.
    # Is that the right thing to do? I'm not sure.

    # print(s, t)

    _, (end_y, end_x) = score_matrix(alphabet, scoring_matrix, s, t, sublinear=True)

    # print(end_y, end_x)

    s = s[:end_y]
    t = t[:end_x]
    #
    # print(s, t)

    s = rev(s)
    t = rev(t)

    # print(s, t)

    _, (start_y, start_x) = score_matrix(alphabet, scoring_matrix, s, t, sublinear=True)

    # print(start_y, start_x)

    s = rev(s[:start_y])
    t = rev(t[:start_x])

    # print(s, t)

    s_align, t_align = hirschberg(alphabet, scoring_matrix, s, t)
    score = align_score(alphabet, scoring_matrix, s_align, t_align)

    s_matches, t_matches = get_alignment_indices(s_align, t_align)

    s_matches += end_y - start_y
    t_matches += end_x - start_x

    return score, s_matches, t_matches, s_align, t_align

# file:///C:/Users/user/Documents/MEGA/University/Year%203/CCS/Bioinformatics/HeuristicAlign.pdf
def heuralign(alphabet, scoring_matrix, s, t):

    #linspace return evenly spaced numbers over a specified interval
    # numpy ufuncs have a reduceat function.

    # s (database), |s| >> |t| is the y-axis
    # t (query), is the x-axis

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)
    ktup = 4
    k = 9

    lookup = build_lookup(alphabet, s, ktup)    # Build lookup table of database string s
    matches = build_matches(alphabet, scoring_matrix, s, t, ktup, lookup)   # Build list of extended alignments from query string t
    diagonals = build_diagonals(k, matches) # Sort into diagonals.

    for key, matches in diagonals.items():

        start_x, start_y = len(t) - 1, len(s) - 1
        end_x, end_y = 0, 0

        for x, y, length in matches:

            if x < start_x:

                start_x = x

            if x + length > end_x:

                end_x = x + length

            if y < start_y:

                start_y = y

            if y + length > end_y:

                end_y = y + length

        s = s[start_y:end_y]
        t = t[start_x:end_x]

        D = banded(alphabet, scoring_matrix, s, t, k)

        return traceback(alphabet, scoring_matrix, s, t, D)

# https://en.wikipedia.org/wiki/Hirschberg's_algorithm
def hirschberg(alphabet, scoring_matrix, s, t):

    s_align, t_align = "", ""
    # We can calculate the score in here, it seems. Not very well, mind you: I'd rather not risk breaking anything.

    if len(s) == 0:

        for i in range(len(t)):

            s_align += "_"
            t_align += t[i]

    elif len(t) == 0:

        for i in range(len(s)):

            s_align += s[i]
            t_align += "_"

    elif len(s) == 1 or len(t) == 1:

        _, _, _, s_align, t_align = dynprog(alphabet, scoring_matrix, s, t)

    else:

        s_mid = len(s) // 2
        score_l, _ = score_matrix(alphabet, scoring_matrix, s[:s_mid], t, local=False, sublinear=True)

        # print("Score_l: ", score_l)

        score_r, _ = score_matrix(alphabet, scoring_matrix, rev(s[s_mid:]), rev(t), local=False, sublinear=True)

        # print("Score_r", score_r)

        # Truth value of an array is ambigous

        t_mid = np.argmax(score_l + rev(score_r))

        # print("t_mid", t_mid)

        z_l, w_l = hirschberg(alphabet, scoring_matrix, s[:s_mid], t[:t_mid])
        z_r, w_r = hirschberg(alphabet, scoring_matrix, s[s_mid:], t[t_mid:])

        s_align = z_l + z_r
        t_align = w_l + w_r

    return s_align, t_align


def build_lookup(alphabet, s, ktup):

    lookup = {"".join(i): [] for i in product(alphabet, repeat=ktup)}

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]
        lookup[word].append(y)

    return lookup

def build_matches(alphabet, scoring_matrix, s, t, ktup, lookup):

    matches = []

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        for y in lookup[word]:

            extended = extend(alphabet, scoring_matrix, s, t, x, y, ktup)

            matches.append(extended)

    return matches

def extend(alphabet, scoring_matrix, s, t, x, y, length):

    up, down = True, True

    while up or down:

        if x == 0 or y == 0 or cost(alphabet, scoring_matrix, t[x - 1], s[y - 1]) < 0:

            up = False

        else:

            x -= 1
            y -= 1
            length += 1

        if x + length == len(t) - 1 or y + length == len(s) - 1 or cost(alphabet, scoring_matrix, t[x + 1], s[y + 1]) < 0:

            down = False

        else:

            length += 1

    return x, y, length

def build_diagonals(k, matches):

    diagonals = dict()

    for match in matches:

        x, y, length = match

        i = k * ((y - x) // k)

        if i not in diagonals:

            diagonals[i] = [match]

        else:

            diagonals[i].append(match)

    return diagonals

def score_matrix(alphabet, scoring_matrix, s, t, local=True, sublinear=False):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    shape = (2, size_x) if sublinear else (size_y, size_x)

    D = np.zeros(shape, dtype="int16")

    max_i = None
    max_score = None

    for y in range(size_y):  # y

        for x in range(size_x):

            if not local and y == 0 and x > 0:    # First row: the cost of matching t with all gaps

                D[0, x] = D[0, x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_")

            elif not local and x == 0 and y > 0:  # First column: the cost of matching s with all gaps

                D_y = 1 if sublinear else y

                D[D_y, 0] = D[D_y - 1, 0] + cost(alphabet, scoring_matrix, s[y - 1], "_")

            elif x > 0 and y > 0:

                D_y = 1 if sublinear else y

                # Ah, right.

                D[D_y, x] = max(
                    D[D_y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),  # The cost of matching two characters
                    D[D_y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"), # The cost of matching a gap in s with a character in t
                    D[D_y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_") # The cost of matching a gap in t with a character in s
                )

                D[D_y, x] = 0 if (local and D[D_y, x] < 0) else D[D_y, x]

                if max_i is None or D[D_y, x] > max_score:

                    max_i = y, x
                    max_score = D[D_y, x]

        if y > 0 and sublinear:

            D[0] = D[1].copy()

    return (D[1] if sublinear else D), max_i

def banded(alphabet, scoring_matrix, s, t, k):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    D = np.zeros((size_y, size_x), dtype="int8")
    width = int(k / 2)

    for y in range(size_y):

        for x in range(y - width, y + width + 1):

            if 0 <= x < size_x and y > 0:

                D[y, x] = max(
                    0,  # For local alignment
                    D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]), # Cost of matching two characters
                    D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"), # The cost of matching a gap in s with a character in t
                    D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_") # The cost of matching a gap in t with a character in s
                )

    return D

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

def get_alignment_indices(s_align, t_align):

    assert(len(s_align) == len(t_align))

    s_matches, t_matches = [], []
    s_point, t_point = 0, 0

    for i in range(len(s_align)):

        if s_align[i] != "_" and t_align[i] != "_":

            s_matches.append(s_point)
            t_matches.append(t_point)

            s_point += 1
            t_point += 1

        if s_align[i] != "_" and t_align[i] == "_":

            s_point += 1

        if s_align[i] == "_" and t_align[i] != "_":

            t_point += 1

    return np.array(s_matches), np.array(t_matches)

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