import numpy as np
from itertools import product

# http://biorecipes.com/DynProgBasic/code.html
def dynprog(alphabet, scoring_matrix, s, t):

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)

    D = smith_waterman(alphabet, scoring_matrix, s, t)

    return traceback(alphabet, scoring_matrix, s, t, D)

# https://en.wikipedia.org/wiki/Hirschberg's_algorithm
def dynproglin(alphabet, scoring_matrix, s, t):

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)

    def recurse(s, t):

        if len(s) == 0 or len(t) == 0:

            s_align, t_align = "", ""

        elif len(s) == 1 or len(t) == 1:

            result = dynprog(alphabet, scoring_matrix, s, t)
            s_align = result[3]
            t_align = result[4]

        else:

            s_mid = int(len(s) / 2)
            score_l = hirschberg(alphabet, scoring_matrix, s[:s_mid], t)
            score_r = hirschberg(alphabet, scoring_matrix, rev(s[s_mid:]), rev(t))
            t_mid = np.argmax(score_l + rev(score_r))

            left = recurse(s[:s_mid], t[:t_mid])
            right = recurse(s[s_mid:], t[t_mid:])

            s_align = left[0] + right[0]
            t_align = left[1] + right[1]

        return s_align, t_align

    # Save node i - j as part of the solution.

    s_align, t_align = recurse(s, t)

    # Should we be able to calculate it here?
    # I feel like yes.

    # Hm. I'm trying to avoid having those auxiliary alignment methods.
    # But I don't understand how well enough how it works, lmao.

    score = align_score(alphabet, scoring_matrix, s_align, t_align)
    s_matches, t_matches = get_alignment_indices(s_align, t_align)

    return score, s_matches, t_matches, s_align, t_align

# file:///C:/Users/user/Documents/MEGA/University/Year%203/CCS/Bioinformatics/HeuristicAlign.pdf
def heuralign(alphabet, scoring_matrix, s, t):

    # s (database), |s| >> |t| is the y-axis
    # t (query), is the x-axis

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)
    ktup = 3
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

def smith_waterman(alphabet, scoring_matrix, s, t):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    D = np.zeros((size_y, size_x), dtype="int8")

    for y in range(1, size_y):  # y

        for x in range(1, size_x):

            D[y, x] = max(
                0,  # For local alignment
                D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]), # Cost of matching two characters
                D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"),  # The cost of matching a gap in s with a character in t
                D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_"),   # The cost of matching a gap in t with a character in s
            )

    return D

def hirschberg(alphabet, scoring_matrix, s, t):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    D = np.zeros((2, size_x), dtype="int8")

    for y in range(1, size_y):  # y

        for x in range(1, size_x):

            D[1, x] = max(
                0,  # For local alignment
                D[0][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),   # The cost of matching two characters
                D[1][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"),    # The cost of matching a gap in s with a character in t
                D[0][x] + cost(alphabet, scoring_matrix, s[y - 1], "_") # The cost of matching a gap in t with a character in s
            )

        D[0] = D[1].copy()

    return D[1]

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

def traceback(alphabet, scoring_matrix, s, t, D):

    score = np.amax(D)
    y, x = np.unravel_index(D.argmax(), D.shape)

    s_align, t_align = "", ""
    s_matches, t_matches = [], []

    while y != 0 or x != 0:

        current = D[y][x]

        if current == 0:  # The end of the best local alignment

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

    return score, s_matches[::-1], t_matches[::-1], s_align, t_align

def get_alignment_indices(s_align, t_align):

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

    return s_matches, t_matches


def align_score(alphabet, scoring_matrix, s_align, t_align):

    score = 0

    for i in range(len(s_align)):

        score += cost(alphabet, scoring_matrix, s_align[i], t_align[i])

    return score

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def rev(l):

    return l[::-1]