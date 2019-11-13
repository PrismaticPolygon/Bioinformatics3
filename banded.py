from itertools import product
import numpy as np

# Real-life matches often contain long strings with gap-less matches.
# Heuristics try to find significant gap-less matches and extend them.

def build_score_matrix_banded(s, t, k):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    D = np.zeros((size_y, size_x), dtype="int8")
    width = int(k / 2)

    for y in range(size_y):

        for x in range(y - width, y + width + 1):

            if 0 <= x < size_x and y > 0:

                D[y, x] = max(
                    0,  # For local alignment
                    match(D, y, x, s, y, t),  # The cost of matching two characters
                    insert_gap_into_s(D, y, x, t),  # The cost of matching a gap in s with a character in t
                    insert_gap_into_t(D, y, x, s, y)  # The cost of matching a gap in t with a character in s
                )

    return D


def match(D, D_y, x, s, s_y, t):   # Conceptually D

    return D[D_y - 1][x - 1] + cost(s[s_y - 1], t[x - 1])

def insert_gap_into_s(D, D_y, x, t):  # s is the y-axis string: conceptually L

    return D[D_y][x - 1] + cost(t[x - 1], "_")

def insert_gap_into_t(D, D_y, x, s, s_y):  # t is the x-axis string: conceptually U

    return D[D_y - 1][x] + cost(s[s_y - 1], "_")



def banded(s, t, matches):

    diagonals = dict()
    k = 9

    # I'm not sure how. Actually, that's untrue. It'll be the nearest (under) multiple of k.
    # So then I guess part of the heuristic is identifying the diagonals! I'll do the most simple kind of banded for now:

    for match in matches:

        x, y, length = match

        i = k * ((y - x) // k)

        if i not in diagonals:

            diagonals[i] = [match]

        else:

            diagonals[i].append(match)

    for key, value in diagonals.items():

        print(value)

        start_x = min(value, key=lambda m:m[0])[0]
        start_y = min(value, key=lambda m:m[1])[1]
        end_x_match = max(value, key=lambda m:m[0] + m[2])
        end_y_match = max(value, key=lambda m:m[1] + m[2])

        end_x = end_x_match[0] + end_x_match[2]
        end_y = end_y_match[1] + end_y_match[2]

        print(start_x, start_y, end_x, end_y)

        print(s, t)

        s = s[start_y:end_y]
        t = t[start_x:end_x]

        D = build_score_matrix_banded(s, t, k)

        score, s_align, t_align, s_matches, t_matches = traceback(D, s, t)

        print(score, t_align, s_align, s_matches, t_matches)


def score_match(s, t, match):

    x, y, length = match

    score = 0

    for i in range(length):

        score += cost(s[y + i], t[x + i])

    return score

def extend(s, t, match):

    x, y, length = match

    while x > 0 and y > 0:

        if cost(t[x - 1], s[y - 1]) > 0:

            x -= 1
            y -= 1
            length += 1

        else:

            break

    while x + length + 1 < len(t) and y + length + 1 < len(s):

        if cost(t[x + 1], s[y + 1]) > 0:

            length += 1

        else:

            break

    return x, y, length

def heuralign(alphabet, scoring_matrix, s, t):

    # s (database), |s| >> |t| is the y-axis
    # t (query), is the x-axis

    global SCORING_MATRIX, ALPAHBET

    ktup = 3
    SCORING_MATRIX = scoring_matrix
    lookup = {"".join(i): [] for i in product(alphabet, repeat=ktup)}
    ALPAHBET = alphabet + "_"

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]
        lookup[word].append(y)

    matches = []

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        for y in lookup[word]:

            extended = extend(s, t, (x, y, ktup))

            matches.append(extended)

    banded(s, t, matches)

def cost(c1, c2):

    global ALPHABET

    return scoring_matrix[ALPHABET.index(c1), ALPHABET.index(c2)]

def traceback(D, s, t):

    score = np.amax(D)
    y, x = np.unravel_index(D.argmax(), D.shape)

    s_align, t_align = "", ""
    s_matches = []
    t_matches = []

    while y != 0 or x != 0:

        current = D[y][x]

        if current == 0:  # The end of the best local alignment

            return score, s_align, t_align, s_matches[::-1], t_matches[::-1]

        elif current == match(D, y, x, s, y, t): # D

            x -= 1
            y -= 1
            s_align = s[y] + s_align
            t_align = t[x] + t_align

            s_matches.append(y)
            t_matches.append(x)

        elif current == insert_gap_into_s(D, y, x, t):  # L

            x -= 1
            s_align = "_" + s_align
            t_align = t[x] + t_align


        elif current == insert_gap_into_t(D, y, x, s, y):  # U

            y -= 1
            s_align = s[y] + s_align
            t_align = "_" + t_align

        else:

            raise ValueError("Something's fucked!")

    return score, s_align, t_align, s_matches[::-1], t_matches[::-1]


alphabet = "ABC"

scoring_matrix = np.array([
    [1, -1, -2, -1],
    [-1, 2, -4, -1],
    [-2, -4, 3, -2],
    [-1, -1, -2, 0]
])


fasta("ABCCCABABACABCABCABCBAABABCCCAAACBCBCBABCABCBABBBCABCA", "AAACCBACBAC", 3)