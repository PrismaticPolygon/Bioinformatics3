import numpy as np


def dynprog(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    D = np.zeros((size_y, size_x), dtype="int8") # Implicitly fill top-left corner as 0.

    for y in range(size_y): # y

        for x in range(size_x): # x

            if y == 0 and x > 0:  # First row (the cost of matching t with all gaps)

                D[0, x] = D[0, x - 1] + cost(t[x - 1], "_")

            elif x == 0 and y > 0:  # First column (the cost of matching s with all gaps)

                D[y, 0] = D[y - 1, 0] + cost(s[y - 1], "_")

            elif y != 0 and x != 0:

                D[y, x] = max(
                    match(D, y, x, s, t),           # The cost of matching two characters
                    insert_gap_into_s(D, y, x, t),  # The cost of matching a gap in s with a character in t
                    insert_gap_into_t(D, y, x, s)   # The cost of matching a gap in t with a character in s
                )

    score = D[-1][-1]
    s_align, t_align, s_matches, t_matches = traceback(D, s, t)

    print(s_align)
    print(t_align)

    return score, s_matches, t_matches

# s is the y-axis string
def insert_gap_into_s(D, y, x, t):  # Conceptually L

    return D[y][x - 1] + cost(t[x - 1], "_")

# t is the x-axis string
def insert_gap_into_t(D, y, x, s):  # Conceptually U

    return D[y - 1][x] + cost(s[y - 1], "_")


def match(D, y, x, s, t):   # Conceptually D

    return D[y - 1][x - 1] + cost(s[y - 1], t[x - 1])


def traceback(D, s, t):

    y, x = len(s), len(t)

    s_align, t_align = "", ""
    s_matches = []
    t_matches = []

    while y != 0 or x != 0:

        current = D[y][x]

        if current == match(D, y, x, s, t): # D

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


        elif current == insert_gap_into_t(D, y, x, s):  # U

            y -= 1
            s_align = s[y] + s_align
            t_align = "_" + t_align

        else:

            raise ValueError("Something's fucked!")

    return s_align, t_align, s_matches[::-1], t_matches[::-1]


def cost(c1, c2):

    global ALPHABET, SCORING_MATRIX

    i = ALPHABET.index(c1)
    j = ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]