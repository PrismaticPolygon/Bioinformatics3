import numpy as np
from dynprog import dynprog_align

# That might fuck up the costs. Oh well.
# It should be all in one, of course.
# Though that will be gross....
# For the sake of code duplication, I'm going to do that now. 

def dynproglin(alphabet, scoring_matrix, s, t):

    # print(s, t)

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    z = ""
    w = ""

    # We care about the score, though.
    # But we can calculate that later.

    if len(s) == 0:

        for i in range(len(t)):

            z += "_"
            w += t[i]

    elif len(t) == 0:

        for i in range(len(s)):

            z += s[i]
            w += "_"

    elif len(s) == 1 or len(t) == 1:

        z, w = dynprog_align(alphabet, scoring_matrix, s, t)

    else:

        s_mid = int(len(s) / 2)
        score_l = build_score_matrix(s[:s_mid], t)

        # print("Score_l: ", score_l)

        score_r = build_score_matrix(rev(s[s_mid:]), rev(t))

        # print("Score_r", score_r)

        t_mid = np.argmax(score_l + rev(score_r))

        # print("t_mid", t_mid)

        z_l, w_l = dynproglin(alphabet, scoring_matrix, s[:s_mid], t[:t_mid])
        z_r, w_r = dynproglin(alphabet, scoring_matrix, s[s_mid:], t[t_mid:])

        z = z_l + z_r
        w = w_l + w_r

    k = score(z, w)

    # Different, incorrect, score.
    # And that'll be HARD to pin down.

    # Close but not perfect. But that's to be expected: we're rounding
    # up midpoints and whatnot.
    # Okay. Now I need to consolidate this.
    # Or perhaps I should test on an odd-length string

    # Shame that it is, genuinely, an inferior score.

    # print(z)
    # print(w)
    # print(k)

    return z, w

def score(s_align, t_align):

    k = 0

    for i in range(len(s_align)):

        k += cost(s_align[i], t_align[i])

    return k

def rev(l):

    return l[::-1]

def build_score_matrix(s, t):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    D = np.zeros((2, size_x), dtype="int8")

    for y in range(size_y):  # y

        for x in range(size_x):  # x

            if y == 0 and x > 0:  # First row (the cost of matching t with all gaps)

                D[0, x] = D[0, x - 1] + cost(t[x - 1], "_")

            elif x == 0 and y > 0:  # First column (the cost of matching s with all gaps)

                D[1, 0] = D[0, 0] + cost(s[y - 1], "_")

            elif y != 0 and x != 0:

                D[1, x] = max(
                    match(D, 0, x, s, t),  # The cost of matching two characters
                    insert_gap_into_s(D, 0, x, t),  # The cost of matching a gap in s with a character in t
                    insert_gap_into_t(D, 0, x, s)  # The cost of matching a gap in t with a character in s
                )

        D[0, :] = D[1, :]

    return D[1]

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