import numpy as np

# http://biorecipes.com/DynProgBasic/code.html
def dynprog(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    D = build_score_matrix(s, t)

    score, s_align, t_align, s_matches, t_matches = traceback(D, s, t)

    return score, s_matches, t_matches, s_align, t_align

# Okay, we're in.
# I need to figure out why dynproglin isn't working. Should be easy, right? No.
# I suspect it's in the nuances of the algorithm itself.
# If either equals 0, then there is no best local alignment.

#https://en.wikipedia.org/wiki/Hirschberg's_algorithm
def dynproglin(alphabet, scoring_matrix, s, t):

    # Snd there should never be any gaps matched with others gaps: it's always less than 0.
    # That's not the issue, mind.
    # If the length of either is 1.. that's fine. Provided the letter is in the other string, we can just return that.

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    def recurse(s, t):

        s_align, t_align = "", ""

        if len(s) == 0 or len(t) == 0:

            s_align = ""
            t_align = ""

        elif len(s) == 1 or len(t) == 1:

            result = dynprog(alphabet, scoring_matrix, s, t)
            s_align = result[3]
            t_align = result[4]

        else:

            s_mid = int(len(s) / 2)
            score_l = build_score_matrix(s[:s_mid], t, sublinear=True)
            score_r = build_score_matrix(rev(s[s_mid:]), rev(t), sublinear=True)
            t_mid = np.argmax(score_l + rev(score_r))

            z_l, w_l = recurse(s[:s_mid], t[:t_mid])
            z_r, w_r = recurse(s[s_mid:], t[t_mid:])

            s_align = z_l + z_r
            t_align = w_l + w_r

        return s_align, t_align

    s_align, t_align = recurse(s, t)

    # So it's BACA, BACA. This doesn't exist.
    # If the scoring matrix is the same... let's just double-check again that it is. For my peace of mind, at least.
    # maybe that's part of the complexity?
    # Ah. It's a different matrix. On the second-last row, mind, which indicates to me it's an off-by-one.
    # The alignment, score, and matches are different. Of course it's all of them. But more interestingly, we're appended characters that shouldn't eist.

    print(s_align, t_align)

    score = align_score(s_align, t_align)
    s_matches, t_matches = get_alignment_indices(s_align, t_align)

    # No... it's more than that.
    # We're actually a row short.

    return score, s_matches, t_matches, s_align, t_align

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

def build_score_matrix(s, t, sublinear=False):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    shape = (2, size_x) if sublinear else (size_y, size_x)
    D = np.zeros(shape, dtype="int8")

    for y in range(1, size_y):  # y

        D_y = 1 if sublinear else y

        for x in range(1, size_x):

            D[D_y, x] = max(
                0,  # For local alignment
                match(D, D_y, x, s, y, t),  # The cost of matching two characters
                insert_gap_into_s(D, D_y, x, t),  # The cost of matching a gap in s with a character in t
                insert_gap_into_t(D, D_y, x, s, y)  # The cost of matching a gap in t with a character in s
            )

        print(D[0])

        if sublinear and y < size_y - 1: # Copy the 1st row onto the second row unless it's the final iteration

            D[0] = D[1].copy()
            D[1] = 0

    return D[1] if sublinear else D

def match(D, D_y, x, s, s_y, t):   # Conceptually D

    return D[D_y - 1][x - 1] + cost(s[s_y - 1], t[x - 1])

def insert_gap_into_s(D, D_y, x, t):  # s is the y-axis string: conceptually L

    return D[D_y][x - 1] + cost(t[x - 1], "_")

def insert_gap_into_t(D, D_y, x, s, s_y):  # t is the x-axis string: conceptually U

    return D[D_y - 1][x] + cost(s[s_y - 1], "_")

def traceback(D, s, t):

    score = np.amax(D)
    y, x = np.unravel_index(D.argmax(), D.shape)

    # We don't need the traceback, right?

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

def cost(c1, c2):

    global ALPHABET, SCORING_MATRIX

    return SCORING_MATRIX[ALPHABET.index(c1), ALPHABET.index(c2)]

def align_score(s_align, t_align):

    score = 0

    for i in range(len(s_align)):

        score += cost(s_align[i], t_align[i])

    return score

def rev(l):

    return l[::-1]

#
# ALPHABET = "ABC_"
# SCORING_MATRIX = np.array([
#     [1, -1, -2, -1],
#     [-1, 2, -4, -1],
#     [-2, -4, 3, -2],
#     [-1, -1, -2, 0]
# ])
# s = "AABBAACA"
# t = "CBACCCBA"
#
# print(build_score_matrix(s, t))
#
# print(build_score_matrix(s, t, sublinear=True))