import numpy as np

def dynprog(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    D = build_score_matrix(s, t)

    score = D[-1][-1]
    s_align, t_align, s_matches, t_matches = traceback(D, s, t)

    print(s_align)
    print(t_align)

    return score, s_matches, t_matches


def dynproglin_recursive(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    s_align, t_align = "", ""

    if len(s) == 0:

        for i in range(len(t)):
            s_align += "_"
            t_align += t[i]

    elif len(t) == 0:

        for i in range(len(s)):
            s_align += s[i]
            t_align += "_"

    elif len(s) == 1 or len(t) == 1:

        s_align, t_align = dynprog_align(alphabet, scoring_matrix, s, t)

    else:

        s_mid = int(len(s) / 2)
        score_l = build_score_matrix(s[:s_mid], t, sublinear=True)

        # print("Score_l: ", score_l)

        score_r = build_score_matrix(rev(s[s_mid:]), rev(t), sublinear=True)

        # print("Score_r", score_r)

        t_mid = np.argmax(score_l + rev(score_r))

        # print("t_mid", t_mid)

        z_l, w_l = dynproglin_recursive(alphabet, scoring_matrix, s[:s_mid], t[:t_mid])
        z_r, w_r = dynproglin_recursive(alphabet, scoring_matrix, s[s_mid:], t[t_mid:])

        s_align = z_l + z_r
        t_align = w_l + w_r

    return s_align, t_align



def dynproglin(alphabet, scoring_matrix, s, t):

    s_align, t_align = dynproglin_recursive(alphabet, scoring_matrix, s, t)
    score = align_score(s_align, t_align)
    s_matches, t_matches = [], []

    print(s_align)
    print(t_align)

    for i in range(len(s_align)):

        if s_align[i] != "_":

            t_matches.append(i)

        if t_align[i] != "_":

            s_matches.append(i)

    return score, s_matches, t_matches


def build_score_matrix(s, t, sublinear=False):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis
    shape = (2, size_x) if sublinear else (size_y, size_x)
    D = np.zeros(shape, dtype="int8")

    for y in range(size_y):  # y

        for x in range(size_x):  # x

            if y == 0 and x > 0:  # First row (the cost of matching t with all gaps)

                D[0, x] = D[0, x - 1] + cost(t[x - 1], "_")

            elif x == 0 and y > 0:  # First column (the cost of matching s with all gaps)

                new_cell = (1, 0) if sublinear else (y, 0)
                old_cell = (0, 0) if sublinear else (y - 1, 0)
                D[new_cell] = D[old_cell] + cost(s[y - 1], "_")

            elif y != 0 and x != 0:

                if sublinear:

                    D[1, x] = max(
                        match(D, 0, x, s, t),  # The cost of matching two characters
                        insert_gap_into_s(D, 0, x, t),  # The cost of matching a gap in s with a character in t
                        insert_gap_into_t(D, 0, x, s)  # The cost of matching a gap in t with a character in s
                    )

                else:

                    D[y, x] = max(
                        match(D, y, x, s, t),  # The cost of matching two characters
                        insert_gap_into_s(D, y, x, t),  # The cost of matching a gap in s with a character in t
                        insert_gap_into_t(D, y, x, s)  # The cost of matching a gap in t with a character in s
                    )

        if sublinear:

            D[0, :] = D[1, :]

    return D[1] if sublinear else D


def dynprog_align(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    D = build_score_matrix(s, t)

    score = D[-1][-1]
    s_align, t_align, s_matches, t_matches = traceback(D, s, t)

    return s_align, t_align


def insert_gap_into_s(D, y, x, t):  # s is the y-axis string: conceptually L

    return D[y][x - 1] + cost(t[x - 1], "_")


def insert_gap_into_t(D, y, x, s):  # t is the x-axis string: conceptually U

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

    i, j = ALPHABET.index(c1), ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]

def align_score(s_align, t_align):

    assert len(s_align) == len(t_align)

    score = 0

    for i in range(len(s_align)):

        score += cost(s_align[i], t_align[i])

    return score

def rev(l):

    return l[::-1]