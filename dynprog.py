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

                match_score = match(D, y, x, s, t)
                insert_gap_into_s_cost = insert_gap_into_s(D, y, x, t)
                insert_gap_into_t_cost = insert_gap_into_t(D, y, x, s)

                # print((y, x))
                # print("Match: ", match_score)
                # print("Insert into s: ", insert_gap_into_s_cost)
                # print("Insert into t:", insert_gap_into_t_cost)

                D[y, x] = max(
                    match(D, y, x, s, t),           # The cost of matching two characters
                    insert_gap_into_s(D, y, x, t),  # The cost of matching a gap in s with a character in t
                    insert_gap_into_t(D, y, x, s)   # The cost of matching a gap in t with a character in s
                )

            print(D)

            # So it goes wrong on the last bit of the third row. We choose a value which is too high.
            # Logically, we must get a 4 from somewhere.
            # But there are no 4s!

    print(D)

    # This isn't actually right.

    score = D[-1][-1]
    print(score)
    s_align, t_align, s_matches, t_matches = traceback(D, s, t)

    print(s_align)
    print(t_align)
    print(s_matches)
    print(t_matches)

    return score, s_matches, t_matches

    #
    # print("{} -> {}".format(s, s_align))
    # print("{} -> {}".format(t, t_align))
    # print(score)

# Going up insert into x-axis string

def insert_gap_into_s(D, y, x, t):  # Conceptually L

    return D[y - 1][x] + cost(t[x - 1], "_")


def insert_gap_into_t(D, y, x, s):  # Conceptually U

    return D[y][x - 1] + cost(s[y - 1], "_")


def match(D, y, x, s, t):   # Conceptually D

    return D[y - 1][x - 1] + cost(s[y - 1], t[x - 1])


def traceback(D, s, t):

    y, x = len(s), len(t)

    # print(y, x, D.shape)
    print(s, t)

    s_align, t_align = "", ""
    s_matches = []
    t_matches = []

    while y != 0 or x != 0:

        current = D[y][x]

        # print("({}, {})".format(y, x))

        if current == match(D, y, x, s, t): #D

            # print("Match")

            x -= 1
            y -= 1
            s_align = s[y] + s_align
            t_align = t[x] + t_align

            s_matches.append(y)
            t_matches.append(x)

        elif current == insert_gap_into_s(D, y, x, t):  # L

            # print("Insert gap into t")

            y -= 1
            s_align = s[y] + s_align
            t_align = "_" + t_align

        elif current == insert_gap_into_t(D, y, x, s):  # U

            # print("Insert gap into s")

            x -= 1
            s_align = "_" + s_align
            t_align = t[x] + t_align

        else:

            raise ValueError("Something's fucked!")

        # print(s_align, t_align)

    return s_align, t_align, s_matches[::-1], t_matches[::-1]


def cost(c1, c2):

    global ALPHABET, SCORING_MATRIX

    i = ALPHABET.index(c1)
    j = ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]