import numpy as np
from dynprog import dynprog


def dynproglin(alphabet, scoring_matrix, x, y):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    z = ""
    w = ""

    if len(x) == 0:

        for i in range(1, len(y)):

            z += "_"
            w += y[i]

    elif len(y) == 0:

        for i in range(1, len(x)):

            z += x[i]
            w += "_"

    elif len(x) == 1 or len(y) == 1:

        z, w = dynprog(alphabet, scoring_matrix, x, y)

    else:

        print(nw_score(x, y))

        x_mid = int(len(x) / 2)

        print(x[:x_mid], y)

        score_l = nw_score(x[:x_mid], y)

        print(score_l)

        score_r = nw_score(x[x_mid:][::-1], y[::-1])

        print(score_r)

        y_mid = np.amax(score_l + score_r[::-1])

        print(y_mid)

        # Okay. Now we get recursive.

        z_1, z_2 = dynproglin(alphabet, scoring_matrix, x[:x_mid], y[:y_mid])
        z_3, z_4 = dynproglin(alphabet, scoring_matrix, x[x_mid:], y[y_mid:])

        z = z_1 + z_3
        w = z_2 + z_4

    print(z, w)

    return z, w


def c(c1, c2):

    global ALPHABET, SCORING_MATRIX

    i = ALPHABET.index(c1)
    j = ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]


def nw_score(s, t):

    backtrack_scoring_matrix = np.zeros((2, len(t) + 1), dtype="int8")

    for i in range(len(s) + 1):  # y

        for j in range(len(t) + 1):  # x

            if i == 0 and j > 0:  # First row

                backtrack_scoring_matrix[0, j] = backtrack_scoring_matrix[0, j - 1] + c(t[j - 1], "_")

            elif j == 0 and i > 0:  # First column

                backtrack_scoring_matrix[1, 0] = backtrack_scoring_matrix[0, 0] + c(s[i - 1], "_")

            elif i != 0 and j != 0:

                match = backtrack_scoring_matrix[0][j - 1] + c(s[i - 1], t[j - 1])
                insert = backtrack_scoring_matrix[0][j] + c(t[j - 1], "_")
                delete = backtrack_scoring_matrix[1][j - 1] + c(s[i - 1], "_")

                backtrack_scoring_matrix[1, j] = max(match, delete, insert)

        backtrack_scoring_matrix[0, :] = backtrack_scoring_matrix[1, :]

    return backtrack_scoring_matrix[1] # The last_line of the NW score matrix.


def needleman_wunsch(x, y):

    return x, y


def hirschberg(x, y):

    z = ""
    w = ""

    if len(x) == 0:

        for i in range(1, len(y)):

            z += "_"
            w += y[i]

    elif len(y) == 0:

        for i in range(1, len(x)):

            z += x[i]
            w += "_"

    elif len(x) == 1 or len(y) == 1:

        z, w = needleman_wunsch(x, y)

    else:

        print(nw_score(x, y))

        x_mid = int(len(x) / 2)

        print(x[:x_mid], y)

        score_l = nw_score(x[:x_mid], y)

        print(score_l)

        score_r = nw_score(x[x_mid:][::-1], y[::-1])

        print(score_r)

        y_mid = np.amax(score_l + score_r[::-1])

        # Okay. Now we get recursive.

        # z_1, z_2 = hirschberg(x[:x_mid], y[:y_mid])
        # z_3, z_4 = hirschberg(x[x_mid:], y[y_mid:])
        #
        # z, w = z_1 + z_3, z_2 + z_4

    return z, w



