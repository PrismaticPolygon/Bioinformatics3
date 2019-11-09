import numpy as np


def c(c1, c2):

    global ALPHABET, SCORING_MATRIX

    i = ALPHABET.index(c1)
    j = ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]


def dynprog(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    backtrack_matrix = np.full((len(s) + 1, len(t) + 1), b" ", dtype="string_")
    backtrack_scoring_matrix = np.zeros((len(s) + 1, len(t) + 1), dtype="int8")

    for i in range(len(s) + 1): # y

        for j in range(len(t) + 1): # x

            if i == 0 and j == 0:   # Top left

                backtrack_scoring_matrix[0, 0] = 0
                backtrack_matrix[0, 0] = "E"

            elif i == 0 and j > 0:  # First row

                backtrack_scoring_matrix[0, j] = backtrack_scoring_matrix[0, j - 1] + c(t[j - 1], "_")
                backtrack_matrix[0, j] = "L"

            elif j == 0 and i > 0:  # First column

                backtrack_scoring_matrix[i, 0] = backtrack_scoring_matrix[i - 1, 0] + c(s[i - 1], "_")
                backtrack_matrix[i, 0] = "U"

            else:

                match = backtrack_scoring_matrix[i - 1][j - 1] + c(s[i - 1], t[j - 1])
                delete = backtrack_scoring_matrix[i][j - 1] + c(s[i - 1], "_")
                insert = backtrack_scoring_matrix[i - 1][j] + c(t[j - 1], "_")

                direction, score = max_score(match, delete, insert)

                backtrack_matrix[i, j] = direction
                backtrack_scoring_matrix[i, j] = score

    best_score = backtrack_scoring_matrix[-1][-1]
    seq1, seq2 = track_back(backtrack_matrix, s, t)

    matches = []
    mismatches = []

    for i in range(len(seq1)):

        if seq1[i] == seq2[i]:

            matches.append(i)

        else:

            mismatches.append(i)

    print(matches, mismatches)

    return best_score, matches, mismatches


def track_back(backtrack_matrix, s, t):

    i, j = len(s), len(t)
    seq1, seq2 = "", ""

    current = backtrack_matrix[-1][-1]

    while current != b'E':

        if current == b'D':

            j -= 1
            i -= 1
            seq1 = s[i] + seq1
            seq2 = t[j] + seq2

        elif current == b'L':

            i -= 1
            seq1 = s[i] + seq1
            seq2 = '-' + seq2

        elif current == b'U':

            j -= 1
            seq1 = '-' + seq1
            seq2 = t[j] + seq2

        current = backtrack_matrix[j][i]

    return seq1, seq2


def max_score(match, insert, delete):

    if match >= insert and match >= delete:

        return "D", match

    elif insert >= delete:

        return "U", insert

    return "L", delete