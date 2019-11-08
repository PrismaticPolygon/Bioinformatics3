import numpy as np

# string, 2D matrix, string, string.

# Ah. I'm passed in the scoring_matrix.
# So the cost of a gap is the last element of the scoring matrix.

# Save the optimal scores for the solution of every subproblem instead of recalculating them
# For two strings, s of length m and t of length n, D[i, j] is the best score of aligning the two substrings
# s[1..j] and t[1...i]
# For global alignment... the best score is the last value in the table.

def c(c1, c2):

    global ALPHABET, SCORING_MATRIX

    i = ALPHABET.index(c1)
    j = ALPHABET.index(c2)

    return SCORING_MATRIX[i, j]


def dynprog(alphabet, scoring_matrix, s, t):

    global ALPHABET, SCORING_MATRIX

    ALPHABET = alphabet + "_"
    SCORING_MATRIX = np.array(scoring_matrix)

    backtrack_scoring_matrix = build_backtrack_scoring_matrix(s, t)
    backtrack_matrix = build_backtrack_matrix(len(s), len(t))

    for i in range(1, len(s) + 1):

        for j in range(1, len(t) + 1):

            match = backtrack_scoring_matrix[i - 1][j - 1] + c(s[i - 1], t[j - 1])
            gap_s = backtrack_scoring_matrix[i][j - 1] + c(s[i - 1], "_")
            gap_t = backtrack_scoring_matrix[i - 1][j] + c(t[j - 1], "_")

            direction, score = max_value(match, gap_s, gap_t)

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

    print(seq1, seq2)

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


def max_value(score_d, score_u, score_l):

    if score_d >= score_u and score_d >= score_l:

        return "D", score_d

    elif score_u >= score_l:

        return "U", score_u

    return "L", score_l


def build_backtrack_scoring_matrix(s, t):

    i = len(s) + 1
    j = len(t) + 1
    backtrack_scoring_matrix = np.zeros((i, j), dtype="int8")

    for k in range(1, max(i, j)):

        if k < j:

            backtrack_scoring_matrix[0, k] = backtrack_scoring_matrix[0, k - 1] + c(t[k - 1], "_")

        if k < i:

            backtrack_scoring_matrix[k, 0] = backtrack_scoring_matrix[k - 1, 0] + c(s[k - 1], "_")

    return backtrack_scoring_matrix


def build_backtrack_matrix(i, j):

    backtrack_matrix = np.full((i + 1, j + 1), b" ", dtype="string_")

    backtrack_matrix[:, 0] = "U"
    backtrack_matrix[0] = "L"
    backtrack_matrix[0, 0] = "E"

    return backtrack_matrix