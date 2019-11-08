import numpy as np

# string, 2D matrix, string, string.

# Ah. I'm passed in the scoring_matrix.
# So the cost of a gap is the last element of the scoring matrix.

# Save the optimal scores for the solution of every subproblem instead of recalculating them
# For two strings, s of length m and t of length n, D[i, j] is the best score of aligning the two substrings
# s[1..j] and t[1...i]
# For global alignment... the best score is the last value in the table.

def c(c1, c2, alphabet, scoring_matrix):

    i = alphabet.index(c1)
    j = alphabet.index(c2)

    return scoring_matrix[i, j]


def dynprog(alphabet, scoring_matrix, seq1, seq2):

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)

    backtrack_scoring_matrix = build_backtrack_scoring_matrix(seq1, seq2, alphabet, scoring_matrix)
    backtrack_matrix = build_backtrack_matrix(len(seq1), len(seq2))

    print(scoring_matrix)
    print(backtrack_scoring_matrix)
    print(backtrack_matrix)

    print(scoring_matrix.shape)
    print(backtrack_scoring_matrix.shape)
    print(backtrack_matrix.shape)

    # So my backtrack matric is now too small.
    # Ah. It depends on the sequence character. Logically, then, we pass it
    # through and calculate it that way. Let's try that.
    # Where am I getting this -2 from?
    # Ah. The values of D[i, j] must be explicitly stated. When there are no matrix elements above and left.
    # This corresponds to aligning strings with gaps....

    # Okay. Now we're in business. I can't imagine that this still works, mind.

    # for i in range(1, len(seq1) + 1):
    #
    #     for j in range(1, len(seq2) + 1):
    #
    #         score = max_value(
    #             c(seq1[i - 1], seq2[j - 1], alphabet, scoring_matrix) + backtrack_scoring_matrix[i - 1][j - 1],
    #             backtrack_scoring_matrix[i - 1][j] - 2, backtrack_scoring_matrix[i][j - 1] - 2,
    #             backtrack_matrix, i, j
    #         )
    #
    #         backtrack_scoring_matrix[i, j] = score
    #
    # best_score = backtrack_scoring_matrix[-1][-1]
    # #
    # print(backtrack_matrix)
    # print(backtrack_scoring_matrix)
    # print(best_score)
    #
    # seq1, seq2 = track_back(backtrack_matrix, seq1, seq2)
    #
    # print(seq1, seq2)
    # # print(best_score)
    #
    return 0, 1, 2



def track_back(backtrack_matrix, seq1, seq2):

    # Okay. Now things ARE getting interesting.

    j = len(seq1)
    i = len(seq2)

    fin1 = ''
    fin2 = ''

    # Ah. Basically, a back matrix.

    current = backtrack_matrix[-1][-1]

    while current != "E":

        if current == "D":

            i -= 1
            j -= 1
            fin1 = seq1[j] + fin1
            fin2 = seq2[i] + fin2

        elif current == "L":

            j -= 1
            fin1 = seq1[j] + fin1
            fin2 = '-' + fin2

        elif current == "U":

            i -= 1
            fin1 = '-' + fin1
            fin2 = seq2[i] + fin2

        current = backtrack_matrix[i][j]

    return fin1, fin2


def max_value(score_d, score_u, score_l, backtrack_matrix, i, j):

    if score_d >= score_u and score_d >= score_l:

        backtrack_matrix[i, j] = "D"

        return score_d

    elif score_u >= score_l:

        backtrack_matrix[i, j] = "U"

        return score_u

    backtrack_matrix[i, j] = "L"

    return score_l



# Now, can I actually remember how this works? Charlie hasn't even done this in NumPy!
# And it'd be more of a dataframe anywhere. But I guess behind the scenes one if just a proxy for the other
# anyway.

# Nice. Which dimension is which?

def build_backtrack_scoring_matrix(seq1, seq2, alphabet, scoring_matrix):

    i = len(seq1)
    j = len(seq2)

    backtrack_scoring_matrix = np.zeros((i + 1, j + 1), dtype=np.int8)

    value = 0

    for k in range(i):

        backtrack_scoring_matrix[0, k] = value

        value += c(seq1[k], "_", alphabet, scoring_matrix)

    value += c("_", "_", alphabet, scoring_matrix)

    backtrack_scoring_matrix[0, -1] = value

    value = 0

    for k in range(j):

        backtrack_scoring_matrix[k, 0] = value

        value += c(seq2[k], "_", alphabet, scoring_matrix)

    value += c("_", "_", alphabet, scoring_matrix)

    backtrack_scoring_matrix[-1, 0] = value

    return backtrack_scoring_matrix


def build_backtrack_matrix(i, j):
    # Vectorised string operations in https://docs.scipy.org/doc/numpy/reference/routines.char.html#module-numpy.char

    backtrack_matrix = np.full((i + 1, j + 1), b" ", dtype="string_")

    backtrack_matrix[:, 0] = "U"
    backtrack_matrix[0] = "L"
    backtrack_matrix[0, 0] = "E"

    return backtrack_matrix


# seq1 = "ABCACA"
# seq2 = "BAACBC"
# backtrack_matrix = [[]]
# backtrack_scoring_matrix = build_backtrack_scoring_matrix(len(seq1), len(seq2))
# backtrack_matrix = build_backtrack_matrix(len(seq1), len(seq2))


# # So we call populate_matrices and it goes from there. That effectively IS our dynprog
#
# best_alignment = []
# best1, best2, best_score = populate_matrices(seq1, seq2, scoring_matrix, backtrack_matrix)
# best_alignment.append(best1)
# best_alignment.append(best2)