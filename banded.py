from itertools import product
from main import dynprog
import numpy as np

# Real-life matches often contain long strings with gap-less matches.
# Heuristics try to find significant gap-less matches and extend them.

def banded(s, t, matches):

    k = 10

    # print(matches)

    matches.sort(key=lambda x: x[1] - x[0])

    matches = list(filter(lambda x: (x[1] - x[0] < k), matches))

    # Get the top-left corner.

    top_left = matches[0][0]
    bottom_right = k + 4    # Should really be worth length

    s = s[top_left: bottom_right]
    t = t[top_left: bottom_right]

    # But then it's still only local. NOT global. Pe

    score, s_align, t_align, s_matches, t_matches = dynprog("ABC", scoring_matrix, s, t)

    s_align = [align + top_left for align in s_align]
    t_align = [align + top_left for align in t_align]

    print(score)
    print(s_align)
    print(t_align)
    print(s_matches)
    print(t_matches)

    # This seems pretty good to me at the moment.

    # So now we search the diagonal band of the matrix.
    # We've already chosen k.
    # Plus we know how big it needs to be: (10 + ktup) * (10 + ktup) at maximum. Very nice.
    # So now we backtrack on that, using dynprog

    # I could do my own simple sort.



    # If the optimal alignment of s and t has few gaps, then the path of the alignment will be close to some diagonal.
    # To find such a path, it is sufficient to search only a diagaonl band of the matrix.
    # If the diagonal band consists of k diagonals (i.e. is of width k) then dynamic programming takes O(kn).


    # So we only have to build our scoring matrix for that narrow band.
    # What do we treat the rest as? O?

    # For diagonals, i and j are constant. Then, for each, we modify k.

    # V[i, i + k/2]         Out of range.
    # V[i + 1, i + k / 2]   V[i + 1, i + k / 2 + 1]


    # But then we have another problem. Where is the banded diagonal?
    # It doesn't need to be the main diagonal when looking for a good local alignment (or when m and n are very different)

    # So: heuristically find potential diagonals and evaluate them using Banded DP.
    # So I guess that's the first bit.
    # Such a heurst


def cost(c1, c2):

    global alphabet, scoring_matrix

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

alphabet = "ABC"

scoring_matrix = np.array([
    [1, -1, -2, -1],
    [-1, 2, -4, -1],
    [-2, -4, 3, -2],
    [-1, -1, -2, 0]
])

#
# Max
# size
# of
# lookup is | alphabet | ^ ktup.For
# small
# ktup, the
# entire
# table is stored;


# for large ktup, only entries for tuples actually found in the database should be kept, so that the index table
# size is at most linear. In this case hashing is needed (why?)

# The index table is prepared for each database sequence ahead of users' matching requests, at compilation time
# So matching time is O(|T| x max (row_length))

# Then we expand to banded DP. But how? I have a list of matches. Each has an i, j position dictating where the match
# starts.

def score_match(s, t, match):

    x_start, y_start, length = match

    score = 0

    for i in range(length):

        score += cost(s[y_start + i], t[x_start + i])

    return score


def fasta(s, t, ktup):

    # s is the database string. s is the y-axis
    # t is the query string. t is the x-axis
    # ktup: seed size.

    lookup = {"".join(i): [] for i in product(alphabet, repeat=ktup)}

    for i in range(len(s) - ktup + 1):

        word = s[i:i + ktup]
        lookup[word].append(i)

    matches = []

    for i in range(len(t) - ktup + 1):

        word = t[i:i + ktup]

        for j in lookup[word]:

            matches.append((i, j, ktup))

    def extend(match):

        # Extend along the diagonal as long as the score improves (and never below some score value)

        x_start, y_start, length = match
        x_end, y_end = x_start + length, y_start + length

        while x_start > 0 and y_start > 0:

            if cost(t[x_start], s[y_start]) > 0:

                x_start -= 1
                y_start -= 1

            else:

                break

        while x_end < len(t) and y_end < len(s):

            if cost(t[x_end], s[y_end]) > 0:

                x_end += 1
                y_end += 1

            else:

                break

        return x_start, y_start, x_end - x_start

    # print(lookup)
    # print(matches)
    # print([score_match(s, t, match) for match in matches])

    extended = [extend(match) for match in matches]

    banded(s, t, extended)

    # Identifying potential diagonals: high-scoring gap-less diagonals contain several 'seeds' of perfect matches.
    # Since the alginment is gapless, all perfect match regions reside on the same diagonal.

    # Let ktup be a parameter denoting the seed length of interest.
    # Assume that s is the database and t is the query string, e.g. |s| >> |t|

    # Prepare an index table which, for any sequence x of length ktup, contains the list of x's positions in S.


    pass



fasta("ABCCCABABACABCABCABCBAABABCCCAAACBCBCBABCABCBABBBCABCA", "AAACCBACBAC", 3)