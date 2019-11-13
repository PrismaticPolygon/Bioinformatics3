from itertools import product
from main import dynprog
import numpy as np

# Real-life matches often contain long strings with gap-less matches.
# Heuristics try to find significant gap-less matches and extend them.


# for large ktup, only entries for tuples actually found in the database should be kept, so that the index table
# size is at most linear. In this case hashing is needed (why?)

# The index table is prepared for each database sequence ahead of users' matching requests, at compilation time
# So matching time is O(|T| x max (row_length))

# Then we expand to banded DP. But how? I have a list of matches. Each has an i, j position dictating where the match
# starts.

def align(s, t, matches, k):

    start = min(matches, key=lambda m: min(m[0], m[1]))
    end = max(matches, key=lambda m: max(m[0], m[1]))

    s_align = ""
    t_align = ""
    score = 0

    print(s)
    print(t)
    print(matches)

    match = matches[0]
    #
    # Out of bounds? Because I'm supposed to match with a 0?
    # This should be dynamic programming for sure: I just can't figure out how.
    # Actually, is it that hard?
    for i in range(match[2]):

        s_i = s[i + match[0]]
        t_i = t[i + match[1]]

        score += cost(s_i, t_i)

        s_align += s_i
        t_align += t_i
    #
    #
    # print(start, end)
    # print(s_align, t_align, score)

    return 0, 0


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

        if sublinear and y < size_y - 1: # Copy the 1st row onto the second row unless it's the final iteration

            D[0] = D[1].copy()
            D[1] = 0

    return D[1] if sublinear else D


def banded(s, t, matches):

    diagonals = dict()
    k = 9

    # I'm not sure how. Actually, that's untrue. It'll be the nearest (under) multiple of k.
    # So then I guess part of the heuristic is identifying the diagonals! I'll do the most simple kind of banded for now:

    for match in matches:

        x, y, length = match

        i = k * ((y - x) // k)

        if i not in diagonals:

            diagonals[i] = [match]

        else:

            diagonals[i].append(match)

    for key, value in diagonals.items():

        s_align, t_align = align(s, t, value, k)

        break


    # And now I do gapped alignment of them. It doesn't actually have to use dynprog either.
    # I could just use all of them.
    # Again: it has to be sub-quadratic, and results have to be decent.

    # So it's the y-offset.
    # Yes. (0, 29, 10), (4, 36, 6) and (7, 36, 3) all lie in the same
    # bad, which lies between (0, 27) and (0, 36).
    # So k is like the offset now.
    # But then: what does a value of 27 mean?
    # It has to be taken from... not even (0, 0), actually.
    # It could start anywhere.
    # It's irrelevant, actually.
    # They're essentially gradients. There we go. THAT's how to imagine it.

    # Well, I don't have to for now.



    # And then we sort by score.

    # Long strings with gapless sections. But I've got the k wrong here.
    # If it's out-of-bounds, count it as a 0, local-alignemtn style. That's what we're doing, after all.

    # k = 9

    # Heuristically find potential diagonals and evaluate them using Banded DP.
    # So we should actually group the matches.
    # We can again create buckets.

    # print(matches)

    # matches.sort(key=lambda x: x[1] - x[0])
    #
    # matches = list(filter(lambda x: (x[1] - x[0] < k), matches))

    # Get the top-left corner.

    # We can't use dynprog. It has to be a modified version of our current scoring_matrix
    # whereby we check for out of bounds.
    # And so k has to be an odd number.

    # top_left = matches[0][0]
    # bottom_right = k + 4    # Should really be worth length
    #
    # s = s[top_left: bottom_right]
    # t = t[top_left: bottom_right]
    #
    # # But then it's still only local. NOT global. Pe
    #
    # score, s_align, t_align, s_matches, t_matches = dynprog("ABC", scoring_matrix, s, t)
    #
    # s_align = [align + top_left for align in s_align]
    # t_align = [align + top_left for align in t_align]
    #
    # print(score)
    # print(s_align)
    # print(t_align)
    # print(s_matches)
    # print(t_matches)

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


def score(s, t, match):

    x, y, length = match

    score = 0

    for i in range(length):

        score += cost(s[y + i], t[x + i])

    return score

def extend(s, t, match):

    x, y, length = match

    while x > 0 and y > 0:

        if cost(t[x - 1], s[y - 1]) > 0:

            x -= 1
            y -= 1
            length += 1

        else:

            break

    while x + length + 1 < len(t) and y + length + 1 < len(s):

        if cost(t[x + 1], s[y + 1]) > 0:

            length += 1

        else:

            break

    return x, y, length

def fasta(s, t, ktup):

    # s (database), |s| >> |t| is the y-axis
    # t (query), is the x-axis

    lookup = {"".join(i): [] for i in product(alphabet, repeat=ktup)}

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]
        lookup[word].append(y)

    matches = []

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        for y in lookup[word]:

            matches.append((x, y, ktup))

    extended = [extend(s, t, match) for match in matches]

    banded(s, t, extended)

def cost(c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

alphabet = "ABC"

scoring_matrix = np.array([
    [1, -1, -2, -1],
    [-1, 2, -4, -1],
    [-2, -4, 3, -2],
    [-1, -1, -2, 0]
])


fasta("ABCCCABABACABCABCABCBAABABCCCAAACBCBCBABCABCBABBBCABCA", "AAACCBACBAC", 3)