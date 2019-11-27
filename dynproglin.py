import numpy as np
from main import dynprog

# Okay. It still doesn't work.
# Can I verify...

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def rev(l):

    return l[::-1]

# For global alignment.
def hirschberg(alphabet, scoring_matrix, s, t):

    print("Hirschberg", s, t)

    # s == X
    # t == Y
    # s_align = Z
    # t_align = W

    s_align, t_align = "", ""

    if len(s) == 0: # Ah. So we fall through to here several times. That's why we're matching with gaps.

        for i in range(len(t)):

            s_align += "_"
            t_align += t[i]

    elif len(t) == 0:

        for i in range(len(s)):

            s_align += s[i]
            t_align += "_"

    elif len(s) == 1 or len(t) == 1:

        print("dynprog", s, t)

        _, _, _, s_align, t_align = dynprog(alphabet, scoring_matrix, s, t, local=False)

        print("dynprog_align", s_align, t_align)    # That's wrong. That is FUCKING WRONG.

    else:

        s_mid = int(len(s) / 2)

        print("s", s)
        print("left", s[:s_mid])
        print("right", s[s_mid:])

        print("t", t)

        score_l, _ = global_sublinear(alphabet, scoring_matrix, s[:s_mid], t)

        # print("Score_l: ", score_l)
        # No. We're deleting a character and not adding it back.

        score_r, _ = global_sublinear(alphabet, scoring_matrix, rev(s[s_mid:]), rev(t))

        # print("Score_r", score_r)

        t_mid = np.argmax(score_l + rev(score_r))

        # print("t_mid", t_mid)

        # Okay. I am happy that this follows the algorithm exactly.
        z_l, w_l = hirschberg(alphabet, scoring_matrix, s[:s_mid], t[:t_mid])
        z_r, w_r = hirschberg(alphabet, scoring_matrix, s[s_mid:], t[t_mid:])

        s_align = z_l + z_r
        t_align = w_l + w_r

        # print("s_align", s_align)
        # print("t_align", t_align)

    return s_align, t_align


def global_sublinear(alphabet, scoring_matrix, s, t):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    D = np.zeros((2, size_x), dtype="int16")

    max_i = None
    max_score = None

    for y in range(size_y):  # y

        for x in range(size_x):

            if y == 0 and x > 0:    # First row: the cost of matching t with all gaps

                D[0, x] = D[0, x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_")

            elif x == 0 and y > 0:  # First column: the cost of matching s with all gaps

                D[1, 0] = D[0, 0] + cost(alphabet, scoring_matrix, s[y - 1], "_")

            elif 0 < x < size_x and y > 0:

                D[1, x] = max(
                    D[1 - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),  # The cost of matching two characters (sub)
                    D[1 - 1][x] + cost(alphabet, scoring_matrix, s[y - 1],
                                         "_"), # The cost of matching a gap in t with a character in s (del)
                    D[1][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_") # The cost of matching a gap in s with a character in t (ins)
                )

                if max_i is None or D[1, x] > max_score:

                    max_i = y, x
                    max_score = D[1, x]

        if y > 0:

            D[0] = D[1].copy()

    return D[1], max_i



def global_quadratic(alphabet, scoring_matrix, s, t):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    D = np.zeros((size_y, size_x), dtype="int16")

    max_i = None
    max_score = None

    for y in range(size_y):  # y

        for x in range(size_x):

            if y == 0 and x > 0:    # First row: the cost of matching t with all gaps

                D[0, x] = D[0, x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_")

            elif x == 0 and y > 0:  # First column: the cost of matching s with all gaps

                D[y, 0] = D[y - 1, 0] + cost(alphabet, scoring_matrix, s[y - 1], "_")

            elif 0 < x < size_x and y > 0:

                D[y, x] = max(
                    D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),  # The cost of matching two characters (sub)
                    D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1],
                                         "_"), # The cost of matching a gap in t with a character in s (del)
                    D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_") # The cost of matching a gap in s with a character in t (ins)
                )

                if max_i is None or D[y, x] > max_score:

                    max_i = y, x
                    max_score = D[y, x]

    return D, max_i


















#
# alphabet = "ABC_"
# scoring_matrix = np.array([
#             [1, -1, -2, -1],
#             [-1, 2, -4, -1],
#             [-2, -4, 3, -2],
#             [-1, -1, -2, 0]
#         ])
# s = "BAAC"
# t = "BAC"

# s = "BA"
# t = "B"

alphabet = "ABCD_"
scoring_matrix = np.array([
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]])

s = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
t = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"


s_align, t_align = hirschberg(alphabet, scoring_matrix, s, t)

print(s_align)
print(t_align)

score, s_matches, t_matches, s_align, t_align = dynprog(alphabet, scoring_matrix, s, t, local=False)

# Okay. Let's try with some more complex strings.

# So what have I proven?
# Both work, globally.
# I need dynprog to work for both, and by default locally.

# So that works.
# Why doesn't mine?
# An off-by-one error?
# Maybe?

# Aha. If it's two characters, it's wrong.
# Why?
# It's global; we're not allowed to delete things.

print(s_align)
print(t_align)
# print(score)