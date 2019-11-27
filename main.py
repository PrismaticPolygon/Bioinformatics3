import numpy as np

# http://biorecipes.com/DynProgBasic/code.html
def dynprog(alphabet, scoring_matrix, s, t, local=True):

    # print("dynprog", local, s, t)

    alphabet = alphabet
    scoring_matrix = np.array(scoring_matrix)

    D, (end_y, end_x) = global_quadratic(alphabet, scoring_matrix, s, t)

    # print(end_y, end_x)

    # We're not passing through the parameter correctly.
    print(local)

    return traceback(alphabet, scoring_matrix, s, t, D, local=local)



def global_quadratic(alphabet, scoring_matrix, s, t):

    print(s, t)

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

    print(D)

    return D, max_i

# Smith Waterman is used for local
def dynproglin(alphabet, scoring_matrix, s, t):

    # print(s, t)

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)

    _, (end_y, end_x) = score_matrix(alphabet, scoring_matrix, s, t, sublinear=True)

    # print(end_y, end_x)
    #
    # print(s[:end_y], t[:end_x])

    s = rev(s[:end_y])
    t = rev(t[:end_x])

    # print(s, t)

    _, (start_y, start_x) = score_matrix(alphabet, scoring_matrix, s, t, sublinear=True)

    # print(start_y, start_x)

    # print(s[:start_y], t[:start_x])
    # This doesn't work. Must be some sublinear BS. Or my global shit is fucked.

    s = rev(s[:start_y])
    t = rev(t[:start_x])

    # print(s, t)

    # So it's a flaw with my hirschberg, or backtrack; I pass in the right fields.
    # Somehow, we're

    s_align, t_align = hirschberg(alphabet, scoring_matrix, s, t)

    # print(s_align, t_align)

    score = align_score(alphabet, scoring_matrix, s_align, t_align)

    s_matches, t_matches = get_alignment_indices(s_align, t_align)

    s_matches += end_y - start_y
    t_matches += end_x - start_x

    return score, s_matches, t_matches, s_align, t_align

# file:///C:/Users/user/Documents/MEGA/University/Year%203/CCS/Bioinformatics/HeuristicAlign.pdf
def heuralign(alphabet, scoring_matrix, s, t):

    #linspace return evenly spaced numbers over a specified interval
    # numpy ufuncs have a reduceat function.

    # s (database), |s| >> |t| is the y-axis
    # t (query), is the x-axis

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)
    ktup = 4
    k = 9

    lookup = build_lookup(alphabet, s, ktup)    # Build lookup table of database string s
    matches = build_matches(alphabet, scoring_matrix, s, t, ktup, lookup)   # Build list of extended alignments from query string t
    diagonals = build_diagonals(k, matches) # Sort into diagonals.

    for key, matches in diagonals.items():

        start_x, start_y = len(t) - 1, len(s) - 1
        end_x, end_y = 0, 0

        for x, y, length in matches:

            if x < start_x:

                start_x = x

            if x + length > end_x:

                end_x = x + length

            if y < start_y:

                start_y = y

            if y + length > end_y:

                end_y = y + length

        s = s[start_y:end_y]
        t = t[start_x:end_x]

        D = banded(alphabet, scoring_matrix, s, t, k)

        return traceback(alphabet, scoring_matrix, s, t, D)

# Perhaps it never worked.

# https://en.wikipedia.org/wiki/Hirschberg's_algorithm
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

    else:

        s_mid = int(len(s) / 2)

        print("s", s)
        print("left", s[:s_mid])
        print("right", s[s_mid:])

        print("t", t)

        score_l, _ = score_matrix(alphabet, scoring_matrix, s[:s_mid], t, local=False, sublinear=False)

        # print("Score_l: ", score_l)
        # No. We're deleting a character and not adding it back.

        score_r, _ = score_matrix(alphabet, scoring_matrix, rev(s[s_mid:]), rev(t), local=False, sublinear=False)

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

def score_matrix(alphabet, scoring_matrix, s, t, local=True, sublinear=False):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    shape = (2, size_x) if sublinear else (size_y, size_x)

    D = np.zeros(shape, dtype="int16")

    max_i = None
    max_score = None

    for y in range(size_y):  # y

        for x in range(size_x):

            if not local and y == 0 and x > 0:    # First row: the cost of matching t with all gaps

                D[0, x] = D[0, x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_")

            elif not local and x == 0 and y > 0:  # First column: the cost of matching s with all gaps

                D_y = 1 if sublinear else y

                D[D_y, 0] = D[D_y - 1, 0] + cost(alphabet, scoring_matrix, s[y - 1], "_")

            elif 0 < x < size_x and y > 0:

                D_y = 1 if sublinear else y

                D[D_y, x] = max(
                    D[D_y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),  # The cost of matching two characters (sub)
                    D[D_y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1],
                                         "_"), # The cost of matching a gap in t with a character in s (del)
                    D[D_y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_") # The cost of matching a gap in s with a character in t (ins)
                )

                D[D_y, x] = 0 if (local and D[D_y, x] < 0) else D[D_y, x]

                if max_i is None or D[D_y, x] > max_score:

                    max_i = y, x
                    max_score = D[D_y, x]

        if y > 0 and sublinear:

            D[0] = D[1].copy()

    return (D[1] if sublinear else D), max_i

def traceback(alphabet, scoring_matrix, s, t, D, local=True):

    score = np.amax(D) if local else D[-1][-1]
    y, x = np.unravel_index(D.argmax(), D.shape) if local else (len(s), len(t))

    s_align, t_align = "", ""
    s_matches, t_matches = [], []

    while y != 0 or x != 0:

        current = D[y][x]

        if local and current == 0:  # The end of the best local alignment

            break

        elif current == D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]): # Match (D)

            x -= 1
            y -= 1

            s_align = s[y] + s_align
            t_align = t[x] + t_align

            s_matches.append(y)
            t_matches.append(x)

        elif current == D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"):  # Matching a gap in s with a character in t (L)

            x -= 1

            s_align = "_" + s_align
            t_align = t[x] + t_align

        elif current == D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_"):  # Matching a gap in t with a character in s (U)

            y -= 1

            s_align = s[y] + s_align
            t_align = "_" + t_align

        else:

            raise ValueError("Something's fucked!")

    s_matches = np.array(rev(s_matches))
    t_matches = np.array(rev(t_matches))

    return score, s_matches, t_matches, s_align, t_align

def get_alignment_indices(s_align, t_align):

    assert(len(s_align) == len(t_align))

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

    return np.array(s_matches), np.array(t_matches)

def align_score(alphabet, scoring_matrix, s_align, t_align):

    assert(len(s_align) == len(t_align))

    score = 0

    for i in range(len(s_align)):

        score += cost(alphabet, scoring_matrix, s_align[i], t_align[i])

    return score

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def rev(l):

    return l[::-1]

# alphabet = "ABC_"
# scoring_matrix = np.array([
#             [1, -1, -2, -1],
#             [-1, 2, -4, -1],
#             [-2, -4, 3, -2],
#             [-1, -1, -2, 0]
#         ])
# s = "BAAC"
# t = "BAC"
#
# s_align, t_align = hirschberg(alphabet, scoring_matrix, s, t)
#
# print(s_align)
# print(t_align)

# I'm sick and tired of this.

# score = align_score(alphabet, scoring_matrix, s_align, t_align)
#
# s_matches, t_matches = get_alignment_indices(s_align, t_align)



# Okay. The two strings in question are BAAC and BAC
# I just need to verify that hirschberg works on them