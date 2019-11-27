import numpy as np

# http://biorecipes.com/DynProgBasic/code.html
def dynprog(alphabet, scoring_matrix, s, t, local=True):

    alphabet = alphabet if alphabet[-1] == "_" else alphabet + "_"
    scoring_matrix = np.array(scoring_matrix)

    D, _ = score_matrix(alphabet, scoring_matrix, s, t, local=local)

    return traceback(alphabet, scoring_matrix, s, t, D, local=local)

def dynproglin(alphabet, scoring_matrix, s, t):

    alphabet = alphabet if alphabet[-1] == "_" else alphabet + "_"
    scoring_matrix = np.array(scoring_matrix)

    _, (end_y, end_x) = score_matrix(alphabet, scoring_matrix, s, t, local=True, linear=True)

    s = rev(s[:end_y])
    t = rev(t[:end_x])

    _, (start_y, start_x) = score_matrix(alphabet, scoring_matrix, s, t, local=True, linear=True)

    s = rev(s[:start_y])
    t = rev(t[:start_x])

    s_align, t_align = hirschberg(alphabet, scoring_matrix, s, t)

    score = 0
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

        score += cost(alphabet, scoring_matrix, s_align[i], t_align[i])

    s_matches = np.array(s_matches) + end_y - start_y
    t_matches = np.array(t_matches) + end_x - start_x

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

        _, _, _, s_align, t_align = dynprog(alphabet, scoring_matrix, s, t, local=False)

    else:

        s_mid = int(len(s) / 2)

        score_l, _ = score_matrix(alphabet, scoring_matrix, s[:s_mid], t, linear=True, local=False)

        score_r, _ = score_matrix(alphabet, scoring_matrix, rev(s[s_mid:]), rev(t), linear=True, local=False)

        t_mid = np.argmax(score_l + rev(score_r))

        z_l, w_l = hirschberg(alphabet, scoring_matrix, s[:s_mid], t[:t_mid])
        z_r, w_r = hirschberg(alphabet, scoring_matrix, s[s_mid:], t[t_mid:])

        s_align = z_l + z_r
        t_align = w_l + w_r

    return s_align, t_align

def score_matrix(alphabet, scoring_matrix, s, t, linear=False, local=True):

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    shape = (2, size_x) if linear else (size_y, size_x)

    D = np.zeros(shape, dtype="int16")

    max_i = None
    max_score = None

    for y in range(size_y):  # y

        for x in range(size_x):

            if not local and y == 0 and x > 0:    # First row: the cost of matching t with all gaps

                D[0, x] = D[0, x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_")

            elif not local and x == 0 and y > 0:  # First column: the cost of matching s with all gaps

                D_y = 1 if linear else y

                D[D_y, 0] = D[D_y - 1, 0] + cost(alphabet, scoring_matrix, s[y - 1], "_")

            elif 0 < x < size_x and y > 0:

                D_y = 1 if linear else y

                D[D_y, x] = max(
                    D[D_y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),  # The cost of matching two characters (sub)
                    D[D_y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_"), # The cost of matching a gap in t with a character in s (del)
                    D[D_y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_") # The cost of matching a gap in s with a character in t (ins)
                )

                D[D_y, x] = 0 if (local and D[D_y, x] < 0) else D[D_y, x]

                if max_i is None or D[D_y, x] > max_score:

                    max_i = y, x
                    max_score = D[D_y, x]

        if linear and y > 0:

            D[0] = D[1].copy()

    return (D[1] if linear else D), max_i

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

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def rev(l):

    return l[::-1]