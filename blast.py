import numpy as np
from itertools import product

def get_threshold(alphabet, scoring_matrix, s, t):

    frequencies = {i: (s.count(i) + t.count(i)) / len(s + t) for i in alphabet}
    pairs = list(product(alphabet, alphabet))
    weights = [frequencies[x] + frequencies[y] for (x, y) in pairs]

    return np.average([cost(alphabet, scoring_matrix, x, y) for (x, y) in pairs], weights=weights)

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]


def heuralign(alphabet, scoring_matrix, s, t):

    w = 2   # Word length.
    threshold = get_threshold(alphabet, scoring_matrix, s, t)
    lookup = get_lookup(s, w, threshold)

    print(lookup)

    segments = get_segments(alphabet, scoring_matrix, s, t, w, lookup, threshold)

    print(segments)

# Search each diagonal for two w length words such that score >= t.
# Okay, nice.

def get_lookup(s, w, threshold):

    lookup = dict()

    for y in range(len(s) - w + 1):

        word = s[y:y + w]

        if word not in lookup:

            lookup[word] = [y]

        else:

            lookup[word].append(y)

    return lookup


def get_segments(alphabet, scoring_matrix, s, t, w, lookup, threshold):

    segments = []

    for x in range(len(t) - w + 1):

        word = t[x:x + w]

        if word in lookup:

            for y in lookup[word]:

                segment = extend(alphabet, scoring_matrix, s, t, x, y, w, threshold)

                segments.append(segment)

    return segments

def score(alphabet, scoring_matrix, s, t, x, y, w):

    score = 0

    for i in range(w):

        score += cost(alphabet, scoring_matrix, s[y + i], t[x + i])

    return score


def extend(alphabet, scoring_matrix, s, t, x, y, w, threshold):

    segment_score = score(alphabet, scoring_matrix, s, t, x, y, w)

    up, down = 1, 1

    while up or down:

        if x == 0 or y == 0:

            up = 0

        else:

            new_score = segment_score + cost(alphabet, scoring_matrix, t[x - 1], s[y - 1])

            if new_score > segment_score:

                x -= 1
                y -= 1
                w += 1

                segment_score = new_score

            else:

                up = 0

        if x + w == len(t) or y + w == len(s):

            down = 0

        else:

            new_score = segment_score + cost(alphabet, scoring_matrix, t[x + 1], s[y + 1])

            if new_score > threshold:

                w += 1
                segment_score = new_score

            else:

                down = 0

    return x, y, w


alphabet = "ABC_"
scoring_matrix = np.array([
            [1, -1, -2, -1],
            [-1, 2, -4, -1],
            [-2, -4, 3, -2],
            [-1, -1, -2, 0]
])
s = "AABBAACA"
t = "CBACCCBA"

# t = threshold(alphabet, scoring_matrix, s, t)
#
# print(t)
#
#
heuralign(alphabet, scoring_matrix, s, t)