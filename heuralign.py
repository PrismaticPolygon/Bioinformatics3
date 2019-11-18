import numpy as np

# Real-life matches often contain long strings with gap-less matches.

# And use ord

# Use RE to find string matches: interesting idea!

# Polynomial rolling hashing function
def hash(s):

    # https://cp-algorithms.com/string/string-hashing.html
    # p should be roughly equal to the number of characters in the input alphabet
    p = 31

    # m should be a large number: the probability of two strings colliding is roughly 1 / m.
    m = 10 ** 9 + 9

    hash_value = 0
    p_pow = 1

    for char in s:

        hash_value = (hash_value + ord(char) * p_pow) % m
        p_pow = (p_pow * p) % m

    return hash_value


def chain():

    pass
    # Find a collection of non-contradicting sub-alignments that maximise some scoring function.



def heuralign(alphabet, scoring_matrix, s, t):

    scoring_matrix = np.array(scoring_matrix)
    ktup = 3
    lookup = dict()
    size_x = len(t)
    size_y = len(s)

    for y in range(len(s) - ktup + 1):

        hash_word = hash(s[y:y + ktup])

        lookup[hash_word] = [y] if hash_word not in lookup else lookup[hash_word] + [y]

    matches = dict()

    for x in range(len(t) - ktup + 1):

        hash_word = hash(t[x:x + ktup])

        if hash_word in lookup:

            for y in lookup[hash_word]:

                diagonal = y - x

                if diagonal not in matches:

                    matches[diagonal] = [(x, y)]

                else:

                    matches[diagonal].append((x, y))

    for match in matches:

        print(match)


def extend(alphabet, scoring_matrix, s, t, x, y, length):

    up, down = 1, 1

    while up or down:

        if x == 0 or y == 0 or cost(alphabet, scoring_matrix, t[x - 1], s[y - 1]) < 0:

            up = 0

        else:

            x -= 1
            y -= 1
            length += 1

        if x + length == len(t) or y + length == len(s) or cost(alphabet, scoring_matrix, t[x + 1], s[y + 1]) < 0:

            down = 0

        else:

            length += 1

    return x, y, length

def cost(alphabet, scoring_matrix, c1, c2):

    # print(c1, c2)

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

alphabet = "ABCD"
scoring_matrix = [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]]
s = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
t = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

x, y, length = 3, 3, 3

# extended = extend(alphabet, scoring_matrix, s, t, x, y, length)
# print(extended)
# We want more than the i

heuralign(alphabet, scoring_matrix, s, t)


# The main modification is a shearing algorithm to save space.
# Looks like I'm going to be doing everything simultaneously.

