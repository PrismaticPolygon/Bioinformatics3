from dynprog import dynprog
from dynproglin import dynproglin

# Duo
alphabet = "ABC"
scoring_matrix = [
    [1,-1,-2,-1],
    [-1,2,-4,-1],
    [-2,-4,3,-2],
    [-1,-1,-2,0]
]
s = "ABCACA"
t = "BAACB"

# test.py
# alphabet = "ABC"
# scoring_matrix = [
#     [1, -1, -2, -1],
#     [-1, 2, -4, -1],
#     [-2, -4, 3, -2],
#     [-1, -1, -2, 0]
# ]
#
# s = "AABBAACA"
# t = "CBACCCBA"

# Biorecipes
# alphabet = "ACGT"
# scoring_matrix = [
#     [2,-1,1,-1,-2],
#     [-1,2,-1,1,-2],
#     [1,-1,2,-1,-2],
#     [-1,1,-1,2,-2],
#     [-2,-2,-2,-2,0]
# ]
# s = "CCTAAG"
# t = "ACGGTAG"


# print("***** DYNPROG *****\n")

# sequence1 = "A"
# sequence2 = "CBACC"

a = dynprog(alphabet, scoring_matrix, s, t)

# print("Score:   ", a[0])
# print("Indices: ", a[1], a[2])

# print("\n**** DYNPROGLIN ****\n")

# a = dynproglin(alphabet, scoring_matrix, sequence1, sequence2)

# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])
