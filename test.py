from main import dynprog, dynproglin
from heuralign import heuralign

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

#BLAST
# alphabet = "ABC"
# scoring_matrix = [
#     [1,-1,-2,-1],
#     [-1,2,-4,-1],
#     [-2,-4,3,-2],
#     [-1,-1,-2,0]
# ]
# s = "AABCACA"
# t = "CBACACA"

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

#
# print("***** DYNPROG *****\n")
#
# a = dynprog(alphabet, scoring_matrix, s, t)
#
# print("Score:   ", a[0])
# print("Indices: ", a[1], a[2])
# print("")
#
print("**** DYNPROGLIN ****\n")

a = dynproglin(alphabet, scoring_matrix, s, t)

print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
print("")

# print("***** HEURALIGN *****\n")

# a = heuralign(alphabet, scoring_matrix, s, t)
#
# print("Score:   ", a[0])
# print("Indices: ", a[1], a[2])
# print("")