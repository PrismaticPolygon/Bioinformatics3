from dynprog import dynprog
from dynproglin import dynproglin

#put ALL your code here
alphabet = "ABC"
scoring_matrix = [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]]
sequence1 = "AABBAACA"
sequence2 = "CBACCCBA"
#
# print("***** DYNPROG *****\n")

# a = dynprog(alphabet, scoring_matrix, sequence1, sequence2)

# print("Score:   ", a[0])
# print("Indices: ", a[1], a[2])

print("\n**** DYNPROGLIN ****\n")

a = dynproglin(alphabet, scoring_matrix, sequence1, sequence2)

# print("Score:   ", a[0])
# print("Indices: ", a[1],a[2])
