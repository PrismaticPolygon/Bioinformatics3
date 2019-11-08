import numpy as np
from main import dynprog


#put ALL your code here


a = dynprog("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
# a = dynprog("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "ABCACA", "BAACB")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
