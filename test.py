from main import dynprog, dynproglin
from heuralign import heuralign


tests = {
    "Duo": ("ABC", [
            [1,-1,-2,-1],
            [-1,2,-4,-1],
            [-2,-4,3,-2],
            [-1,-1,-2,0]
            ], "ABCACA", "BAACB"),
    "test.py": ("ABC", [
            [1, -1, -2, -1],
            [-1, 2, -4, -1],
            [-2, -4, 3, -2],
            [-1, -1, -2, 0]
        ], "AABBAACA", "CBACCCBA"),
    "biorecipes": ("ACGT", [
            [2,-1,1,-1,-2],
            [-1,2,-1,1,-2],
            [1,-1,2,-1,-2],
            [-1,1,-1,2,-2],
            [-2,-2,-2,-2,0]],
                   "CCTAAG", "ACGGTAG")
}

# Expected: 5 [3, 5, 6] [1, 2, 3]. Should be the same for both.

functions = [dynprog, dynproglin]

# Yeah, I don't trust that one bit.
# It'll be very hard to pin down, though.

for function in functions:

    print("**** " + function.__name__.upper() + " ****\n")

    a = function(*tests["test.py"])

    print(a[3])
    print(a[4])
    print("Score:   ", a[0])
    print("Indices: ", a[1], a[2])
    print("")