from main import dynprog, dynproglin, heuralign

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
                   "CCTAAG", "ACGGTAG"),
    "fasta": ("ABC", [
            [1, -1, -2, -1],
            [-1, 2, -4, -1],
            [-2, -4, 3, -2],
            [-1, -1, -2, 0]],
              "ABCCCABABACABCABCABCBAABABCCCAAACBCBCBABCABCBABBBCABCA", "AAACCBACBAC")
}

functions = [dynprog]

for function in functions:

    print("**** " + function.__name__.upper() + " ****\n")

    a = function(*tests["test.py"])

    print(a)

    print(a[3])
    print(a[4])
    print("Score:   ", a[0])
    print("Indices: ", a[1], a[2])
    print("")