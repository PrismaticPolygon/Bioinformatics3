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
              "ABCCCABABACABCABCABCBAABABCCCAAACBCBCBABCABCBABBBCABCA", "BCCCABACBACBCABCA"),
    "heuristic1": ("ABCD", [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]],
                   "AAAAACCDDCCDDAAAAACC", "CCAAADDAAAACCAAADDCCAAAA"),
    "heuristic2": ("ABCD", [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]],
                   "AACAAADAAAACAADAADAAA", "CDCDDD"),
    "heuristic3": ("ABCD", [
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]],
                   "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
                   "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
}

functions = [dynprog, dynproglin, heuralign]

# for key, items in tests.items():
# Definitely not perfect.
# Okay. Peter's is better. Let's groom his into my format.

for key, variables in tests.items():

    for function in functions:

        variables = tests[key]

        print("**** " + function.__name__.upper() + " ****\n")

        # print("")
        # It is definitely faster. No I groom it to match my style of code.

        a = function(*variables)

        # print(a)
        #
        # print(a[3])
        # print(a[4])
        print("Score:   ", a[0])
        print("Indices: ", a[1], a[2])
        print("")