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
                   "CCTAAG", "ACGGTAG"),
    "fasta": ("ABC", [
            [1, -1, -2, -1],
            [-1, 2, -4, -1],
            [-2, -4, 3, -2],
            [-1, -1, -2, 0]],
              "ABCCCABABACABCABCABCBAABABCCCAAACBCBCBABCABCBABBBCABCA", "AAACCBACBAC"),
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

functions = [dynprog, heuralign]

for function in functions:

    print("**** " + function.__name__.upper() + " ****\n")

    # for key, items in tests.items():

    a = function(*tests["heuristic3"])
    #
    print(a[3])
    print(a[4])
    print("Score:   ", a[0])
    print("Indices: ", a[1], a[2])
    print("")