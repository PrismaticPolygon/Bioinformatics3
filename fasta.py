import numpy as np
from itertools import product
from main import traceback

# Compare a query string against a single text string
# Based on the assumption that good local alignment is likely to have some exact matching subsequences
# I can't do this right now.
# I'll have a look at CV.

# Lots of space. Might be

# Distance between hot spots is negative and length dependant. Find and store the 10 best diagonal rnus.
# I'm fairly sure that my banded DP works.

def extend(alphabet, scoring_matrix, s, t, seed):

    x, y, length = seed

    # So. For each seed, we evaluate the score of that diagonal.
    # Now what?
    # But then we're back to the problem of combining them again!
    # And runs of identities shorter than the size of a word.
    # The ends of each region are trimmed to inclde only hose resides that contribute
    # to the highest score the region.


    up, down = 1, 1

    while up or down:

        if x == 0 or y == 0:

            up = 0

        else:

            cost(alphabet, scoring_matrix, t[x - 1], s[y - 1])

            if cost(alphabet, scoring_matrix, t[x - 1], s[y - 1]) > 0:

                x -= 1
                y -= 1
                length += 1

            else:

                up = 0

        if x + length == len(t) or y + length == len(s):

            down = 0

        else:

            if cost(alphabet, scoring_matrix, t[x + 1], s[y + 1]) > 0:

                length += 1

            else:

                down = 0

    return x, y, length

def get_threshold(alphabet, scoring_matrix, s, t):

    frequencies = {i: (s.count(i) + t.count(i)) / len(s + t) for i in alphabet}
    pairs = list(product(alphabet, alphabet))
    weights = [frequencies[x] + frequencies[y] for (x, y) in pairs]

    return np.average([cost(alphabet, scoring_matrix, x, y) for (x, y) in pairs], weights=weights)

def cost(alphabet, scoring_matrix, c1, c2):

    return scoring_matrix[alphabet.index(c1), alphabet.index(c2)]

def get_score(alphabet, scoring_matrix, s, t):

    score = 0

    for i in range(min(len(s), len(t))):

        score += cost(alphabet, scoring_matrix, s[i], t[i])

    return score

def get_lookup(s, ktup):

    lookup = dict()

    for y in range(len(s) - ktup + 1):

        word = s[y:y + ktup]

        if word not in lookup:

            lookup[word] = [y]

        else:

            lookup[word].append(y)

    return lookup
#
# def get_diagonals(alphabet, scoring_matrix, s, t, ktup, lookup):
#
#     threshold = get_threshold(alphabet, scoring_matrix, s, t)
#     diagonals = dict()
#     scores = dict()
#
#     for x in range(len(t) - ktup + 1):
#
#         word = t[x:x + ktup]
#
#         if word in lookup:
#
#             for y in lookup[word]:
#
#                 diagonal = x - y
#
#                 score = get_score(alphabet, scoring_matrix, word, word)
#
#                 if diagonal not in diagonals:
#
#                     diagonals[diagonal] = [(x, y, ktup)]
#                     scores[diagonal] = score
#
#                 else:
#
#                     prev = diagonals[diagonal][-1]
#
#                     diagonals[diagonal].append((x, y, ktup))
#
#                     # Cause threshold is negative.
#                     # Ah. But we do want my previous combinatorial code.
#                     # Oof. I'm close to giving up.
#                     # Solution? Abs
#
#                     scores[diagonal] += score - abs((x - prev[0] + ktup) * threshold)
#
#
#     print(diagonals)
#     print(scores)
#
#     return diagonals

def build_diagonals(alphabet, scoring_matrix, s, t, matches, ktup):

    diagonals = dict()

    for diagonal, match in matches:

        start = match[0]
        end = match[-1]

        diag_start = extend(alphabet, scoring_matrix, s, t, (*start, ktup))
        end_x, end_y, length = extend(alphabet, scoring_matrix, s, t, (*end, ktup))

        # Like so. We both have all of the parts; it's just putting them together.
        # Also has the threshold val bollocks.

        diag_end = [end_x + length, end_y + length]

        diag_score = score_diagonal(alphabet, scoring_matrix, s, t, diag_start, diag_end, diag_start[1])

    return diagonals


def score_diagonal(alphabet, scoring_matrix, s, t, diag_start, diag_end, q_start) :

    diag_score = 0
    index = 0

    while index + diag_start < diag_end:

        diag_score += cost(alphabet, scoring_matrix, s[diag_start + index], t[q_start + index])

        index += 1

    return diag_score



def get_hotspots(t, ktup, lookup):

    hotspots = dict()

    for x in range(len(t) - ktup + 1):

        word = t[x:x + ktup]

        if word in lookup:

            for y in lookup[word]:

                diagonal = y - x

                if word not in hotspots:

                    hotspots[word] = [(x, y)]

                else:

                    hotspots[word].append((x, y))

    return hotspots

def score_seed_II(alphabet, scoring_matrix, s, t, seed):

    start_x, start_y, end_x, end_y = seed
    diagonal = start_y - start_x

    # print("Scoring", seed)

    # start = max(0, start_x, start_y)
    # end = min(len(s) - start, len(t) - start, end_x - start, end_y - start)
    #
    #
    score = 0

    for i in range(start_x, end_x):

        try:

            score += cost(alphabet, scoring_matrix, s[diagonal + i], t[diagonal + i])

        except IndexError:

            return score

    return score


def score_seed(alphabet, scoring_matrix, s, t, seed):

    x, y, length = seed
    score = 0

    for i in range(length):

        score += cost(alphabet, scoring_matrix, s[y + i], t[x + i])

    return score

def get_ktup(alphabet):

    return max(1, int(-0.25 * len(alphabet)) + 7)

    # ktup is 6 for an alphabet of size 4 (DNA), and 2 for an alphabet for size 20 (proteins)
    # y = -0.25x + 7

def get_diagonals(hotspots):

    diagonal_runs = dict()

    for word, positions in hotspots.items():

        # print(word, positions)

        for x, y in positions:

            diagonal = y - x

            if diagonal not in diagonal_runs:

                diagonal_runs[diagonal] = [(x, y)]

            else:

                diagonal_runs[diagonal].append((x, y))

    return diagonal_runs

# I don't care any more.
# It works, right?

def FASTA(alphabet, scoring_matrix, s, t):

    alphabet += "_"
    scoring_matrix = np.array(scoring_matrix)
    ktup = get_ktup(alphabet)

    # s (the y-axis) is the database string, |s| >> |t|. t (the x-axis), is the query string. If t > s, switch them.
    if len(s) < len(t):

        print("Switching s and t")

        x = s
        s = t
        t = x

    # 1. Identify hot-spots. If we can't find any matches with our current ktup, reduce it until we do. If ktup = 0...
    while 1:

        print("ktup", ktup)

        lookup = get_lookup(s, ktup)

        print("Lookup", lookup)

        hotspots = get_hotspots(t, ktup, lookup)

        print("Hotspots", hotspots)

        if len(hotspots) == 0:

            ktup -= 1

        else:

            break

    # 2. Build the diagonals. Score them, and get the top 10.
    diagonal_runs = get_diagonals(hotspots)

    print("Diagonals", diagonal_runs)

    average_match = 3
    average_gap = -2

    scores = dict()

    for diagonal, hotspots_ in diagonal_runs.items():

        if len(hotspots_) == 1:

            scores[diagonal] = average_match * ktup    # Oh. Because they are all 1. So that's actually somewhat misleading.

        else:

            diagonal_score = 1

            for i in range(1, len(hotspots_)):

                prev_x, prev_y = hotspots_[i - 1]
                cur_x, cur_y = hotspots_[i]

                if prev_x + ktup > cur_x:

                    distance = cur_x - prev_x

                else:

                    distance = cur_x - prev_x - ktup

                print("Distance", distance)

                diagonal_score += average_match * ktup + distance * average_gap

            scores[diagonal] = diagonal_score

    print("Scores", scores)

    sorted_diagonals = sorted(diagonal_runs.keys(), key=lambda x: scores[x], reverse=True)[:10]

    print("Sorted diagonals", sorted_diagonals)

    # 3. Evaluate using a substitution matrix, and identify initial regions.

    D = np.full((len(s) + 1, len(t) + 1), fill_value=" ")

    for diagonal, hotspots in diagonal_runs.items():

        for x, y in hotspots:

            for i in range(ktup):

                D[y + i, x + i] = "X"

    print(D)

    initial_regions = dict()

    for diagonal in sorted_diagonals:

        hotspots_ = diagonal_runs[diagonal]
        best_score = 0
        initial_region = None

        for x, y in hotspots_:

            hotspot = extend(alphabet, scoring_matrix, s, t, (x, y, ktup))

            score = score_seed(alphabet, scoring_matrix, s, t, hotspot)

            if initial_region is None or score > best_score:

                initial_region = hotspot

        initial_regions[diagonal] = initial_region

    # What's this bit?

    initial_regions = list(initial_regions.values())

    print(initial_regions)

    initial_regions = [(start_x, start_y, start_x + length, start_y + length) for (start_x, start_y, length)in initial_regions]

    print("Initial regions:", initial_regions)
    #
    def manhattan(region_a, region_b):

        a, b, _, _ = region_a
        x, y, _, _ = region_b

        return abs(a - x) + abs(b - y)

    def join(initial_regions):

        # print("Joining: ", initial_regions)

        for i in range(len(initial_regions)):

            region_a = initial_regions[i]
            score_a = score_seed_II(alphabet, scoring_matrix, s, t, region_a)

            for j in range(i + 1, len(initial_regions)):

                region_b = initial_regions[j]
                score_b = score_seed_II(alphabet, scoring_matrix, s, t, region_b)

                manhattan_distance = manhattan(region_a, region_b)

                # print("Region a: {} ({})".format(region_a, score_a))
                # print("Region b: {} ({})".format(region_b, score_b))
                # print("Manhattan: {}".format(manhattan_distance))

                if score_a + score_b - manhattan(region_a, region_b) * average_gap > 0:

                    # So now we get start AND end points. No need for length.

                    start_x = min(region_a[0], region_b[0])
                    start_y = min(region_a[1], region_b[1])

                    end_x = max(region_a[2], region_b[2])
                    end_y = max(region_a[3], region_b[3])

                    region_q = (start_x, start_y, end_x, end_y)

                    # Tada. Is this how DP works?
                    # We're not breaking it into subproblems recursively.
                    # It should almost be like merge merge in that regard, lmao.
                    #

                    # print("Region q: {} ({})".format(region_q, score_seed_II(alphabet, scoring_matrix, s, t, region_q)))

                    # But what is up to 0 anyway?
                    # Hm. That I can't do this simply indicates

                    regions = [region_q] + [region for x, region in enumerate(initial_regions) if x != i and x != j]

                    # Which is wrong anyway lol

                    return join(regions)

        return initial_regions

    alignments = join(initial_regions)

    # print(alignments)

    # for alignment in alignments:
    #
    #     print(alignment, "({})".format(score_seed_II(alphabet, scoring_matrix, s, t, alignment)))

    # Those are atrocious scores lmao.
    # Oh well.

    start_x, start_y, end_x, end_y = max(alignments, key=lambda x: score_seed_II(alphabet, scoring_matrix, s, t, x))

    # print(start_x, start_y, end_x, end_y)

    # We don't trim it though. We just use the region.
    # The diagonal, that is.
    # Okay. I haven't done that before.

    diagonal = start_y - start_x

    if diagonal > 0:

        y = 0
        x = diagonal

        end = len(t) + 1


    else:

        y = abs(diagonal)
        x = 0

        end = len(s) + 1

    print(x, y)


    s = s[y:end]
    t = t[x:end]

    # Should truncate as well.

    print(s)
    print(t)

    D, _ = score_matrix(alphabet, scoring_matrix, s, t, band_width=8)

    print(D)

    # Empty sequence on heuristic2
    # And nothing for FASTA. Shame. I guess that we go through again and rectify any mistakes.
    # For now, I'm going to have a shower. I will get this finished.

    return traceback(alphabet, scoring_matrix, s, t, D, local=True)

    # Very respectable.


    # A list of regions: (start_x, start_y, end_x, end_y).

    # print(alignments)
    # print(len(alignments))

    # So we get 11.
    #

                # If joining them would not let the score decrease below a certain threshold?






    # DP works either bottom-up or top-down.
    # Top-down is simpler, and uses a dictionary to store results
    # Bottom up is....

    # Wait. So we keep passing through a list of the joined regions
    # We keep going until no regions can be be joined.
    # If we can join a region, we call the method again.

    # for diagonal, hotspot in initial_regions.
    # "Can be rapidly computing using DP...
    # Okay. So we join a region and then.... store the result.
    # Optimally, mind you.
    # We could banded DP for all of them.
    # Right? There's a max
    # Create an alignment of the trimmined initial regions using DP.
    # How would DP work here?
    # Only non-overlapping regions may be joined. Why? Cause we can't go 'backwards' in the alignment.
    # No, that's untrue. We can go left, but definitely not up.

    # Let's think.
    # 1. An alignment could start at any region.
    # 2. An alignment could end at any region.

    # But then we see why it's recursive / DP, right?
    # This is a key step.
    # I could just fudge it. But I don't want to.
    # This is the heuristic step that I keep going on about.




    # Given the locations of the initial regions, their respective scores, and a "joining penalty" (i.e. a gap penalty)
    # FASTA calculates an optimal alignment of initial regions as a combination of compatible regions.
    # By allowing the introduction of gaps... how is this not just DP?
    # Dijkstra's shortest-path algorithm is O(E + VlogV time). Very respectable.
    # Do I want the shortest path?
    # That assumes they all have to be joined up. Or I could do it in

    # Graphically, is this so hard?

    # But then it's TSP!

    # Vertices = initial regions

    # Vertex has a dictionary to say which vertices it is connected to.


    # Calculate initial edges, right?

    # Edges

    # How can I make this a graph?

    # Combine the diagonals into one big diagonal. Next determine if any of the initial regions from different
    # diagonals may be joined together to form an approximate alignment with gaps.

    # Okay.
    # Travelling horizontally or vertically corresponds to an indel. So we're gap matching.
    #


    # Okay. How hard can this be?

    #
    # Construct a directed weighted graph, where vertices are the runs, and edge weights represent gap penalties.
    # Find the best path through the graph: init_n.

    # print("Initial regions", initial_regions)


        #
        # if diagonal < 0:
        #
        #     start_x = 0
        #     start_y = abs(diagonal)
        #
        # else:
        #
        #     start_x = diagonal
        #     start_y = 0
        #
        # print("Hotspots", hotspots_)
        #
        # # start_x, start_y = hotspots_[0]
        # # Hm. Maybe not.
        # # end_x, end_y = hotspots_[-1][0] + ktup, hotspots_[-1][1] + ktup
        #
        # print("Diagonal, start_x, start_y")
        # print(diagonal, start_x, start_y)

        # Ah, fairs.
        # It sounds like we're running local alignment here.
        # Which would be fair enough.


        # After hashing... are RESCORED using a scoring matrix that allows conservative replacements.
        # Conservative replacements are characters matching that do not impose significant cost.

        # In the usual way. Optimally scoring sub-regions are identified.
        # I think I can DO that that.
        # That's basically the extension thing.
        # So we extend along the diagonal for each hotspot.
        # Right?
        #
        #So we have our diagonal runs. We find the optimal local alignment on these diagonals. They are called
        # initial regions. They'd mention using Banded SW if that were the case.
        # I think.
        # I can't quite decide.
        # We have to stick within the diagonal anyway.
        # So we look for the optimal local alignment within just this diagonal.
        # Interesting. That way we can identify a subregion with maximal score.
        # And that's all we want.
        # So how do I run local alignment just in that diagonal?
        # It corresponds to dynprog


    # diagonals = get_diagonals(alphabet, scoring_matrix, s, t, ktup, lookup)

    # hot_spot: matching ktup length substrings.
    # consecutive hot-spots are located along the diagonal
    # diagonal run:
    # a sequence of nearby hotspots on the same diagonal
    # init_1: the best scoring run
    # init_n: the best local alignment

    # print(lookup)

    # return "00000"

# scoring_matrix = np.array([
#     [2,-1,1,-1,-2],
#     [-1,2,-1,1,-2],
#     [1,-1,2,-1,-2],
#     [-1,1,-1,2,-2],
#     [-2,-2,-2,-2,0]
# ])
# alphabet = "ACGT_"
# s = "TACCAA"
# t = "ACTACG"

alphabet = "ABCD_"
scoring_matrix = np.array([
            [ 1,-5,-5,-5,-1],
            [-5, 1,-5,-5,-1],
            [-5,-5, 5,-5,-4],
            [-5,-5,-5, 6,-4],
            [-1,-1,-4,-4,-9]])
s = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
t = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"
#
# print("len(s)", len(s))
# print("len(t)", len(t))


# FASTA(alphabet, scoring_matrix, s, t)
def score_matrix(alphabet, scoring_matrix, s, t, band_width=None, threshold=0):

    print(s, t)

    size_y = len(s) + 1  # y-axis
    size_x = len(t) + 1  # x-axis

    D = np.zeros((size_y, size_x), dtype="int16")

    max_i = None  # The index of the cell containing the maximum score
    max_score = None  # The maximum score observed

    for y in range(size_y):  # y

        min_x = 0 if band_width is None else y - band_width  # If banded, the minimum x that we explore
        max_x = size_x if band_width is None else y + band_width + 1  # If banded, the maximum x that we explore.

        for x in range(min_x, max_x):

            if 0 < x < size_x and 0 < y:

                D[y, x] = max(
                    D[y - 1][x - 1] + cost(alphabet, scoring_matrix, s[y - 1], t[x - 1]),
                    D[y][x - 1] + cost(alphabet, scoring_matrix, t[x - 1], "_"),
                    D[y - 1][x] + cost(alphabet, scoring_matrix, s[y - 1], "_"),
                    threshold
                )

                if max_i is None or D[y, x] > max_score:

                    max_i = y, x
                    max_score = D[y, x]


    return D, max_i



hotspots = {'DDCDDC': [(0, 0)], 'DCDDCC': [(1, 1)], 'CDDCCC': [(2, 2)], 'DDCCCD': [(3, 3)], 'DCCCDC': [(4, 4)], 'CCCCDD': [(40, 11)], 'CCCDDD': [(41, 12)], 'CCDDDC': [(42, 13)], 'CDDDCD': [(43, 14)], 'DCDCDC': [(50, 61), (50, 63), (52, 61), (52, 63)], 'CDCDCD': [(51, 62), (51, 64), (53, 62), (53, 64)]}

regions = []

for word, coordinates in hotspots.items():

    for x, y in coordinates:

        regions.append((x, y, x + 6, y + 6))

# Right. This is absurd.
# Aha!

def score(region):

    start_x, start_y, end_x, end_y = region
    score = 0

    # print(s, t)

    for i in range(6):

        try:

            # print(start_x + i)
            # print(start_y + i)

            score += scoring_matrix[alphabet.index(s[start_y + i]), alphabet.index(t[start_x + i])]

        except IndexError:

            return score

    return score


# print(regions)

class Vertex:

    def __init__(self, region):

        self.start_x = region[0]
        self.start_y = region[1]
        self.end_x = region[2]
        self.end_y = region[3]

        self.weight = score(region)
        self.edges = []

        self.id = str(self.start_x) + "-" + str(self.start_y) + "-" + str(self.end_x) + "-" + str(self.end_y)


    def __str__(self):

        return "({}, {}) -> ({}, {}) [{}]".format(self.start_x, self.start_y, self.end_x, self.end_y, self.weight)

    def __eq__(self, other):

        return self.id == other.id

class Edge:

    def __init__(self, start, end, weight):

        self.start = start
        self.end = end
        self.weight = weight

    def __str__(self):

        return str(self.start) + " -> " + str(self.end) + " [" + str(self.weight) + "]"

vertices = [Vertex(region) for region in regions]

# We extend an edge from vertex u to vertex v.
# Makes sense.
print(len(vertices))    # 17.
# Nice. Kinda inefficient; we don't have to loop through them all.

for u in vertices:

    for v in vertices:

        if u != v:

            if u.end_y < v.start_y: # If the subalignment represented by v starts at a lower row than where u ends

                weight = v.start_y - u.end_y

                u.edges.append(Edge(u.id, v.id, weight))

                # An id.
                # Each edge needs a negative weight.

    for edge in u.edges:

        print(edge)


            # print(u, v)


#
#
# for vertex in vertices:
#
#     print(vertex)
#
#


# Constuct a directed weighted graph whose vertices are the subalignment found in the previous stage.
# The weight in each vertex is the score of the subalignment it represents.
# We extend an edge from vertex u to vertex v if the subalignment represented by v starters at a lower
# row than where the subalignment represents by v ends. We give the edge a negative weight which
# depends on the number of gaps that would be created by aligning according to subalignment v followed
# by subalignment u. Then, we find a maximum weight path in this graph. The selected alignment specifies
# a single local alignment between the two strings.
# Okay.






