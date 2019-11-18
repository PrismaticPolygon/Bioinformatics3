# An algorithm for aligning two sequences within a diagonal band that requires only O(NW) time and O(N)
# space, where N is the length of the shorter of the two sequences and W is the width of the band.

# Can be used for EITHER local or global alignment scores.

# We have is O(MW)-space, score-only algorithm for local alignment within a band

# Linear-space, score-only algorithm for local alignment within a band.



# To construct an alignment in linear space, the score-only algorithm determines the point where
# an optimal path crosses the middle row, leaving two subproblems to be solved.

# Local alignments are produced by finding the beginning and end of the best local alignment in the band
# and then applying the global local alignment in between those points.

# Aligning a_i with b_j is permitted
# L = lower diagonal limit, i.e. L <= j - i <= U