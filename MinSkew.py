# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = []  # output variable
    # your code here
    skew = SkewArray(Genome)
    min_ = min(skew)
    for i, gen in enumerate(skew):
        if gen == min_:
            positions.append(i)
    return positions


# Input:  A String Genome
# Output: SkewArray(Genome)
# HINT:   This code should be taken from the last Code Challenge.
def SkewArray(Genome):
    skew = [0] * (len(Genome) + 1)
    for i in range(len(Genome)):
        if Genome[i] == "A" or Genome[i] == "T":
            skew[i + 1] = skew[i]
        elif Genome[i] == "G":
            skew[i + 1] = skew[i] + 1
        elif Genome[i] == "C":
            skew[i + 1] = skew[i] - 1

    return skew

def HammingDistance(p, q):
    dist = 0
    for i, j in zip(p, q):
        if i != j:
            dist += 1

    return dist

Genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
Genome = "GATACACTTCCCGAGTAGGTACTG"
print(MinimumSkew(Genome))
print(HammingDistance("TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC",
                      "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"))