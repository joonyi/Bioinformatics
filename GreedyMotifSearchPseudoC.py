def Score(Motifs):
    # Insert code here
    score = 0
    count = CountWithPseudocounts(Motifs)
    consensus = Consensus(Motifs)
    k = len(Motifs[0])
    for j in range(k):
        symbols = ['A', 'C', 'G', 'T']
        symbols.remove(consensus[j])
        for sym in symbols:
            score += count[sym][j]
    return score


def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}
    for j in range(k):
        for i in range(t):
            count[Motifs[i][j]][j] += 1
    return count


def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {'A': [0] * k, 'C': [0] * k, 'G': [0] * k, 'T': [0] * k}  # output variable
    t_norm = t + 4
    count = CountWithPseudocounts(Motifs)
    for sym in count:
        for j in range(k):
            profile[sym][j] = count[sym][j] / t_norm

    return profile


def ProfileMostProbableKmer(text, k, profile):
    p = -1
    kmer = ''
    for i in range(len(text) - k + 1):
        q = Pr(text[i:i + k], profile)
        if p < q:
            p = q
            kmer = text[i:i + k]

    return kmer


def Pr(Text, Profile):
    # insert your code here
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]

    return p


# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs





Dna = [
    ['G', 'G', 'C', 'G', 'T', 'T', 'C', 'A', 'G', 'G', 'C', 'A'],
    ['A', 'A', 'G', 'A', 'A', 'T', 'C', 'A', 'G', 'T', 'C', 'A'],
    ['C', 'A', 'A', 'G', 'G', 'A', 'G', 'T', 'T', 'C', 'G', 'C'],
    ['C', 'A', 'C', 'G', 'T', 'C', 'A', 'A', 'T', 'C', 'A', 'C'],
    ['C', 'A', 'A', 'T', 'A', 'A', 'T', 'A', 'T', 'T', 'C', 'G']
]
k, t = 3, 5
print(GreedyMotifSearchWithPseudocounts(Dna, k, t))
