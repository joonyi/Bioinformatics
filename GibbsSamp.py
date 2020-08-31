import random


def Score(Motifs):
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

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}
    for j in range(k):
        for i in range(t):
            count[Motifs[i][j]][j] += 1
    return count

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

def Normalize(Probabilities):
    # your code here
    total = sum(Probabilities.values())
    for key, val in Probabilities.items():
        Probabilities[key] = val / total

    return Probabilities

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]

    return p

def WeightedDie(Probabilities):
    kmer = '' # output variable
    n = random.uniform(0, 1)
    for kmer in Probabilities:
        n -= Probabilities[kmer]
        if n <= 0:
            return kmer

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def RandomMotifs(Dna, k, t):
    t = len(Dna)
    l = len(Dna[0])
    random_motifs = []
    for i in range(t):
        r = random.randint(1, l-k)
        random_motifs.append(Dna[i][r:r+k])
    return random_motifs


# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs.copy()   # output variable
    for j in range(N):
        exclude = random.randint(0, t - 1)
        profile = ProfileWithPseudocounts([BestMotifs[i] for i in range(len(BestMotifs)) if i != exclude])
        Motifs[exclude] = ProfileGeneratedString(Dna[exclude], profile, k)

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


#############################################
p = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
p = Normalize(p)
print(WeightedDie(p))

############################################
Text = 'AAACCCAAACCC'
profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
k = 2
print(ProfileGeneratedString(Text, profile, k))

###########################################
k, t, N = 8, 5, 100
Dna = [
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA',
]

k, t, N = 4, 5, 100
Dna = [
    'TTACCTTAAC',
    'GATGTCTGTC',
    'CCGGCGTTAG',
    'CACTAACGAG',
    'CGTCAGAGGT'
]
print(GibbsSampler(Dna, k, t, N))