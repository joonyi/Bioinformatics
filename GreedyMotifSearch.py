# Copy your Score(Motifs), Count(Motifs), Profile(Motifs), and Consensus(Motifs) functions here.
def Score(Motifs):
    # Insert code here
    score = 0
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    k = len(Motifs[0])
    for j in range(k):
        symbols = ['A', 'C', 'G', 'T']
        symbols.remove(consensus[j])
        for sym in symbols:
            score += count[sym][j]
    return score

def Count(Motifs):
    k = len(Motifs[0])
    count = {'A':[0] * k, 'C':[0] * k, 'G':[0] * k, 'T':[0] * k} # initializing the count dictionary
    # your code here
    for j in range(k):
        for i in range(len(Motifs)):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count



def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
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

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {'A':[0] * k, 'C':[0] * k, 'G':[0] * k, 'T':[0] * k}
    # your code here
    for j in range(k):
        for i in range(t):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1 / t

    return profile

# Then copy your ProfileMostProbableKmer(Text, k, Profile) and Pr(Text, Profile) functions here.
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
def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

Dna = [
    ['G','G','C','G','T','T','C','A','G','G','C','A'],
    ['A','A','G','A','A','T','C','A','G','T','C','A'],
    ['C','A','A','G','G','A','G','T','T','C','G','C'],
    ['C','A','C','G','T','C','A','A','T','C','A','C'],
    ['C','A','A','T','A','A','T','A','T','T','C','G']
]
k, t = 3, 5
print(GreedyMotifSearch(Dna, k, t))
    