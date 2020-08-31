# Copy your Consensus(Motifs) function here.
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

# Copy your Count(Motifs) function here.
def Count(Motifs):
    k = len(Motifs[0])
    count = {'A': [0] * k, 'C': [0] * k, 'G': [0] * k, 'T': [0] * k}  # initializing the count dictionary
    # your code here
    for j in range(k):
        for i in range(len(Motifs)):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
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


Motifs = [['A','A','C','G','T','A'],
          ['C','C','C','G','T','T'],
          ['C','A','C','C','T','T'],
          ['G','G','A','T','T','A'],
          ['T','T','C','C','G','G'],
          ]
print(Score(Motifs))