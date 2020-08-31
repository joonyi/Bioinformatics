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

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]

    return p

def ProfileMostProbableKmer(text, k, profile):
    p = -1
    kmer = ''
    for i in range(len(text) - k + 1):
        q = Pr(text[i:i + k], profile)
        if p < q:
            p = q
            kmer = text[i:i + k]

    return kmer

def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = len(Profile['A'])
    for i in range(0, t):
        motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))

    return motifs

Dna = [
    'TGACGTTC',
    'TAAGAGTT',
    'GGACGAAA',
    'CTGTTCGC'
]
rand_motif = [
    'TGA',
    'GTT',
    'GAA',
    'TGT'
]
k, t = 3, len(Dna)

print(Motifs(ProfileWithPseudocounts(rand_motif), Dna))
