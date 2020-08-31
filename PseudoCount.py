# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
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


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
# HINT:   You need to use CountWithPseudocounts as a subroutine of ProfileWithPseudocounts
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}
    for j in range(k):
        for i in range(t):
            count[Motifs[i][j]][j] += 1
    return count


Motifs = ['AACGTA',
          'CCCGTT',
          'CACCTT',
          'GGATTA',
          'TTCCGG'
          ]
print(ProfileWithPseudocounts(Motifs))
