# Insert your Pr(text, profile) function here from Motifs.py.
def Pr(Text, Profile):
    # insert your code here
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]

    return p


# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbableKmer(text, k, profile):
    p = -1
    kmer = ''
    for i in range(len(text) - k + 1):
        q = Pr(text[i:i + k], profile)
        if p < q:
            p = q
            kmer = text[i:i + k]

    return kmer

text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
k = 5
profile = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]
}
print(ProfileMostProbableKmer(text, k, profile))

# Quiz
profile = {
    'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
}
Text = 'GAGCTA'
print(Pr(Text, profile))