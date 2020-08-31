# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2], symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1

    return array

# Input:  Strings Text and Pattern
# Output: The number of times Pattern appears in Text
# HINT:   This code should be identical to when you last implemented PatternCount
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Pattern) - len(Text) + 1):
        if Pattern[i:i + len(Text)] == Text:
            count += 1
    return count


with open('e_coli.txt') as file:
    e_coli = file.read()

# genome, symbol = e_coli, 'C'
genome, symbol = 'AAAAGGGG', 'A'
array = FasterSymbolArray(genome, symbol)

import matplotlib.pyplot as plt
plt.plot(*zip(*sorted(array.items())))
plt.show()