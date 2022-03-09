from Bio.Seq import Seq
seq1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
#       (((((^((((((((............))))))))^)))))...........
pos1 = 5

seq3 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
#       ...(((..........................................))^
pos3 = 0
    
def checkRight(sequence, pos):
    lo, hi = pos, pos+1
    right = pos
    if (pos == len(sequence) - 1):
        return 1
    while (str(Seq(sequence[lo:hi]).reverse_complement()) in sequence):
        new_str = sequence[:pos] + sequence[right + 1:]
        if (str(Seq(sequence[lo : hi + 1]).reverse_complement()) in new_str and (hi < len(new_str))):
            hi += 1
            right += 1
        if (str(Seq(sequence[lo : hi + 1]).reverse_complement()) not in new_str or (hi >= len(new_str))):
            return (len(str(Seq(sequence[lo:hi]).reverse_complement())))
    #return (hi-lo)

def checkLeft(sequence, pos):
    lo, hi = pos, pos+1
    left = pos
    if (pos == 0):
        return 1
    while (str(Seq(sequence[lo:hi]).reverse_complement()) in sequence):
        new_str = sequence[:left] + sequence[pos + 1:]
        if (str(Seq(sequence[lo - 1: hi]).reverse_complement()) in new_str and lo > -1):
            lo -= 1
            left -= 1
        if (str(Seq(sequence[lo-1: hi]).reverse_complement()) not in new_str or (lo <= 0)):
            return (len(str(Seq(sequence[lo:hi]).reverse_complement())))
    #return (hi-lo)
# print(pos3)
# print(pos3+1)
print(checkLeft(seq1, pos1))
print(checkRight(seq1, pos1))
#