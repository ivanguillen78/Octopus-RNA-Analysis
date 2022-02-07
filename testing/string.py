from Bio.Seq import Seq

teststring = "AGCTAGCTAGCTAGCT"

pos = 3
length = 7

seq = Seq(teststring[pos-(length//2):pos+(length//2+1)])
rev_comp = seq.reverse_complement()

print(seq+'\n')
print(rev_comp)