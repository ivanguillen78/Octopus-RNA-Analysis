from Bio.Seq import Seq
#TAGCTAC
test_seq = 'ATCGATCGATCGATCGATCGGCACGTTCGATCCGTACTTAGAGCTAGATCCCCGCTGAGCTGGATCGTTTTTTGATCTAGCTCTAAGTACGGATCGAAGACATGAGGGGGATAGAATAAGACATAGACGATCGATCGATCGAT'
#                ^edit site
pos = 36

#reverse = test_seq[pos]

lo = pos
hi = pos + 1

#print(test_seq[lo:hi]+'\n')
# first line inside while loop is GCT
# converts GCT to string
# decrements lo
# incrememnts hi
# prints
# GCT is in 
# reverse = str(Seq(test_seq[lo:hi]).reverse_complement())
for i in range (len(test_seq)):
    reverse = str(Seq(test_seq[lo:hi]).reverse_complement())
    if (reverse not in test_seq):
        break;
    # if ((str(Seq(test_seq[lo-1:hi]).reverse_complement()) in test_seq) and (reverse not in test_seq[lo:hi])):
    #     lo -= 1
    # if ((str(Seq(test_seq[lo:hi+1]).reverse_complement()) in test_seq) and (reverse not in test_seq[lo:hi])):
    #     hi += 1
    if ((reverse in test_seq) and (reverse not in test_seq[lo:hi])):
        lo -= 1
        hi += 1
    print(len(reverse))

#print(Seq('TAGCTA').reverse_complement())
# ACTGATACTAGCT