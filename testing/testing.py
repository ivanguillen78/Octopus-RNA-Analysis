import fastaparser
import csv
from Bio.Seq import Seq

edit_dict = {}

with open("../data/aes_profile.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if (row['orf'] not in edit_dict):
            edit_dict[row['orf']] = [int(row['pos'])]
        else:
            edit_dict[row['orf']].append(int(row['pos']))

def secondary_structure(id, pos_list, sequence, length):
    for pos in pos_list:
        rev_comp = sequence[pos-(length//2):pos+(length//2+1)].reverse_complement()
        if rev_comp in sequence:
            print(id)

with open("../data/swissprotORF.fasta") as fasta_file:
    parser = fastaparser.Reader(fasta_file)
    for seq in parser:
        if (seq.id in edit_dict):
            secondary_structure(seq.id, edit_dict[seq.id], Seq(seq.sequence_as_string()), 15) 
            
