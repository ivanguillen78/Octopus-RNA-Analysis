import fastaparser
import csv
from Bio.Seq import Seq

edit_dict = {}
score_dict = {}

with open("../data/aes_profile.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if (row['orf'] not in edit_dict):
            edit_dict[row['orf']] = [int(row['pos'])]
        else:
            edit_dict[row['orf']].append(int(row['pos']))

def secondary_structure(pos_list, sequence):
    len_list = []
    for pos in pos_list:
        lo = pos
        hi = pos + 1
        for i in range(len(sequence)):
           rev_comp = str(Seq(sequence[lo:hi]).reverse_complement())
           if (rev_comp not in sequence):
               break
           if ((rev_comp in sequence) and (rev_comp not in sequence[lo:hi])):
               lo -= 1
               hi += 1
        len_list.append(len(rev_comp))
    return len_list

with open("../data/swissprotORF.fasta") as fasta_file:
    parser = fastaparser.Reader(fasta_file)
    for seq in parser:
        if (seq.id in edit_dict):
            score_dict[seq.id] = secondary_structure(edit_dict[seq.id], seq.sequence_as_string())

with open('score_list.csv', 'w', newline='') as csvfile:
    fieldnames = ['id', 'score']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for id, score in score_dict.items():
        writer.writerow({'id': id, 'score': score})