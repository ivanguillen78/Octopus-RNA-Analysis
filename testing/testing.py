'''Identifying secondary structures given
 - FASTA file of genetic sequences
 - CSV file with potential edit sites
'''
import csv
import fastaparser
from Bio.Seq import Seq

edit_dict = {}
score_dict = {}

with open("../data/aes_profile.csv", newline='', encoding="utf8") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['orf'] not in edit_dict:
            edit_dict[row['orf']] = [int(row['pos'])]
        else:
            edit_dict[row['orf']].append(int(row['pos']))

def secondary_structure(pos_list, sequence):
    '''
    Takes in:
        - pos_list: list of potential edit positions
        - sequence: genetic sequence corresponding to pos_list
    Returns:
        - len_list: list of secondary structure lengths for pos_list
    '''
    len_list = []
    for pos in pos_list:
        lo_ = pos
        hi_ = pos + 1
        for i in range(len(sequence)): # pylint: disable=unused-variable
            rev_comp = str(Seq(sequence[lo_:hi_]).reverse_complement())
            if rev_comp not in sequence:
                break
            if ((rev_comp in sequence) and (rev_comp not in sequence[lo_:hi_])):
                lo_ -= 1
                hi_ += 1
        len_list.append(len(rev_comp))
    return len_list

with open("../data/swissprotORF.fasta", encoding="utf8") as fasta_file:
    parser = fastaparser.Reader(fasta_file)
    for seq in parser:
        if seq.id in edit_dict:
            score_dict[seq.id] = secondary_structure(edit_dict[seq.id], seq.sequence_as_string())

with open('score_list.csv', 'w', newline='', encoding="utf8") as csvfile:
    fieldnames = ['id', 'score']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for id_, score in score_dict.items():
        writer.writerow({'id': id_, 'score': score})
