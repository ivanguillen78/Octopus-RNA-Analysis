import fastaparser
import csv
import pandas as pd

seq_list = []
edit_list = []

with open("../data/swissprotORF.fasta") as fasta_file:
    parser = fastaparser.Reader(fasta_file)
    for seq in parser:
        seq_list.append({'id': seq.id, 'seq': seq.sequence_as_string()})


with open("../data/aes_profile.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        edit_list.append({'id': row['orf'], 'pos': [int(row['pos'])]})

sortedList = sorted(edit_list, key=lambda x: x['id'])

for x in sortedList:
    for y in seq_list:
        if x['id'] == y['id']:
            x['seq'] = y['seq']

for i in range(20):
    print(sortedList[i])
