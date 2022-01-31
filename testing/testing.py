import fastaparser
import csv
import pandas as pd

seq_list = []
edit_list = []

# with open("../data/swissprotORF.fasta") as fasta_file:
#     parser = fastaparser.Reader(fasta_file)
#     for seq in parser:
#         seq_list.append({'id': seq.id, 'seq': seq.sequence_as_string()})


with open("test.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        edit_list.append({'id': row['orf'], 'pos': [int(row['pos'])]})

df = pd.DataFrame(edit_list)
print(type(df['id'][0]))
# for x in edit_list:
#     for y in seq_list:
#         if x['id'] == y['id']:
#             x['seq'] = y['seq']

# for i in range(len(edit_list)):
#     if (edit_list[i]['id'] == edit_list[i-1]['id']):
#         edit_list[i]['flag'] = 'duplicate'
#         #del edit_list[i-1]


# for i in range(len(edit_list)):
#     print(edit_list[i])
