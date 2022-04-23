#!/usr/bin/python
"""
Identifying secondary structures given
    - FASTA file of genetic sequences
    - CSV file with potential edit sites
"""
from alive_progress import alive_bar
import csv
import fastaparser
from Bio.Seq import Seq


def create_edit_dict(csv_file):
    """
    Returns a dictionary with ids as keys and list of edits as values.
    """
    edit_dict = {}
    with open(csv_file, newline="", encoding="utf8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            edit_dict[row["orf"]] = int(row["pos"])
    csvfile.close()
    return edit_dict


def return_leftmost_index(sequence, pos):
    """
    Starts at given position and returns leftmost index.
    """
    if (pos == 0):
        return pos
    left, right = pos, pos + 1
    substr = str(Seq(sequence[left:right]).reverse_complement())
    searchstring = sequence[:left] + sequence[right:]
    while substr in searchstring and left >= 0:
        left -= 1
        substr = str(Seq(sequence[left:right]).reverse_complement())
        searchstring = sequence[:left] + sequence[right:]
    return left + 1

def return_longest_rev_comp(sequence, original_pos, left_index):
    """
    Starts at leftmost index and returns longest reverse complement.
    """
    substr_list = []
    for index in range(left_index, original_pos + 1):
        left, right = index, index + 1
        substr = str(Seq(sequence[left:right]).reverse_complement())
        searchstring = sequence[:left] + sequence[right:]
        while substr in searchstring and right < len(sequence):
            right += 1
            substr = str(Seq(sequence[left:right]).reverse_complement())
            searchstring = sequence[:left] + sequence[right:]
            if str(Seq(sequence[left:right+1]).reverse_complement()) not in searchstring:
                break
        substr_list.append(str(Seq(sequence[left:right]).reverse_complement()))
    return max(enumerate(substr_list), key=lambda x: len(x[1]))


def secondary_structure(sequence, pos):
    """
    Returns all information needed to create output csv file.
    """
    messedup = False
    leftindex = return_leftmost_index(sequence, pos)
    rev_comp = return_longest_rev_comp(sequence, pos, leftindex)
    length = len(rev_comp[1])
    base_string = sequence[leftindex+rev_comp[0]:leftindex + rev_comp[0] + length]
    base_string_loc = [leftindex+rev_comp[0], leftindex+rev_comp[0]+length-1]
    rev_comp_loc_start = sequence.find(rev_comp[1])
    rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    if (
         rev_comp_loc[0] > base_string_loc[0] and rev_comp_loc[0] < base_string_loc[1]
     ) or (
         rev_comp_loc[1] > base_string_loc[0] and rev_comp_loc[1] < base_string_loc[1]
     ):
         messedup = True

    if messedup:
        rev_comp_loc_start = sequence.find(rev_comp[1], base_string_loc[1]+1, len(sequence)-1)
        rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    return length, base_string, base_string_loc, rev_comp[1], rev_comp_loc


def find_secondary_structures(edit_dict, fasta_file):
    """
    Iterates through edit dictionary and finds longest reverse complement
    for each edit.
    """
    output_list = []
    with open(fasta_file, encoding="utf8") as fastafile:
        size = int(len(fastafile.readlines()) / 2)
        parser = fastaparser.Reader(fastafile)
        with alive_bar(size, bar='smooth', spinner='fish2') as bar:
            for seq in parser:
                bar()
                if seq.id in edit_dict:
                    length, base, base_loc, rev, rev_loc = secondary_structure(
                        seq.sequence_as_string(), edit_dict[seq.id]
                    )
                    output_list.append(
                        {
                            "id": seq.id,
                            "position": edit_dict[seq.id],
                            "length": length,
                            "base": base,
                            "base_location": base_loc,
                            "rev_comp": rev,
                            "rev_comp_location": rev_loc,
                        }
                    )
    fastafile.close()
    return output_list


def create_output_csv(file_name, score_list, minLength=5):
    """
    Creates output csv file.
    """
    with open(file_name, "w", newline="", encoding="utf8") as csvfile:
        fieldnames = [
            "id",
            "position",
            "length",
            "base_string",
            "base_string_loc",
            "rev_comp",
            "rev_comp_loc",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for item in score_list:
            if item["length"] >= minLength:
                writer.writerow(
                    {
                        "id": item["id"],
                        "position": item["position"],
                        "length": item["length"],
                        "base_string": item["base"],
                        "base_string_loc": item["base_location"],
                        "rev_comp": item["rev_comp"],
                        "rev_comp_loc": item["rev_comp_location"],
                    }
                )
    csvfile.close()
