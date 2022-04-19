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
    Takes in:
        - csv file containing ids and edits.
    Returns:
        - edit_dict: ids as keys and list of edits as values.
    """
    edit_dict = {}
    with open(csv_file, newline="", encoding="utf8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            edit_dict[row["orf"]] = int(row["pos"])
    csvfile.close()
    return edit_dict


def checkRight(sequence, pos):
    """
    Takes in:
        - sequence (string)
        - pos (int)
    Returns:
        - Largest reverse complement starting at pos and moving to the right (string)
    """
    lo, hi = pos, pos + 1
    right = pos
    if pos == len(sequence) - 1:
        return sequence[pos]
    while str(Seq(sequence[lo:hi]).reverse_complement()) in sequence:
        new_str = sequence[:pos] + sequence[right + 1 :]
        if str(Seq(sequence[lo : hi + 1]).reverse_complement()) in new_str and (
            hi < len(new_str)
        ):
            hi += 1
            right += 1
        if str(Seq(sequence[lo : hi + 1]).reverse_complement()) not in new_str or (
            hi >= len(new_str)
        ):
            return str(Seq(sequence[lo:hi]).reverse_complement())
    return ""


def checkLeft(sequence, pos):
    """
    Takes in:
        - sequence (string)
        - pos (int)
    Returns:
        - Largest reverse complement starting at pos and moving to the left (string)
    """
    lo, hi = pos, pos + 1
    left = pos
    rightSequence = checkRight(sequence, pos)
    if pos == 0:
        return sequence[pos]
    while (
        rightSequence[0 : len(rightSequence) - 1]
        + str(Seq(sequence[lo:hi]).reverse_complement())
        in sequence
    ):
        new_str = sequence[:left] + sequence[pos + 1 :]
        if (
            rightSequence[0 : len(rightSequence) - 1]
            + str(Seq(sequence[lo - 1 : hi]).reverse_complement())
            in new_str
            and lo > -1
        ):
            lo -= 1
            left -= 1
        if rightSequence[0 : len(rightSequence) - 1] + str(
            Seq(sequence[lo - 1 : hi]).reverse_complement()
        ) not in new_str or (lo <= 0):
            return str(Seq(sequence[lo:hi]).reverse_complement())
    return ""


def secondary_structure(sequence, pos):
    """
    Takes in:
        - pos: edit site
        - sequence: genetic sequence corresponding to pos_list
    Returns:
        - length of reverse complement
        - sequence around edit site + location
        - reverse complement + location
    """
    right = checkRight(sequence, pos)
    left = checkLeft(sequence, pos)
    length = len(right + left) - 1
    rev_comp = right[0 : len(right) - 1] + left
    base_string = str(Seq(rev_comp).reverse_complement())
    rev_comp_loc_start = sequence.find(rev_comp)
    rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]
    base_string_loc = [pos - len(left) + 1, pos + len(right) - 1]
    return length, base_string, base_string_loc, rev_comp, rev_comp_loc


def find_secondary_structures(edit_dict, fasta_file):
    """
    Takes in:
        - edit_dict: ids as keys and list of edits as values
        - fasta_file: file containing all genetic sequences
    Returns:
        - score_dict: ids as keys and scores (length for now) as values
    """
    output_list = []
    with open(fasta_file, encoding="utf8") as fastafile:
        size = int(len(fastafile.readlines()) / 2)
        parser = fastaparser.Reader(fastafile)
        with alive_bar(size) as bar:
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
    Takes in:
        - file_name: name of output csv file to be created
        - score_dict: ids as keys and scores (length for now) as values
    Creates:
        - output csv file
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
