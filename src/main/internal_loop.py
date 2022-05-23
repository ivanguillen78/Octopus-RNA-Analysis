#!/usr/bin/python
"""
Identifying secondary structures given
    - FASTA file of genetic sequences
    - CSV file with potential edit sites
"""
from Bio.Seq import Seq
import re


def return_longest_rev_comp_loop(
    sequence, original_pos, left_index, maxLengthOfLoop=1, numLoops=1
):
    """
    Starts at leftmost index and returns longest reverse complement
    """
    substr_list = []
    for index in range(left_index, original_pos + 1):
        left, right = index, index + 1
        loopsLeft = numLoops
        substr = sequence[left:right]
        if right > len(sequence) - 1:
            substr_list.append(substr)
            break
        searchstring = sequence[:left] + sequence[right:]
        addition = sequence[right]
        while re.search(
            str(Seq(substr).reverse_complement()), searchstring
        ) and right < len(sequence):
            checking = checkInternalLoop(
                substr, maxLengthOfLoop, sequence, right, searchstring, "right"
            )
            if re.search(
                str(Seq(substr + addition).reverse_complement()), searchstring
            ):
                right += 1
                substr += addition
                if right == len(sequence):
                    break
                searchstring = sequence[:left] + sequence[right:]
                addition = sequence[right]
            elif checking[0] and loopsLeft != 0:
                right += checking[1]
                addition = sequence[right]
                substr = substr + ("." * checking[1]) + addition
                searchstring = sequence[:left] + sequence[right + checking[1] :]
                loopsLeft -= 1
                if right < len(sequence) - 1:
                    right += 1
                    addition = sequence[right]
            else:
                break
        substr_list.append(substr)
    return max(enumerate(substr_list), key=lambda x: len(x[1]))


def return_leftmost_index_loop(sequence, pos, maxLengthOfLoop=1, numLoops=1):
    """
    Starts at given position and returns leftmost index for internal loops
    """
    if pos == 0:
        return pos
    index, right = pos - 1, pos + 1
    substr = sequence[index + 1 : right]
    searchstring = sequence[: index + 1] + sequence[right:]
    loopsLeft = numLoops
    addition = sequence[index]
    while re.search(str(Seq(substr).reverse_complement()), searchstring) and index >= 0:
        checking = checkInternalLoop(
            substr, maxLengthOfLoop, sequence, index, searchstring, "left"
        )
        if re.search(str(Seq(addition + substr).reverse_complement()), searchstring):
            index -= 1
            substr = addition + substr
            if index == 0:
                return 0
            searchstring = sequence[:index] + sequence[right:]
            addition = sequence[index]
        elif checking[0] and loopsLeft != 0:
            index -= checking[1]
            addition = sequence[index]
            substr = addition + ("." * checking[1]) + substr
            searchstring = sequence[:index] + sequence[right:]
            loopsLeft -= 1
            index -= 1
            addition = sequence[index]
        else:
            break
    return index + 1


def checkInternalLoop(substr, maxLength, sequence, index, searchstring, type):
    if type == "right":
        for i in range(1, maxLength + 1):
            if index + i < len(sequence):
                if re.search(str(Seq(substr + ("." * i) + sequence[index + i]).reverse_complement()), searchstring):
                    return (True, i)
    if type == "left":
        for i in range(1, maxLength + 1):
            if index - i >= 0:
                if re.search(str(Seq(sequence[index - i] + ("." * i) + substr).reverse_complement()), searchstring):
                    return (True, i)
    return (False, -1)


def create_secondary_structure_loop(sequence, pos):
    """
    Returns all information needed to create output csv file (internal loops)
    """
    leftindex = return_leftmost_index_loop(sequence, pos)
    structure = return_longest_rev_comp_loop(sequence, pos, leftindex)
    rev_comp = str(Seq(structure[1]).reverse_complement())
    length = len(rev_comp)
    base_string = structure[1]
    base_string_loc = [leftindex + structure[0], leftindex + structure[0] + length - 1]
    rev_comp_loc = [0, 0]
    loc = re.search(rev_comp, sequence)
    if loc:
        rev_comp_loc = [loc.start(), loc.end()]
    # if (
    #     rev_comp_loc[0] >= base_string_loc[0] and rev_comp_loc[0] <= base_string_loc[1]
    # ) or (
    #     rev_comp_loc[1] >= base_string_loc[0] and rev_comp_loc[1] <= base_string_loc[1]
    # ):
    #     rev_comp_loc_start = sequence.find(
    #         rev_comp[1], base_string_loc[1] + 1, len(sequence) - 1
    #     )
    #     rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    return length, base_string, base_string_loc, rev_comp, rev_comp_loc
