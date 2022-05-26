import re

from Bio.Seq import Seq

teststring = "TTTTCGTCTTTTGTACGTTTT"
#                 ((((    ).)))

position = 6


def return_longest_rev_comp_bulge(
    sequence, original_pos, left_index, maxLengthOfBulge=1, numBulges=1
):
    """
    Starts at leftmost index and returns longest reverse complement
    """
    substr_list = []
    rstring_list = []
    for index in range(left_index, original_pos + 1):
        left, right = index, index + 1
        bulgesLeft = numBulges
        substr = sequence[left:right]
        rstring = substr
        if right > len(sequence) - 1:
            substr_list.append(substr)
            break
        searchstring = sequence[:left] + sequence[right:]
        addition = sequence[right]
        while re.search(
            str(Seq(substr).reverse_complement()), searchstring
        ) and right < len(sequence):
            basecheck = check_bulge_base(
                substr, maxLengthOfBulge, sequence, right, searchstring, "right"
            )
            rev_comp_check = check_bulge_rev_comp(
                substr, maxLengthOfBulge, sequence, right, searchstring, "right"
            )
            if re.search(
                str(Seq(substr + addition).reverse_complement()), searchstring
            ):
                right += 1
                substr += addition
                rstring += addition
                if right == len(sequence):
                    break
                searchstring = sequence[:left] + sequence[right:]
                addition = sequence[right]
            elif basecheck[0] and bulgesLeft != 0:
                addition = sequence[right]
                rstring = rstring + addition
                substr = substr + ("." * basecheck[1]) + addition
                searchstring = sequence[:left] + sequence[right:]
                bulgesLeft -= 1
                if right < len(sequence) - 1:
                    right += 1
                    addition = sequence[right]
            elif rev_comp_check[0] and bulgesLeft != 0:
                right += rev_comp_check[1]
                addition = sequence[right]
                rstring = rstring + ("." * rev_comp_check[1]) + addition
                substr = substr + addition
                searchstring = sequence[:left] + sequence[right:]
                bulgesLeft -= 1
                if right < len(sequence) - 1:
                    right += 1
                    addition = sequence[right]
            else:
                break
        substr_list.append(substr)
        rstring_list.append(rstring)
    return max(enumerate(substr_list), key=lambda x: len(x[1])), max(
        enumerate(rstring_list), key=lambda x: len(x[1])
    )


def return_leftmost_index_bulge(sequence, pos, maxLengthOfBulge=1, numBulges=1):
    """
    Starts at given position and returns leftmost index for bulges
    """
    if pos == 0:
        return pos
    index, right = pos - 1, pos + 1
    substr = sequence[index + 1 : right]
    searchstring = sequence[: index + 1] + sequence[right:]
    bulgesLeft = numBulges
    addition = sequence[index]
    while re.search(str(Seq(substr).reverse_complement()), searchstring) and index >= 0:
        basecheck = check_bulge_base(
            substr, maxLengthOfBulge, sequence, index, searchstring, "left"
        )
        rev_comp_check = check_bulge_rev_comp(
            substr, maxLengthOfBulge, sequence, index, searchstring, "left"
        )
        if re.search(str(Seq(addition + substr).reverse_complement()), searchstring):
            index -= 1
            substr = addition + substr
            if index == 0:
                return 0
            searchstring = sequence[:index] + sequence[right:]
            addition = sequence[index]
        elif basecheck[0] and bulgesLeft != 0:
            addition = sequence[index]
            substr = addition + ("." * basecheck[1]) + substr
            searchstring = sequence[:index] + sequence[right:]
            bulgesLeft -= 1
            index -= 1
            addition = sequence[index]
        elif rev_comp_check[0] and bulgesLeft != 0:
            index -= rev_comp_check[1]
            addition = sequence[index]
            substr = addition + substr
            searchstring = sequence[:index] + sequence[right:]
            bulgesLeft -= 1
            index -= 1
            addition = sequence[index]
        else:
            break
    return index + 1


def check_bulge_base(substr, maxLength, sequence, index, searchstring, type):
    if type == "right":
        for i in range(1, maxLength + 1):
            if index + i < len(sequence):
                if re.search(
                    str(Seq(substr + ("." * i) + sequence[index]).reverse_complement()),
                    searchstring,
                ):
                    return (True, i)
    if type == "left":
        for i in range(1, maxLength + 1):
            if index - i >= 0:
                if re.search(
                    str(Seq(sequence[index] + ("." * i) + substr).reverse_complement()),
                    searchstring,
                ):
                    return (True, i)
    return (False, -1)


def check_bulge_rev_comp(substr, maxLength, sequence, index, searchstring, type):
    if type == "right":
        for i in range(1, maxLength + 1):
            if index + i < len(sequence):
                if re.search(
                    str(Seq(substr + sequence[index + i]).reverse_complement()),
                    searchstring,
                ):
                    return (True, i)
    if type == "left":
        for i in range(1, maxLength + 1):
            if index - i >= 0:
                if re.search(
                    str(Seq(sequence[index - i] + substr).reverse_complement()),
                    searchstring,
                ):
                    return (True, i)
    return (False, -1)


def create_secondary_structure_bulge(sequence, pos, maxLengthOfBulge=1, numBulges=1):
    """
    Returns all information needed to create output csv file (internal loops)
    """
    leftindex = return_leftmost_index_bulge(sequence, pos, maxLengthOfBulge, numBulges)
    structure, structure2 = return_longest_rev_comp_bulge(
        sequence, pos, leftindex, maxLengthOfBulge, numBulges
    )
    baselength = len(structure2[1])
    length = len(structure[1].replace(".", ""))
    base_string = structure2[1]
    rev_comp = str(Seq(structure[1]).reverse_complement())
    base_string_loc = [
        leftindex + structure2[0],
        leftindex + structure2[0] + baselength - 1,
    ]
    rev_comp_loc = [0, 0]
    loc = re.search(rev_comp, sequence)
    if loc:
        rev_comp_loc = [loc.start(), loc.end() - 1]
    # if (
    #     rev_comp_loc[0] >= base_string_loc[0] and rev_comp_loc[0] <= base_string_loc[1]
    # ) or (
    #     rev_comp_loc[1] >= base_string_loc[0] and rev_comp_loc[1] <= base_string_loc[1]
    # ):
    #     loc = re.search(rev_comp, sequence[pos + 1 :])
    #     if loc:
    #         rev_comp_loc = [loc.start(), loc.end() - 1]

    return length, base_string, base_string_loc, rev_comp, rev_comp_loc
