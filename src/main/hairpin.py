from Bio.Seq import Seq


def return_leftmost_index_hairpin(sequence, pos):
    """
    Starts at given position and returns leftmost index.
    """
    if pos == 0:
        return pos
    left, right = pos, pos + 1
    substr = sequence[left:right]
    searchstring = sequence[:left] + sequence[right:]
    while str(Seq(substr).reverse_complement()) in searchstring and left >= 0:
        if (
            str(Seq(sequence[left - 1 : right]).reverse_complement())
            not in searchstring
            or left == 0
        ):
            break
        left -= 1
        substr = sequence[left:right]
        searchstring = sequence[:left] + sequence[right:]
    return left


def return_longest_rev_comp_hairpin(sequence, position, left_index):
    """
    Starts at leftmost index and returns longest reverse complement.
    """
    substr_list = []
    for index in range(left_index, position + 1):
        left, right = index, index + 1
        substr = str(Seq(sequence[left:right]).reverse_complement())
        searchstring = sequence[:left] + sequence[right:]
        while substr in searchstring and right < len(sequence):
            right += 1
            substr = str(Seq(sequence[left:right]).reverse_complement())
            searchstring = sequence[:left] + sequence[right:]
            if (
                str(Seq(sequence[left : right + 1]).reverse_complement())
                not in searchstring
            ):
                break
        substr_list.append(str(Seq(sequence[left:right]).reverse_complement()))
    return max(enumerate(substr_list), key=lambda x: len(x[1]))


def create_secondary_structure_hairpin(sequence, pos):
    """
    Returns all information needed to create output csv file.
    """
    leftindex = return_leftmost_index_hairpin(sequence, pos)
    rev_comp = return_longest_rev_comp_hairpin(sequence, pos, leftindex)
    length = len(rev_comp[1])
    base_string = sequence[leftindex + rev_comp[0] : leftindex + rev_comp[0] + length]
    base_string_loc = [leftindex + rev_comp[0], leftindex + rev_comp[0] + length - 1]
    rev_comp_loc_start = sequence.find(rev_comp[1])
    rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    if (
        rev_comp_loc[0] >= base_string_loc[0] and rev_comp_loc[0] <= base_string_loc[1]
    ) or (
        rev_comp_loc[1] >= base_string_loc[0] and rev_comp_loc[1] <= base_string_loc[1]
    ):
        rev_comp_loc_start = sequence.find(
            rev_comp[1], base_string_loc[1] + 1, len(sequence) - 1
        )
        rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    return length, base_string, base_string_loc, rev_comp[1], rev_comp_loc
