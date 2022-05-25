from Bio.Seq import Seq


def check_two_strings(substr, string1, string2):
    if substr in string1 or substr in string2:
        return True
    return False


def return_leftmost_index_hairpin(sequence, pos):
    """
    Starts at given position and returns leftmost index.
    """
    if pos == 0:
        return pos
    left, right = pos, pos + 1
    substr = sequence[left:right]
    while (
        check_two_strings(
            str(Seq(substr).reverse_complement()), sequence[:left], sequence[right:]
        )
    ) and left >= 0:
        if (
            not (
                check_two_strings(
                    str(Seq(sequence[left - 1 : right]).reverse_complement()),
                    sequence[:left],
                    sequence[right:],
                )
            )
            or left == 0
        ):
            break
        left -= 1
        substr = sequence[left:right]
    return left


def return_rightmost_index_hairpin(sequence, left, right):
    """
    Starts at given position and returns rightmost index.
    """
    if right == len(sequence) - 1:
        return right
    substr = sequence[left:right]
    while (
        check_two_strings(
            str(Seq(substr).reverse_complement()), sequence[:left], sequence[right:]
        )
    ) and right < len(sequence):
        if (
            not (
                check_two_strings(
                    str(Seq(sequence[left : right + 1]).reverse_complement()),
                    sequence[:left],
                    sequence[right:],
                )
            )
            or right == len(sequence) - 1
        ):
            break
        right += 1
        substr = sequence[left:right]
    return right


def return_longest_hairpin(sequence, position, left_index):
    left = left_index
    right = return_rightmost_index_hairpin(sequence, left_index, position + 1)
    longest = sequence[left_index:right]
    for index in range(left_index + 1, position + 1):
        if right + 1 >= len(sequence):
            break
        if not check_two_strings(
            str(Seq(sequence[index : right + 1]).reverse_complement()),
            sequence[:index],
            sequence[right + 1 :],
        ):
            continue
        right = return_rightmost_index_hairpin(sequence, index, right + 1)
        if len(sequence[index:right]) > len(longest):
            longest = sequence[index:right]
            left = index
    return [str(Seq(longest).reverse_complement()), left]


def create_secondary_structure_hairpin(sequence, pos):
    """
    Returns all information needed to create output csv file.
    """
    leftindex = return_leftmost_index_hairpin(sequence, pos)
    structure = return_longest_hairpin(sequence, pos, leftindex)
    rev_comp = structure[0]
    length = len(rev_comp)
    base_string = sequence[structure[1] : structure[1] + length]
    base_string_loc = [structure[1], structure[1] + length - 1]
    rev_comp_loc_start = sequence.find(rev_comp)
    rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    if (
        rev_comp_loc[0] >= base_string_loc[0] and rev_comp_loc[0] <= base_string_loc[1]
    ) or (
        rev_comp_loc[1] >= base_string_loc[0] and rev_comp_loc[1] <= base_string_loc[1]
    ):
        rev_comp_loc_start = sequence.find(
            rev_comp, base_string_loc[1] + 1, len(sequence) - 1
        )
        rev_comp_loc = [rev_comp_loc_start, rev_comp_loc_start + length - 1]

    return length, base_string, base_string_loc, rev_comp, rev_comp_loc
