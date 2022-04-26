import re
from Bio.Seq import Seq

def return_leftmost_index(sequence, pos):
    """
    Starts at given position and returns leftmost index.
    """
    if (pos == 0):
        return pos
    left, right = pos, pos + 1
    substr = sequence[left:right]
    searchstring = sequence[:left] + sequence[right:]
    while substr in searchstring and left >= 0:
        if sequence[left-1:right] not in searchstring:
            break
        left -= 1
        substr = sequence[left:right]
        searchstring = sequence[:left] + sequence[right:]
    return left

def return_longest_rev_comp(sequence, original_pos, left_index):
    """
    Starts at leftmost index and returns longest reverse complement.
    """
    substr_list = []
    for index in range(left_index, original_pos + 1):
        left, right = index, index + 1
        substr = sequence[left:right]
        searchstring = sequence[:left] + sequence[right:]
        while substr in searchstring and right < len(sequence):
            if sequence[left:right+1] not in searchstring:
                break
            right += 1
            substr = sequence[left:right]
            searchstring = sequence[:left] + sequence[right:]
        substr_list.append(sequence[left:right])
    return max(enumerate(substr_list), key=lambda x: len(x[1]))

def return_longest_rev_comp_bulge(sequence, original_pos, left_index, maxLengthOfBulge, numBulges):
    """
    Starts at leftmost index and returns longest reverse complement
    """
    substr_list = []
    for index in range(left_index, original_pos + 1):
        left, right = index, index + 1
        bulgesLeft = numBulges
        substr = sequence[left:right]
        if right > len(sequence) - 1:
            substr_list.append(substr)
            break
        searchstring = sequence[:left] + sequence[right:]
        addition = sequence[right]
        while re.search(substr, searchstring) and right < len(sequence):
            if re.search(substr + addition, searchstring):
                right += 1
                substr += addition
                if right == len(sequence):
                    break
                searchstring = sequence[:left] + sequence[right:]
                addition = sequence[right]
            elif re.search(substr + ('.' * maxLengthOfBulge) + sequence[right+maxLengthOfBulge], searchstring) and bulgesLeft != 0:
                right += maxLengthOfBulge
                addition = sequence[right]
                substr = substr + ('.' * maxLengthOfBulge) + addition
                searchstring = sequence[:left] + sequence[right+maxLengthOfBulge:]
                bulgesLeft -= 1
                right += maxLengthOfBulge
                addition = sequence[right]
            else:
                break
        substr_list.append(substr)
    #return max(enumerate(substr_list), key=lambda x: len(x[1]))
    return substr_list

teststring = "fsjaodpfknspaseniorprojectdfladkmdfisaneordrodectvjnsad"
position = 13

leftindex = return_leftmost_index(teststring, position)
longest = return_longest_rev_comp_bulge(teststring, position, leftindex, 1, 4)
print(longest)

"""
TODO List:
- Create tests for all functions
- Fix bulge structure identifier to work with bulges larger than size 1
- Allow for more than 1 bulge structure. COMPLETE
- Create left index identifying function for bulges.
- Modify this function, along with return_longest_rev_comp_bulge, to work with genetic sequences
- Create functions for internal loops and do the same modifications
"""