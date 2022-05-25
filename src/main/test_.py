#!/usr/bin/python

# from files import create_edit_dict

from files import create_edit_dict
from hairpin import (
    return_leftmost_index_hairpin,
    return_longest_hairpin,
    create_secondary_structure_hairpin,
)
from internal_loop import (
    return_leftmost_index_loop,
    return_longest_rev_comp_loop,
    checkInternalLoop,
    create_secondary_structure_loop,
)


def test_edit_dict_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {list}


"""
HAIRPIN TESTS
"""


def test_leftmost_index_hairpin_regular():
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    #        ^ pos 0
    pos_1 = 5
    assert return_leftmost_index_hairpin(seq_1, pos_1) == 0


def test_longest_rev_comp_hairpin_regular():
    # test middle
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    pos_1 = 5
    leftindex = return_leftmost_index_hairpin(seq_1, pos_1)
    assert len(return_longest_hairpin(seq_1, pos_1, leftindex)[0]) == 14


def test_create_sec_struct_hairpin_regular():
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    pos_1 = 5
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_hairpin(
        seq_1, pos_1
    )
    assert length == 14
    assert base == "AGCGTAGCTAGCTA"
    assert base_loc == [0, 13]
    assert rev_comp == "TAGCTAGCTACGCT"
    assert rev_comp_loc == [26, 39]


def test_leftmost_index_hairpin_leftedge():
    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    #        ^ pos 0
    pos_2 = 0
    assert return_leftmost_index_hairpin(seq_2, pos_2) == 0


def test_longest_rev_comp_hairpin_leftedge():
    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    pos_2 = 0
    leftindex = return_leftmost_index_hairpin(seq_2, pos_2)
    assert len(return_longest_hairpin(seq_2, pos_2, leftindex)[0]) == 6


def test_create_sec_struct_hairpin_leftedge():
    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    pos_2 = 0
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_hairpin(
        seq_2, pos_2
    )
    assert length == 6
    assert base == "GCTAGC"
    assert base_loc == [0, 5]
    assert rev_comp == "GCTAGC"
    assert rev_comp_loc == [8, 13]


def test_leftmost_index_hairpin_rightedge():
    # test right edge
    seq_3 = "AGTATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    #                                                        ^ pos 48
    pos_3 = 50
    assert return_leftmost_index_hairpin(seq_3, pos_3) == 48


def test_longest_rev_comp_hairpin_rightedge():
    # test right edge
    seq_3 = "AGTATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    pos_3 = 50
    leftindex = return_leftmost_index_hairpin(seq_3, pos_3)
    assert len(return_longest_hairpin(seq_3, pos_3, leftindex)[0]) == 3


def test_create_sec_struct_hairpin_rightedge():
    # test right edge
    seq_3 = "AGTATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    pos_3 = len(seq_3) - 1
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_hairpin(
        seq_3, pos_3
    )
    assert length == 3
    assert base == "CAT"
    assert base_loc == [48, 50]
    assert rev_comp == "ATG"
    assert rev_comp_loc == [3, 5]


"""
INTERNAL LOOP TESTS
"""

# Test 1: Regular internal loop with loop length 1, numloops 1


def test_longest_int_loop1():
    seq4 = "ACGCTCGATGTCAGCCGTACGTACGATCGATCTGGCATCGTATC"
    #       .....(^(((.(((.................))).)))))....
    pos4 = 12
    leftindex = return_leftmost_index_loop(seq4, pos4)
    assert leftindex == 5
    assert len(return_longest_rev_comp_loop(seq4, pos4, leftindex)[1]) == 9


def test_create_int_loop1():
    seq4 = "ACGCTCGATGTCAGCCGTACGTACGATCGATCTGGCATCGTATC"
    #       .....(^(((.(((.................))).)))))....
    pos4 = 12
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq4, pos4
    )

    assert length == 9
    assert base == "CGATG.CAG"
    assert base_loc == [5, 13]
    assert rev_comp == "CTG.CATCG"
    assert rev_comp_loc == [31, 39]


# Test 2: Left internal loop with loop length 1, numloops 1


def test_longest_int_loop1():
    seq4 = "CGATGTCAGCCGTACGTACGATCGATCTGGCATCGTATC"
    #       (^(((.(((.................))).)))))....
    pos4 = 0
    leftindex = return_leftmost_index_loop(seq4, pos4)
    assert leftindex == 0
    assert len(return_longest_rev_comp_loop(seq4, pos4, leftindex)[1]) == 9


def test_create_int_loop1():
    seq4 = "CGATGTCAGCCGTACGTACGATCGATCTGGCATCGTATC"
    #       ^((((.(((.................))).)))))....
    pos4 = 0
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq4, pos4
    )

    assert length == 9
    assert base == "CGATG.CAG"
    assert base_loc == [0, 8]
    assert rev_comp == "CTG.CATCG"
    assert rev_comp_loc == [26, 34]


# Test 3: Right internal loop with loop length 1, numloops 1


def test_longest_int_loop1():
    seq4 = "TGCATCTGGCATCGCCGTACGTACGATCGATCGATGTCAG"
    #       .....(((.(((((.................))))).))^
    pos4 = len(seq4) - 1
    leftindex = return_leftmost_index_loop(seq4, pos4)
    assert leftindex == 31
    assert len(return_longest_rev_comp_loop(seq4, pos4, leftindex)[1]) == 9


def test_create_int_loop1():
    seq4 = "TGCATCTGGCATCGCCGTACGTACGATCGATCGATGTCAG"
    #       .....(((.(((((.................))))).))^
    pos4 = len(seq4) - 1
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq4, pos4
    )

    assert length == 9
    assert base == "CGATG.CAG"
    assert base_loc == [31, 39]
    assert rev_comp == "CTG.CATCG"
    assert rev_comp_loc == [5, 13]
