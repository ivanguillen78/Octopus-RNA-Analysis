#!/usr/bin/python

# from files import create_edit_dict

from argparse import ArgumentParser

from files import create_edit_dict, valid_file
from hairpin import (
    create_secondary_structure_hairpin,
    return_leftmost_index_hairpin,
    return_longest_hairpin,
)
from internal_loop import (
    create_secondary_structure_loop,
    return_leftmost_index_loop,
    return_longest_rev_comp_loop,
)

"""
File tests
"""


def test_edit_dict_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {list}


def test_valid_file():
    testparser = ArgumentParser(
        description="Parser created for testing valid file extensions", add_help=False
    )
    testfile1 = "testfile.csv"
    testfile2 = "testfile.fasta"
    assert valid_file("csv", testfile1, testparser) == testfile1
    assert valid_file("fasta", testfile2, testparser) == testfile2


"""
HAIRPIN TESTS
"""

# Test 1: Regular hairpin loop test


def test_longest_hairpin1():
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    pos_1 = 5
    leftindex = return_leftmost_index_hairpin(seq_1, pos_1)
    assert leftindex == 0
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


# Test 2: Left hairpin loop test


def test_leftmost_index_hairpin_leftedge():
    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    #        ^ pos 0
    pos_2 = 0
    leftindex = return_leftmost_index_hairpin(seq_2, pos_2)
    assert leftindex == 0
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


# Test 3: Right hairpin loop test


def test_leftmost_index_hairpin_rightedge():
    # test right edge
    seq_3 = "AGTATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    #                                                        ^ pos 48
    pos_3 = 50
    leftindex = return_leftmost_index_hairpin(seq_3, pos_3)
    assert leftindex == 48
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


def test_longest_int_loop2():
    seq5 = "CGATGTCAGCCGTACGTACGATCGATCTGGCATCGTATC"
    #       (^(((.(((.................))).)))))....
    pos5 = 0
    leftindex = return_leftmost_index_loop(seq5, pos5)
    assert leftindex == 0
    assert len(return_longest_rev_comp_loop(seq5, pos5, leftindex)[1]) == 9


def test_create_int_loop2():
    seq5 = "CGATGTCAGCCGTACGTACGATCGATCTGGCATCGTATC"
    #       ^((((.(((.................))).)))))....
    pos5 = 0
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq5, pos5
    )

    assert length == 9
    assert base == "CGATG.CAG"
    assert base_loc == [0, 8]
    assert rev_comp == "CTG.CATCG"
    assert rev_comp_loc == [26, 34]


# Test 3: Right internal loop with loop length 1, numloops 1


def test_longest_int_loop3():
    seq6 = "TGCATCTGGCATCGCCGTACGTACGATCGATCGATGTCAG"
    #       .....(((.(((((.................))))).))^
    pos6 = len(seq6) - 1
    leftindex = return_leftmost_index_loop(seq6, pos6)
    assert leftindex == 31
    assert len(return_longest_rev_comp_loop(seq6, pos6, leftindex)[1]) == 9


def test_create_int_loop3():
    seq6 = "TGCATCTGGCATCGCCGTACGTACGATCGATCGATGTCAG"
    #       .....(((.(((((.................))))).))^
    pos6 = len(seq6) - 1
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq6, pos6
    )

    assert length == 9
    assert base == "CGATG.CAG"
    assert base_loc == [31, 39]
    assert rev_comp == "CTG.CATCG"
    assert rev_comp_loc == [5, 13]


# Test 4: Regular internal loop with loop length 1, numloops 2 (multiloop)


def test_longest_int_loop4():
    seq7 = "ACGCTCGATGTCAGCCGTACGTACGATCGATCTGGCACCGAGCTATC"
    #       ..(^(((.((.(((.................))).)).)))))....
    pos7 = 3
    leftindex = return_leftmost_index_loop(seq7, pos7, numLoops=2)
    assert leftindex == 2
    assert len(return_longest_rev_comp_loop(seq7, pos7, leftindex, numLoops=2)[1]) == 12


def test_create_int_loop4():
    seq7 = "ACGCTCGATGTCAGCCGTACGTACGATCGATCTGGCACCGAGCTATC"
    #       ..(^(((.((.(((.................))).)).)))))....
    pos7 = 3
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq7, pos7, numLoops=2
    )

    assert length == 12
    assert base == "GCTCG.TG.CAG"
    assert base_loc == [2, 13]
    assert rev_comp == "CTG.CA.CGAGC"
    assert rev_comp_loc == [31, 42]


# Test 5: Left internal loop with loop length 1, numloops 2 (multiloop)


def test_longest_int_loop5():
    seq8 = "GCTCGATGTCAGCCGTACGTACGATCGATCTGGCACCGAGCTATC"
    #       ^((((.((.(((.................))).)).)))))....
    pos8 = 0
    leftindex = return_leftmost_index_loop(seq8, pos8, numLoops=2)
    assert leftindex == 0
    assert len(return_longest_rev_comp_loop(seq8, pos8, leftindex, numLoops=2)[1]) == 12


def test_create_int_loop5():
    seq8 = "GCTCGATGTCAGCCGTACGTACGATCGATCTGGCACCGAGCTATC"
    #       ^((((.((.(((.................))).)).)))))....
    pos8 = 0
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq8, pos8, numLoops=2
    )

    assert length == 12
    assert base == "GCTCG.TG.CAG"
    assert base_loc == [0, 11]
    assert rev_comp == "CTG.CA.CGAGC"
    assert rev_comp_loc == [29, 40]


# Test 6: Right internal loop with loop length 1, numloops 2 (multiloop)


def test_longest_int_loop6():
    seq8 = "ACGCTCGATGTCAGCCGTACGTACGATCGATCTGGCACCGAGC"
    #       ..(((((.((.(((.................))).)).))))^
    pos8 = len(seq8) - 1
    leftindex = return_leftmost_index_loop(seq8, pos8, numLoops=2)
    assert leftindex == 31
    assert len(return_longest_rev_comp_loop(seq8, pos8, leftindex, numLoops=2)[1]) == 12


def test_create_int_loop6():
    seq8 = "ACGCTCGATGTCAGCCGTACGTACGATCGATCTGGCACCGAGC"
    #       ..(((((.((.(((.................))).)).))))^
    pos8 = len(seq8) - 1
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq8, pos8, numLoops=2
    )

    assert length == 12
    assert base == "CTG.CA.CGAGC"
    assert base_loc == [31, 42]
    assert rev_comp == "GCTCG.TG.CAG"
    assert rev_comp_loc == [2, 13]


# Test 7: Regular internal loop with max loop length 3, numloops 3


def test_longest_int_loop7():
    seq9 = "ACGGTCCGATGTCAAACGTTTATGCGATGTCAGATGCATGGACGCCCGACACCGCTGGCA"
    #       ......((.((^(...(((..((((..........))))..)))...)))).))......
    pos9 = 11
    leftindex = return_leftmost_index_loop(seq9, pos9, maxLengthOfLoop=3, numLoops=3)
    assert leftindex == 6
    assert (
        len(
            return_longest_rev_comp_loop(
                seq9, pos9, leftindex, maxLengthOfLoop=3, numLoops=3
            )[1]
        )
        == 19
    )


def test_create_int_loop7():
    seq9 = "ACGGTCCGATGTCAAACGTTTATGCGATGTCAGATGCATGGACGCCCGACACCGCTGGCA"
    #       ......((.((^(...(((..((((..........))))..)))...)))).))......
    pos9 = 11
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq9, pos9, maxLengthOfLoop=3, numLoops=3
    )

    assert length == 19
    assert base == "CG.TGTC...CGT..ATGC"
    assert base_loc == [6, 24]
    assert rev_comp == "GCAT..ACG...GACA.CG"
    assert rev_comp_loc == [35, 53]


# Test 8: Left internal loop with max loop length 3, numloops 3


def test_longest_int_loop8():
    seq10 = "CGATGTCAAACGTTTATGCGATGTCAGATGCATGGACGCCCGACACCGCTGGCA"
    #        ^(.((((...(((..((((..........))))..)))...)))).))......
    pos10 = 0
    leftindex = return_leftmost_index_loop(seq10, pos10, maxLengthOfLoop=3, numLoops=3)
    assert leftindex == 0
    assert (
        len(
            return_longest_rev_comp_loop(
                seq10, pos10, leftindex, maxLengthOfLoop=3, numLoops=3
            )[1]
        )
        == 19
    )


def test_create_int_loop8():
    seq10 = "CGATGTCAAACGTTTATGCGATGTCAGATGCATGGACGCCCGACACCGCTGGCA"
    #        ^(.((((...(((..((((..........))))..)))...)))).))......
    pos10 = 0
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq10, pos10, maxLengthOfLoop=3, numLoops=3
    )

    assert length == 19
    assert base == "CG.TGTC...CGT..ATGC"
    assert base_loc == [0, 18]
    assert rev_comp == "GCAT..ACG...GACA.CG"
    assert rev_comp_loc == [29, 47]


# Test 9: Right internal loop with loop length 3, numloops 3


def test_longest_int_loop9():
    seq11 = "ACGGTCCGGTGTCAAACGTTTATGCGATGTCAGATGCATGGACGCCCGACACCG"
    #        ......(((((((...(((..((((..........))))..)))...))))))^
    pos11 = len(seq11) - 1
    leftindex = return_leftmost_index_loop(seq11, pos11, maxLengthOfLoop=3, numLoops=3)
    assert leftindex == 35
    # assert len(return_longest_rev_comp_loop(seq11, pos11, leftindex, maxLengthOfLoop=3, numLoops=3)[1]) == 19


def test_create_int_loop9():
    seq11 = "ACGGTCCGGTGTCAAACGTTTATGCGATGTCAGATGCATGGACGCCCGACACCG"
    #        ......(((((((...(((..((((..........))))..)))...))))))^
    pos11 = len(seq11) - 1
    length, base, base_loc, rev_comp, rev_comp_loc = create_secondary_structure_loop(
        seq11, pos11, maxLengthOfLoop=3, numLoops=3
    )

    assert length == 19
    assert base == "GCAT..ACG...GACACCG"
    assert base_loc == [35, 53]
    assert rev_comp == "CGGTGTC...CGT..ATGC"
    assert rev_comp_loc == [6, 24]
