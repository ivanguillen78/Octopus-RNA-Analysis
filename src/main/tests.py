#!/usr/bin/python

from secondary_structure import checkLeft, checkRight, create_edit_dict


def test_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {int}


def test_secondary_structure_regular():
    # test middle
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    pos_1 = 5
    assert len(checkRight(seq_1, pos_1) + checkLeft(seq_1, pos_1)) - 1 == 14


def test_secondary_structure_leftedge():
    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    pos_2 = 0
    assert len(checkRight(seq_2, pos_2) + checkLeft(seq_2, pos_2)) - 1 == 6


def test_secondary_structure_rightedge():
    # test right edge
    seq_3 = "AGTATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    pos_3 = len(seq_3) - 1
    assert len(checkRight(seq_3, pos_3) + checkLeft(seq_3, pos_3)) - 1 == 3


def test_secondary_structure_multiple():
    seq_4 = "AGCTAGTCAGGCGGGACTCAAATCATGCATGAACATGATTTGTTTTCACGA"
    #        .....(^(......))).(^(((((((......))))))))).........
    pos_4a = 6
    pos_4b = 19
    assert len(checkRight(seq_4, pos_4a) + checkLeft(seq_4, pos_4a)) - 1 == 3
    assert len(checkRight(seq_4, pos_4b) + checkLeft(seq_4, pos_4b)) - 1 == 9
