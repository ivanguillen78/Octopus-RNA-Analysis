#!/usr/bin/python
import sys

from secondary_structure import create_edit_dict, secondary_structure


def test_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {list}


def test_secondary_structure():
    # test middle
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    pos_1 = [5]

    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    pos_2 = [0]

    # test right edge
    seq_3 = "AGCATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    pos_3 = [len(seq_3)-1]

    seq_4 = "AGCTAGTCAGGCGGGACTCAAATCATGCATGAACATGATTTGTTTTCACGA"
    #        .....(^(......))).(^(((((((......))))))))).........
    pos_4 = [6, 19]

    assert secondary_structure(pos_1, seq_1) == [14]
    assert secondary_structure(pos_2, seq_2) == [6]
    assert secondary_structure(pos_3, seq_3) == [3]
    assert secondary_structure(pos_4, seq_4) == [3, 9]
