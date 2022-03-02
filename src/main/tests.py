#!/usr/bin/python

from secondary_structure import (
    create_edit_dict,
    find_secondary_structures,
    secondary_structure,
)


def test_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {list}


def test_secondary_structure_regular():
    # test middle
    seq_1 = "AGCGTAGCTAGCTAGCTGACTGCTAGTAGCTAGCTACGCTAGTGCATGCAT"
    #        (((((^((((((((............))))))))))))))...........
    pos_1 = [5]
    assert secondary_structure(pos_1, seq_1) == [14]


def test_secondary_structure_left():
    # test left edge
    seq_2 = "GCTAGCCAGCTAGCGCTACGTAGCATCGATCGTACGATGCATCGATCGATC"
    #        ^(((((..)))))).....................................
    pos_2 = [0]
    assert secondary_structure(pos_2, seq_2) == [6]


def test_secondary_structure_right():
    # test right edge
    seq_3 = "AGCATGGCTAGCTAGCTGACTGCTAAAAGCTAGCTACGCTAGTGCCTGCAT"
    #        ...(((..........................................))^
    pos_3 = [len(seq_3) - 1]
    assert secondary_structure(pos_3, seq_3) == [3]


def test_secondary_structure_multiple():
    seq_4 = "AGCTAGTCAGGCGGGACTCAAATCATGCATGAACATGATTTGTTTTCACGA"
    #        .....(^(......))).(^(((((((......))))))))).........
    pos_4 = [6, 19]
    assert secondary_structure(pos_4, seq_4) == [3, 9]


def test_create_score_dict():
    edit_dict = create_edit_dict("../../data/testcsv.csv")
    score_dict = find_secondary_structures(edit_dict, "../../data/testfasta.fasta")
    test_score_dict = {
        "lcl|TRINITY_DN155104_c0_g1_i1:196-417": [5],
        "lcl|TRINITY_DN1563_c0_g1_i5:188-964": [5],
    }
    assert score_dict == test_score_dict
