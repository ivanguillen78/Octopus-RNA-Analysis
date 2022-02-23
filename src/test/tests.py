#!/usr/bin/python
from secondary_structure import create_edit_dict
import sys

sys.path.append("../main/")


def test_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {list}
