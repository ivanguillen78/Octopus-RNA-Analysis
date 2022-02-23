#!/usr/bin/python
import sys

sys.path.append("../main/")
from secondary_structure import *


def test_types():
    edit_dict = create_edit_dict("../../data/aes_profile.csv")
    assert set(map(type, edit_dict)) == {str}
    assert set(map(type, edit_dict.values())) == {list}
