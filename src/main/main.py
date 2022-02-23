#!/usr/bin/python
"""Main source file"""

import sys

from secondary_structure import create_edit_dict

edit_dict = create_edit_dict(sys.argv[1])

for k, v in edit_dict.items():
    print(k, v)
