#!/usr/bin/python
"""Main source file"""

import secondary_structure
import sys

edit_dict = secondary_structure.create_edit_dict(sys.argv[1])

for k,v in edit_dict.items():
    print(k,v)