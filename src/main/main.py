#!/usr/bin/python
'''Main source file'''

from secondary_structure import *
import sys

edit_dict = create_edit_dict(sys.argv[1])

for k,v in edit_dict.items():
    print(k,v)