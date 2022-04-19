#!/usr/bin/python
"""Main source file"""

from argparse import ArgumentParser, SUPPRESS
import os.path

from secondary_structure import (
    create_edit_dict,
    create_output_csv,
    find_secondary_structures,
)

parser = ArgumentParser(
    description="Identify secondary structures in genetic sequence",
    add_help=False
)
required = parser.add_argument_group('required arguemnts')
optional = parser.add_argument_group('optional arguments')

optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=SUPPRESS,
    help='show this help message and exit'
)

def valid_file(extension, fileName):
    ext = os.path.splitext(fileName)[1][1:]
    if ext != extension:
        parser.error("Your file must have the {} file extension".format(extension))
    return fileName


required.add_argument(
    "-f",
    "--fasta",
    type=lambda file: valid_file("fasta", file),
    metavar="",
    required=True,
    help="fasta file containing genetic sequences",
)
required.add_argument(
    "-c",
    "--csv",
    type=lambda file: valid_file("csv", file),
    metavar="",
    required=True,
    help="csv file containing edit positions",
)
required.add_argument(
    "-of",
    "--outputFile",
    type=lambda file: valid_file("csv", file),
    metavar="",
    required=True,
    help="name of output csv file"
)
optional.add_argument(
    "-ml",
    "--minLength",
    type=int,
    metavar="",
    help="minimum length of reverse complement (Default:5)",
)
args = parser.parse_args()


def main():
    edit_dict = create_edit_dict(args.csv)
    score_dict = find_secondary_structures(edit_dict, args.fasta)
    if args.minLength is not None:
        create_output_csv(args.outputFile, score_dict, args.minLength)
    else:
        create_output_csv(args.outputFile, score_dict)
    print("Output succesfully generated and stored in %s." % args.outputFile)


if __name__ == "__main__":
    main()
