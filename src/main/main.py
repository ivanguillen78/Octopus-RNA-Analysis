#!/usr/bin/python
"""Main source file"""

from argparse import ArgumentParser, SUPPRESS
from alive_progress import alive_bar
import fastaparser


from files import create_edit_dict, create_output_csv, valid_file
from hairpin import create_secondary_structure_hairpin
from internal_loop import create_secondary_structure_loop


parser = ArgumentParser(
    description="Identify secondary structures in genetic sequence", add_help=False
)
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")

optional.add_argument(
    "-h",
    "--help",
    action="help",
    default=SUPPRESS,
    help="show this help message and exit",
)

required.add_argument(
    "-f",
    "--fasta",
    type=lambda file: valid_file("fasta", file, parser),
    metavar="",
    required=True,
    help="fasta file containing genetic sequences",
)
required.add_argument(
    "-c",
    "--csv",
    type=lambda file: valid_file("csv", file, parser),
    metavar="",
    required=True,
    help="csv file containing edit positions",
)
required.add_argument(
    "-of",
    "--outputFile",
    type=lambda file: valid_file("csv", file, parser),
    metavar="",
    required=True,
    help="name of output csv file",
)
optional.add_argument(
    "-ml",
    "--minLength",
    type=int,
    metavar="",
    help="minimum length of reverse complement (Default:5)",
)
optional.add_argument(
    "-ll",
    "--loopLength",
    type=int,
    metavar="",
    help="max num of mismatches in loop (Default:1)",
)
optional.add_argument(
    "-nl",
    "--numLoops",
    type=int,
    metavar="",
    help="max number of loops allowed in loop structure (Default:1)",
)
args = parser.parse_args()


def find_secondary_structures(edit_dict, fasta_file):
    """
    Iterates through edit dictionary and finds longest reverse complement
    for each edit.
    """
    output_list = []
    with open(fasta_file, encoding="utf8") as fastafile:
        edit_dict_size = 0
        for edit in edit_dict:
            edit_dict_size += len(edit_dict[edit])
        numSequences = int(len(fastafile.readlines()) / 2)
        parser = fastaparser.Reader(fastafile)
        sequence_dict = {}
        print("Preparing data")
        with alive_bar(numSequences, bar="smooth", spinner="fish2") as bar:
            for seq in parser:
                bar()
                sequence_dict[seq.id] = seq.sequence_as_string()
        print("Data preparation complete. Locating secondary structures.")
        with alive_bar(edit_dict_size, bar="smooth", spinner="fish2") as bar:
            for edit in edit_dict:
                if edit in sequence_dict:
                    for pos in edit_dict[edit]:
                        bar()
                        (
                            length,
                            base,
                            base_loc,
                            rev,
                            rev_loc,
                        ) = create_secondary_structure_loop(sequence_dict[edit], pos)
                        output_list.append(
                            {
                                "id": edit,
                                "position": pos,
                                "length": length,
                                "base": base,
                                "base_location": base_loc,
                                "rev_comp": rev,
                                "rev_comp_location": rev_loc,
                            }
                        )
    fastafile.close()
    return output_list


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
