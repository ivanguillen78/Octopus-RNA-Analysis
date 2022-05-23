import csv
import os.path


def create_edit_dict(csv_file):
    """
    Returns a dictionary with ids as keys and list of edits as values.
    """
    edit_dict = {}
    with open(csv_file, newline="", encoding="utf8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["orf"] not in edit_dict:
                edit_dict[row["orf"]] = [int(row["pos"])]
            else:
                edit_dict[row["orf"]].append(int(row["pos"]))
    csvfile.close()
    return edit_dict


def create_output_csv(file_name, score_list, minLength=5):
    """
    Creates output csv file.
    """
    with open(file_name, "w", newline="", encoding="utf8") as csvfile:
        fieldnames = [
            "id",
            "position",
            "length",
            "base_string",
            "base_string_loc",
            "rev_comp",
            "rev_comp_loc",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for item in score_list:
            if item["length"] >= minLength:
                writer.writerow(
                    {
                        "id": item["id"],
                        "position": item["position"],
                        "length": item["length"],
                        "base_string": item["base"],
                        "base_string_loc": item["base_location"],
                        "rev_comp": item["rev_comp"],
                        "rev_comp_loc": item["rev_comp_location"],
                    }
                )
    csvfile.close()


def valid_file(extension, fileName, parser):
    ext = os.path.splitext(fileName)[1][1:]
    if ext != extension:
        parser.error(
            "The file "
            + fileName
            + " must have the {} file extension".format(extension)
        )
    return fileName
