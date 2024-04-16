from pathlib import Path
import argparse
import sys
import utils


parser = argparse.ArgumentParser()
parser.add_argument("-f1", type=str, required=True, metavar='',
                    help="Group file 1")
parser.add_argument("-f2", type=str, required=True, metavar="",
                    help="Group file 2")
parser.add_argument("-id", type=str, metavar="", required=True,
                    help="identity_target_ref.tsv")

args = parser.parse_args()

file1 = Path(args.f1)
file2 = Path(args.f2)
file_id = Path(args.id)

if not file1.exists():
    print(f"{__file__}: error: argument -f1 '{file1}' doesn't exist")
    sys.exit(1)
if not file2.exists():
    print(f"{__file__}: error: argument -f1 '{file2}' doesn't exist")
    sys.exit(1)
if not file_id.exists():
    print(f"{__file__}: error: argument -id '{file_id}' doesn't exist")
    sys.exit(1)
        
id_dict = utils.build_comparison_data({}, file1)
id_dict = utils.build_comparison_data(id_dict, file2, empty=False)
id_dict, ref_set = utils.add_ref_data_to_comparison_data(id_dict, file_id)
id_dict = utils.compute_levenshtein(id_dict)
text = utils.build_active_site_checking_file(id_dict, ref_set)
output = Path.cwd().absolute() / "active_site_checking.tsv"

# Hides 3 last columns
text = "\n".join(["\t".join(elem.split("\t")[:-3]) for elem in text.split("\n")])
output.write_text(text)