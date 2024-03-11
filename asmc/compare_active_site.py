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
        
id_dict = utils.read_asmc_output({}, file1)
id_dict = utils.read_asmc_output(id_dict, file2, empty=False)
id_dict, ref_set = utils.read_identity_target_ref(id_dict, file_id)
id_dict = utils.compute_levenshtein(id_dict)
utils.write_output(id_dict, ref_set)