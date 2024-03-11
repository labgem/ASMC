import argparse
import sys
import utils
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, metavar="", required=True,
                    help="tsv file from run_asmc.py")
parser.add_argument("-p", "--position", type=int, metavar="", required=True,
                    help="position where to find the specified amino acid type, e.g: 5")
parser.add_argument("-a", "--aa-type", type=str, metavar="",
                    required=True ,help="amino acid type to search, must be "
                    "either 1-letter amino acid, 'aromatic', 'acidic', 'basic'"
                    ", 'polar' or 'hydrophobic")
parser.add_argument("-g", "--group", type=int, metavar="",
                    help="group id, if not used search in all groups")

args = parser.parse_args()

tsv_file = Path(args.file)
# Check if the tsv file exist
if not tsv_file.exists():
    print(f"error: argument -f/--file '{tsv_file}' doesn't exists")
    sys.exit(1)
    
result = utils.extract_aa(tsv_file, pos=args.position,
                          aa=args.aa_type, group=args.group)

print(result, end="")