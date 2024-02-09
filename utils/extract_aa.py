import argparse
import sys
from pathlib import Path

# Dictionnary containing all possibilities
aa_type = {"F":["F"], "W":["W"], "Y":["Y"], "A":["A"], "V":["V"], "L":["L"],
           "I":["I"], "P":["P"], "M":["M"], "G":["G"], "S":["S"], "T":["T"],
           "C":["C"], "Q":["Q"], "N":["N"], "E":["E"], "D":["D"], "H":["H"],
           "R":["R"], 'K': ["K"], "aromatic":['F','W','Y'], "acidic":['E','D'],
           "basic":['H','K','R'], "polar":['S','T','C','Y','Q','N'],
           "hydrophobic":['A','V','L','I','P','W','F','M','G']}

# Command line Parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, metavar="", required=True,
                    help="tsv file from run_asmc.py")
parser.add_argument("-p", "--position", type=int, metavar="", required=True,
                    help="position where to find the specified amino acid type, e.g: 5")
parser.add_argument("-a", "--aa-type", type=str, metavar="", choices=list(aa_type.keys()),
                    required=True ,help="amino acid type to search, must be "+
                    "either 1-letter amino acid, 'aromatic', 'acidic', 'basic'"+
                    ", 'polar' or 'hydrophobic")
parser.add_argument("-g", "--group", type=int, metavar="",
                    help="group id, if not used search in all groups")

args = parser.parse_args()

pos = args.position
aa_list = aa_type[args.aa_type]


tsv_file = Path(args.file)
# Check if the tsv file exist
if not tsv_file.exists():
    print(f"{__file__}: error: argument -f/--file '{tsv_file}' doesn't exists")
    sys.exit(1)
    
lines_to_show = ""

# No specified group
if args.group is None:
    with tsv_file.open(mode="r") as f:
        for line in f:
            try:
                sequence = line.split("\t")[1]
            except IndexError:
                print(f"{__file__}: error: argument -f/--file '{tsv_file}' "+
                      "seems to not be a tsv file or don't contains at least 2 columns")
                sys.exit(1)
            
            try:
                if sequence[pos-1] in aa_list:
                    lines_to_show += line
            except IndexError:
                print(f"{__file__}: error: argument -p/--position '{pos}' "+
                      f"must be bewteen 1 and {len(sequence)}")
                sys.exit(1)

# Specified group
else:
    with tsv_file.open(mode="r") as f:
        for line in f:
            try:
                sequence = line.split()[1]
                group = line.split()[2]
            except IndexError:
                print(f"{__file__}: error: argument -f/--file '{tsv_file}' "+
                      "seems to not be a tsv file or don't contains at 3 columns")
                sys.exit(1)
            
            
            if str(args.group) == group:
                try:
                    if sequence[pos-1] in aa_list:
                        lines_to_show += line
                except IndexError:
                    print(f"{__file__}: error: argument -p/--position '{pos}' "+
                      f"must be bewteen 1 and {len(sequence)}")
                    sys.exit(1)
            
print(lines_to_show, end="")
