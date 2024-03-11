import argparse
import sys
import utils
from pathlib import Path

parser = argparse.ArgumentParser()
input_opt = parser.add_mutually_exclusive_group(required=True)
input_opt.add_argument("-s", "--seqs", type=str, metavar="",
                       help="multi fasta file")
input_opt.add_argument("-m", "--models", type=str, metavar="",
                       help="file containing all PDB paths")
parser.add_argument("-r", "--ref", type=str, metavar="",
                    help="file contaning the reference structure paths")
args = parser.parse_args()

if args.seqs is not None:
    if not Path(args.seqs).exists():
        print(f"{__file__}: error: argument -s/--seqs '{args.seqs}' doesn't exist")
        sys.exit(1)
    print("Read sequences")
    targets_seq = utils.read_multi_fasta(args.seqs)
elif args.models is not None:
    if not Path(args.models).exists():
        print(f"{__file__}: error: argument -m/--models '{args.models}' doesn't exist")
        sys.exit(1)
    print("Reading models")
    targets_seq = utils.read_models(args.models)
    
print("Reading references")
ref_seq = utils.read_models(args.ref)
    
text = ""
    
print("Global alignment")
print(f"0/{len(targets_seq)}", end="\r")
for i, t in enumerate(targets_seq):
    if t in ref_seq:
        continue
    ref_max, perc_id_max = utils.get_identity(ref_seq, targets_seq[t])
    text += f"{t}\t{ref_max}\t{perc_id_max:.2f}\n"
    print(f"{i+1}/{len(targets_seq)}", end="\r")
    
output = Path.joinpath(Path.cwd().absolute(), "identity_target_ref.tsv")
output.write_text(text)
print(end="\x1b[2K")
print("done")