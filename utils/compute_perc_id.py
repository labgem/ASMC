import argparse
import sys
import re
import warnings
from pathlib import Path
from Bio import Align
from Bio.Align import substitution_matrices


###############
## Functions ##
###############

def get_seq_from_pdb(pdb):
    """Get sequence from pdb file

    Args:
        pdb (pathlib.Path): Pdb file

    Returns:
        str: The sequence
    """

    mapping = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
               'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',
               'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
               'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
    
    seq = ""
    
    with open(pdb, "r") as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM') and \
                line[17:20] != 'HOH':
                    try:
                        seq += mapping[line[17:20]]
                    except KeyError:
                        seq += 'X'
    
    return seq


def read_models(models):
    """Reads the model file

    For each model, add its id as key in a dictionnary and add as value the
    sequence return by get_seq_from_pdb().

    Args:
        models (pathlib.Path): The file containing the model paths

    Returns:
        dict: A dictionnary with each pair of id - seq
    """
    
    all_seq = {}
    
    with open(models, "r") as f:
        for line in f: 
            model = Path(line.strip().split()[0])
            model_id = model.stem
            all_seq[model_id] = get_seq_from_pdb(model)
    
    return all_seq

def read_multi_fasta(fasta):
    """Reads a multi fasta file
    
    Add in a dictionnary all pair id - sequence.

    Args:
        fasta (pathlib.Path): The multi fasta file

    Returns:
        dict: A dictionnary with each pair of id - seq
    """
    
    all_seq = {}
    
    with open(fasta, "r") as f:
        seq_id = ""
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                if "tr|" in seq_id or "sp|" in seq_id:
                    seq_id = re.search("\\|(\\w+)\\|", seq_id).group(1)
                    all_seq[seq_id] = ""
                else:
                    for c in [".", ":", "|"]:
                        seq_id = seq_id.replace(c, "_")
                    all_seq[seq_id] = ""
            else:
                all_seq[seq_id] += line.strip()

    return all_seq

def get_identity(ref_seq, target):
    """Get the % identity between two sequences
    
    For each pair reference - target, build a global alignment and calculates
    the percentage of identity.

    Args:
        ref_seq (dict): Dictionnary with ids of reference as key and their
                        sequences as value
        target (str): The target sequence

    Returns:
        (str, float): The reference id with the best identity
                      The percentage identity
    """
    
    # Set the alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    ref_max = ""
    perc_id_max = 0
    
    # For each ref, align with the target sequence
    for ref in ref_seq:
        
        aln = aligner.align(ref_seq[ref], target)
        score = 0
        for c in str(aln[0]).split()[1]:
            if c == "|":
                score += 1
            
        perc_id = score/len(str(aln[0]).split()[1]) *100
        if perc_id > perc_id_max:
            perc_id_max = perc_id
            ref_max = ref
    
    return ref_max, perc_id_max

##########
## Main ##
##########

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    input_opt = parser.add_mutually_exclusive_group(required=True)
    input_opt.add_argument("-s", "--seqs", type=str, metavar="",
                           help="multi fasta file")
    input_opt.add_argument("-m", "--models", type=str, metavar="",
                           help="file containing all PDB paths")
    parser.add_argument("-r", "--ref", type=str, metavar="",
                        help="file contaning the reference structure paths")
    
    args = parser.parse_args()
    
    warnings.filterwarnings('ignore', module='Bio')
    
    if args.seqs is not None:
        if not Path(args.seqs).exists():
            print(f"{__file__}: error: argument -s/--seqs '{args.seqs}' doesn't exist")
            sys.exit(1)
        print("Read sequences")
        targets_seq = read_multi_fasta(args.seqs)
    elif args.models is not None:
        if not Path(args.models).exists():
            print(f"{__file__}: error: argument -m/--models '{args.models}' doesn't exist")
            sys.exit(1)
        print("Reading models")
        targets_seq = read_models(args.models)
    
    print("Reading references")
    ref_seq = read_models(args.ref)
    
    text = ""
    
    print("Global alignment")
    print(f"0/{len(targets_seq)}", end="\r")
    for i, t in enumerate(targets_seq):
        if t in ref_seq:
            continue
        ref_max, perc_id_max = get_identity(ref_seq, targets_seq[t])
        text += f"{t}\t{ref_max}\t{perc_id_max:.2f}\n"
        print(f"{i+1}/{len(targets_seq)}", end="\r")
    
    output = Path.joinpath(Path.cwd().absolute(), "identity_target_ref.tsv")
    output.write_text(text)
    print(end="\x1b[2K")
    print("done")