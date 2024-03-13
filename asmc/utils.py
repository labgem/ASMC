import sys
import re
import warnings
from pathlib import Path
from Bio import Align
from Bio.Align import substitution_matrices

warnings.filterwarnings('ignore', module='Bio')

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
    
    set_res = set()
    
    with open(pdb, "r") as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM') and \
                line[17:20] != 'HOH':
     
                    res = line[22:26]
                    if res not in set_res:
                        set_res.add(res)
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

def read_asmc_output(id_dict, file, empty=True):
    """Read the ASMC groups.tsv

    Args:
        id_dict (dict): An empty dictionnary or contaning sub dict as value of
                        seqID (key)
        file (pathlib.Path): The ASMC groups.tsv 
        empty (bool, optional): Defaults to True.
        
        For the first file (f1) empty is True so a sub dictionnary is crate and
        associate to the seqID as key. The sequence is add at the 'f1' key in
        the sub dict.
        For the 2nd file (f2) empty is set to False in the main, so if the
        seqID is already in id_dict, the value of 'f2' key in the sub dict is
        update. Otherwise the sub dict is create for the seqID 

    Returns:
        dict: The updated id_dict
    """
    
    with open(file, "r") as f:
        if empty == True:
            for line in f:
                split_line = line.strip().split()
                id_dict[split_line[0]] = {"f1":split_line[1], "f2":None, "d":None,
                                          "ref":{"id":None,
                                                 "seq":None,
                                                 "d1":None,
                                                 "d2":None,
                                                 "pid":None}}
                
        else:
            for line in f:
                split_line = line.strip().split()
                try:
                    id_dict[split_line[0]]["f2"] = split_line[1]
                except KeyError:
                    id_dict[split_line[0]] = {"f1":None, "f2":split_line[1], "d":None,
                                              "ref":{"id":None,
                                                     "seq":None,
                                                     "d1":None,
                                                     "d2":None,
                                                     "pid":None}}
    
    return id_dict


def read_identity_target_ref(id_dict, file):
    """Reads the identity_target_ref.tsv file

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key
        file (pathlib.Path): The identity_target_ref.tsv

    Returns:
        (dict, set): The updated id_dict,
                     Set containing the reference IDs
    """
    
    ref_set = set()
    
    with open(file, "r") as f:
        for line in f:
            split_line = line.strip().split()
            try:
                ref_set.add(split_line[1])
                id_dict[split_line[0]]["ref"]["id"] = split_line[1]
                id_dict[split_line[0]]["ref"]["pid"] = split_line[2]
                id_dict[split_line[0]]["ref"]["seq"] = id_dict[split_line[1]]["f1"]
            except KeyError:
                continue
            
    
    return id_dict, ref_set

def LD_two_rows(s1, s2):
    """Calcultes Levenshtein distance between two strings
    
    Simple implementation of Levenshtein distance based on the two rows
    algorithm.

    Args:
        s1 (str): The first string
        s2 (str): The second string

    Returns:
        int: The Levenshtein/edit distance
    """

    # Switch s1 and s2 for reduce the columns number    
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    
    rowA = [j for j in range(len(s1)+1)]
    for i ,c2 in enumerate(s2):
        rowB = [i + 1]
        
        for j, c1 in enumerate(s1):
            if c1 == c2:
                cost = 0
            else:
                cost = 1
                
            rowB.append(min(rowB[j],
                            rowA[j],
                            rowA[j+1]) + cost)
        
        rowA = rowB
    
    return rowB[-1]

def compute_levenshtein(id_dict):
    """Calculates all possible distances

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key

    Returns:
        dict: The updated id_dict
    """
    
    for key in id_dict:
        seq1 = id_dict[key]["f1"]
        seq2 = id_dict[key]["f2"]
        seq_ref = id_dict[key]["ref"]["seq"]
        
        # SeqID present in f1 and f2 -> Levenshtein distance between their seqs
        if seq1 is not None and seq2 is not None:
            distance = LD_two_rows(seq1, seq2)
            id_dict[key]["d"] = distance
            
            # Levenshtein distance between seq1 - seq_ref and seq2 - seq_ref
            if seq_ref is not None:
                distance = LD_two_rows(seq_ref, seq1)
                id_dict[key]["ref"]["d1"] = distance
                distance = LD_two_rows(seq_ref, seq2)
                id_dict[key]["ref"]["d2"] = distance
        
        # Levenshtein distance between seq_ref - seq1 or seq2
        elif seq_ref is not None:
            if seq1 is not None:
                distance = LD_two_rows(seq_ref, seq1)
                id_dict[key]["ref"]["d1"] = distance
            else:
                distance = LD_two_rows(seq_ref, seq2)
                id_dict[key]["ref"]["d2"] = distance
    
    return id_dict

def build_active_site_checking_file(id_dict, ref_set):
    """format text of active site checking

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key
        ref_set (set): Set containing the reference IDs

    Returns:
        str: The text to write

    """
    
    # Headers
    text = "ID\tF1\tF2\tD\tREF_ID\tREF_SEQ\tD_REF_1\tD_REF_2\tNEAR_REF\n"
    for key in id_dict:
        if key in ref_set:
            continue
        seq1 = id_dict[key]["f1"]
        seq2 = id_dict[key]["f2"]
        d = id_dict[key]["d"]
        ref = id_dict[key]["ref"]["id"]
        seq_ref = id_dict[key]["ref"]["seq"]
        d1 = id_dict[key]["ref"]["d1"]
        d2 = id_dict[key]["ref"]["d2"]
        
        # Add dictionnary items
        text += f"{key}\t{seq1}\t{seq2}\t{d}\t{ref}\t{seq_ref}\t{d1}\t{d2}\t"
        
        # Last column
        if d1 is not None and d2 is not None:
            if d1 < d2:
                text += "f1\n"
            elif d1 > d2:
                text += "f2\n"
            else:
                text += "both\n"
        elif ref is None:
            text += "None\n"
        elif d1 is None and d2 is None:
            if seq1 is None:
                text += "f2\n"
            elif seq2 is None:
                text += "f1\n"
        elif d1 is None:
            text += "f2\n"
        elif d2 is None:
            text += "f1\n"
        
    return text
    
def extract_aa(file, pos, aa, group):
    """Extracts sequences with a specific amino acid type at a given position

    Args:
        file (pathlib.Path): File path
        pos (int): Position in the sequence 
        aa (str): 1-letter amino acid or amino acid type
        group (int): Group id

    Returns:
        str: The extracted lines
    """
    
    aa_type = {"F":["F"], "W":["W"], "Y":["Y"], "A":["A"], "V":["V"],
               "L":["L"], "I":["I"], "P":["P"], "M":["M"], "G":["G"],
               "S":["S"], "T":["T"], "C":["C"], "Q":["Q"], "N":["N"],
               "E":["E"], "D":["D"], "H":["H"], "R":["R"], 'K': ["K"],
               "aromatic":['F','W','Y'], "acidic":['E','D'],
               "basic":['H','K','R'], "polar":['S','T','C','Y','Q','N'],
               "hydrophobic":['A','V','L','I','P','W','F','M','G']}
    
    try:
        aa_list = aa_type[aa]
    except KeyError:
        print(f"error: argument -a/--aa-type isn't a 1-letter amino acid or a "
              "valid amino acid")
        sys.exit(1)
    
    result = ""
    
    with open(file, "r") as f:
        # No specified group
        if group is None:
            for line in f:
                try:
                    sequence = line.split("\t")[1]
                except IndexError:
                    print(f"error: argument -f/--file '{file}' seems to not be "
                          "a tsv file or does not contain at least 2 columns")
                    sys.exit(1)
                
                try:
                    if sequence[pos-1] in aa_list:
                        result += line
                except IndexError:
                    print(f"error: argument -p/--position '{pos}' must be "
                          f"between 1 and {len(sequence)}")
                    sys.exit(1)
        
        # Specified group
        else:
            for line in f:
                try:
                    sequence = line.split()[1]
                    group_line = line.split()[2]
                except IndexError:
                    print(f"error: argument -f/--file '{file}' seems to not be "
                          "a tsv file or does not contain at least 2 columns")
                    sys.exit(1)
                    
                if str(group) == group_line:
                    try:
                        if sequence[pos-1] in aa_list:
                            result += line
                    except IndexError:
                        print(f"error: argument -p/--position '{pos}' must be "
                              f"between 1 and {len(sequence)}")
                        sys.exit(1)
                        
    return result