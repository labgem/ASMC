import re
import warnings
from pathlib import Path
from typing import Dict, Tuple, Set, Union, Optional
from Bio import Align
from Bio.Align import substitution_matrices

warnings.filterwarnings('ignore', module='Bio')


class FileFormatError(Exception):
    """Exception raised for file which not contain 2 columns

    Attribute:
        file (pathlib.Path): file which caused the error
        n (int): number of columns
    """
    def __init__(self, file: Path, n: int) -> None:
        self.file = file
        self.message = f"'{file}' does not appear to contain at least {n} columns"
        super().__init__(self.message)

class PositionError(Exception):
    """Exception raised for a position not between 1 and the sequence size

    Attributes:
        pos (int): the position which caused the error
        limit (int): the sequence size
    """
    
    def __init__(self, pos: int, limit: int) -> None:
        self.pos = pos
        self.limit = limit
        self.message = f"position must be between 1 and {limit}, got '{pos}'"
        super().__init__(self.message)

class AminoAcidTypeError(Exception):
    """Exception raised for Amino Acid which does not correspond to a 1-letter
    code or a valid amino acid type
    
    valid amino acid type : 'aromatic', 'acidic', 'basic', 'polar', 'hydrophobic'

    Attribute:
        aa (str): the amino acid string which caused the error
    """
    
    def __init__(self, aa: str) -> None:
        self.aa = aa
        self.message = "expected 1-letter amino acid or an amino acid type"
        self.message += f", got '{aa}'"
        super().__init__(self.message)

def get_seq_from_pdb(pdb: Path) -> str:
    """Get sequence from pdb file

    Args:
        pdb (pathlib.Path): Pdb file

    Returns:
        sequence (str): The sequence
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

def read_models(models: Path) -> Dict[str, str]:
    """Reads the model file

    For each model, add its id as key in a dictionary and add as value the
    sequence return by get_seq_from_pdb().

    Args:
        models (pathlib.Path): The file containing the model paths

    Returns:
        all_seq (dict): A dictionary with each pair of id - seq
    """
    
    all_seq = {}
    
    with open(models, "r") as f:
        for line in f: 
            model = Path(line.strip().split()[0])
            model_id = model.stem
            all_seq[model_id] = get_seq_from_pdb(model)
    
    return all_seq

def read_multi_fasta(fasta: Path) -> Dict[str, str]:
    """Reads a multi fasta file
    
    Add in a dictionary all pair id - sequence.

    Args:
        fasta (pathlib.Path): The multi fasta file

    Returns:
        all_seq (dict): A dictionary with each pair of id - seq
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

def get_identity(ref_seq: Dict[str, str], target: str) -> Tuple[str, float]:
    """Get the % identity between two sequences
    
    For each pair reference - target, build a global alignment and calculate
    the percentage of identity.

    Args:
        ref_seq (dict): Dictionary with ids of reference as key and their
                        sequences as value
        target (str): The target sequence

    Returns:
        ref_max (str): The reference id with the best identity
        perc_id_max (float): The percentage identity
    """
    
    # Set the alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    ref_max = ""
    perc_id_max = 0.0
    
    # For each ref, align with the target sequence
    for ref in ref_seq:
        
        aln = aligner.align(ref_seq[ref], target)
        score = 0.0
        
        c = aln[0].counts()
        score = c.identities

        perc_id = score / aln[0].shape[1] * 100

        if perc_id > perc_id_max:
            perc_id_max = perc_id
            ref_max = ref
    
    return ref_max, perc_id_max

def build_comparison_data(id_dict: Optional[Dict[str, Union[str, int, None]]],
                          file: Path, 
                          empty=True) -> Dict[str, Union[str, int, None]]:
    """Read the ASMC groups.tsv and load information in a dictionary

    Args:
        id_dict (dict): An empty dictionary or containing sub dict as value of
                        seqID (key)
        file (pathlib.Path): The ASMC groups.tsv 
        empty (bool, optional): Default to True.

    Returns:
        id_dict (dict): The updated id_dict
    """
    
    with open(file, "r") as f:
        if empty is True:
            for line in f:
                split_line = line.strip().split()
                id_dict[split_line[0]] = {"f1":split_line[1], "f2":None,
                                          "d":None,
                                          "ref_id":None,
                                          "ref_seq":None,
                                          "ref_d1":None,
                                          "ref_d2":None,
                                          "ref_pid":None,
                                          'g1':split_line[-1],
                                          'g2':None}
                
        else:
            for line in f:
                split_line = line.strip().split()
                try:
                    id_dict[split_line[0]]["f2"] = split_line[1]
                    id_dict[split_line[0]]["g2"] = split_line[-1]
                except KeyError:
                    id_dict[split_line[0]] = {"f1":None, "f2":split_line[1],
                                              "d":None,
                                              "ref_id":None,
                                              "ref_seq":None,
                                              "ref_d1":None,
                                              "ref_d2":None,
                                              "ref_pid":None,
                                              'g1':None,
                                              'g2':split_line[-1]}
    
    return id_dict


def add_ref_data_to_comparison_data(id_dict:Dict[str, Union[str, int, None]],
                                    file:Path) -> Tuple[Dict[str, Union[str, int, None]], Set[str]]:
    """Reads the identity_targets_refs.tsv file and add information to a 
    comparison data

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key
        file (pathlib.Path): The identity_targets_refs.tsv

    Returns:
        id_dict (dict): The updated id_dict
        ref_set (set): Set containing the reference IDs
    """
    
    ref_set = set()
    
    with open(file, "r") as f:
        for line in f:
            split_line = line.strip().split()
            try:
                ref_set.add(split_line[1])
                id_dict[split_line[0]]["ref_id"] = split_line[1]
                id_dict[split_line[0]]["ref_pid"] = split_line[2]
                id_dict[split_line[0]]["ref_seq"] = id_dict[split_line[1]]["f1"]
            except KeyError:
                continue
            
    
    return id_dict, ref_set

def LD_two_rows(s1: str, s2: str) -> int:
    """Calculates Levenshtein distance between two strings
    
    Simple implementation of Levenshtein distance based on the two rows
    algorithm.

    Args:
        s1 (str): The first string
        s2 (str): The second string

    Returns:
        int: The Levenshtein/edit distance
    """

    # Switch s1 and s2 to reduce the columns number    
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

def compute_levenshtein(id_dict: Dict[str,Union[str, int, None]]) -> Dict[str, Union[str, int, None]]:
    """Calculates all possible distances

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key

    Returns:
        id_dict (dict): The updated id_dict
    """
    
    for key in id_dict:
        seq1 = id_dict[key]["f1"]
        seq2 = id_dict[key]["f2"]
        seq_ref = id_dict[key]["ref_seq"]
        
        # SeqID present in f1 and f2 -> Levenshtein distance between their seqs
        if seq1 is not None and seq2 is not None:
            distance = LD_two_rows(seq1, seq2)
            id_dict[key]["d"] = distance
            
            # Levenshtein distance between seq1 - seq_ref and seq2 - seq_ref
            if seq_ref is not None:
                distance = LD_two_rows(seq_ref, seq1)
                id_dict[key]["ref_d1"] = distance
                distance = LD_two_rows(seq_ref, seq2)
                id_dict[key]["ref_d2"] = distance
        
        # Levenshtein distance between seq_ref - seq1 or seq2
        elif seq_ref is not None:
            if seq1 is not None:
                distance = LD_two_rows(seq_ref, seq1)
                id_dict[key]["ref_d1"] = distance
            else:
                distance = LD_two_rows(seq_ref, seq2)
                id_dict[key]["ref_d2"] = distance
    
    return id_dict

def build_active_site_checking_file(id_dict: Dict[str, Union[str, int, None]],
                                    ref_set: Set[str]) -> str:
    """format text of active site checking

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key
        ref_set (set): Set containing the reference IDs

    Returns:
        text (str): The text to write

    """
    
    # Headers
    text = "ID\tG1\tSEQ1\tG2\tSEQ2\tDIFF\tREF_ID\tPERC_ID\tREF_SEQ\tD_REF_1\tD_REF_2\tNEAR_REF\n"
    for key in id_dict:
        if key in ref_set:
            continue
        seq1 = id_dict[key]["f1"]
        seq2 = id_dict[key]["f2"]
        d = id_dict[key]["d"]
        ref = id_dict[key]["ref_id"]
        seq_ref = id_dict[key]["ref_seq"]
        d1 = id_dict[key]["ref_d1"]
        d2 = id_dict[key]["ref_d2"]
        g1 = id_dict[key]["g1"]
        g2 = id_dict[key]["g2"]
        ref_pid = id_dict[key]["ref_pid"]
        
        # Add dictionary items
        text += f"{key}\t{g1}\t{seq1}\t{g2}\t{seq2}\t{d}\t{ref}\t{ref_pid}\t"
        text += f"{seq_ref}\t{d1}\t{d2}\t"
        
        # Last column
        if d1 is not None and d2 is not None:
            if d1 < d2:
                text += "f1"
            elif d1 > d2:
                text += "f2"
            else:
                text += "both"
        elif ref is None:
            text += "None"
        elif d1 is None and d2 is None:
            if seq1 is None:
                text += "f2"
            elif seq2 is None:
                text += "f1"
            else:
                text += "None"
        elif d1 is None:
            text += "f2"
        elif d2 is None:
            text += "f1"
            
        text += "\n"
        
    return text
    
def extract_aa(file: Path, pos: int, aa: str, group: Optional[int]):
    """Extracts sequences with a specific amino acid type at a given position

    Args:
        file (pathlib.Path): File path
        pos (int): Position in the sequence 
        aa (str): 1-letter amino acid or amino acid type
        group (int, optional): Group id

    Returns:
        result (str): The extracted lines
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
        raise AminoAcidTypeError(aa)
    
    result = ""
    
    with open(file, "r") as f:
        # No specified group
        if group is None:
            for line in f:
                try:
                    sequence = line.split("\t")[1]
                except IndexError:
                    raise FileFormatError(file, 2)
                
                try:
                    if sequence[pos-1] in aa_list:
                        result += line
                except IndexError:
                    raise PositionError(pos, len(sequence))                 
                    
        # Specified group
        else:
            for line in f:
                try:
                    sequence = line.split()[1]
                    group_line = line.split()[2]
                except IndexError:
                    raise FileFormatError(file, 2)
                    
                if str(group) == group_line:
                    try:
                        if sequence[pos-1] in aa_list:
                            result += line
                    except IndexError:
                        raise PositionError(pos, len(sequence))
                        
    return result

def get_unique(group_file: Path) -> Tuple[Dict[str, Tuple[str, Set[str]]],
                                         Dict[str, Tuple[int, int, float]]]:
    """Calculates statistics on the number of unique sequences per group

    Args:
        group_file (Path): ASMC tsv output

    Raises:
        FileFormatError: Raised if the tsv contains less than 3 columns

    Returns:
        unique_seq (dict): Dict with seq as key and a tuple containing
        the group id and a set of sequence ids as values
        
        groups_stats (dict): Dict with group id as key and a tuple of int and
        float as value
    """
    unique_seq = {}
    group_stats = {}
    
    with open(group_file, "r") as f:
        for line in f:
            splitted = line.strip().split()
            if len(splitted) < 3:
                raise FileFormatError(group_file, 3)
            
            if splitted[2] not in unique_seq:
                group_stats[splitted[2]] = (0, 0, 0.0)
            
            if splitted[1] not in unique_seq:
                unique_seq[splitted[1]] = (splitted[2], {splitted[0]})
            else:
                unique_seq[splitted[1]][1].add(splitted[0])
    
    for group in group_stats:
        seq_list = [(s, len(unique_seq[s][1]))
                    for s in unique_seq if unique_seq[s][0] == group]
    
        total = 0
        nb_unique =  len(seq_list)
        diff_list = []
        mean_diff = 0.0
        
        for i in range(nb_unique-1):
            total += seq_list[i][1]
            
            for j in range(i+1, nb_unique):
                d = LD_two_rows(seq_list[i][0], seq_list[j][0])
                diff_list.append(d)
                
        total += seq_list[nb_unique-1][1]
        try:
            mean_diff = sum(diff_list) / len(diff_list)
        except ZeroDivisionError:
            mean_diff = 0.0
        
        group_stats[group] = (nb_unique, total, mean_diff)
        
    return unique_seq, group_stats