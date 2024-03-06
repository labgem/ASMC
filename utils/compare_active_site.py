from pathlib import Path
import argparse
import sys


###############
## Functions ##
###############

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

def write_output(id_dict, ref_set):
    """Write the output file

    Args:
        id_dict (dict): Dict with sub dict as value and seqID as key
        ref_set (set): Set containing the reference IDs

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
        elif d1 is None:
            text += "f2\n"
        elif d2 is None:
            text += "f1\n"
        elif ref is None:
            text += "None\n"
        
    output = Path(Path.cwd().absolute(), "active_site_checking.tsv")
    output.write_text(text)

##########
## Main ##
##########

if __name__ == "__main__":
    
    
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
        
    id_dict = read_asmc_output({}, file1)
    id_dict = read_asmc_output(id_dict, file2, empty=False)
    id_dict, ref_set = read_identity_target_ref(id_dict, file_id)
    id_dict = compute_levenshtein(id_dict)
    write_output(id_dict, ref_set)
