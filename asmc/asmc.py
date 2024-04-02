import sys
import logging
import re
from itertools import groupby
from pathlib import Path
from typing import Dict, List, Tuple, Union, Any, Optional

import numpy as np
from sklearn import preprocessing
from sklearn.cluster import DBSCAN
from PIL import Image, ImageDraw

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    import weblogo
    
###############
## Functions ##
###############

## -------------------------- Pocket detection ------------------------------ ##

def build_ds(ref: Path, outdir: Path, chains: str) -> Tuple[Path, str]:
    """Build dataset file for p2rank

    Args:
        ref (patlib.Path): Path to reference file
        outdir (pathlib.Path): Path to the output directory
        chains (str): String indicating which chain to search

    Returns:
        ds (pathlib.Path): Path to the dataset file
        text (str): Text to write in the dataset file
    """
    
    # Use the 1st reference in the file
    try:
        pdb = ref.read_text().split("\n")[0]
        if not Path(pdb).exists():
            raise FileNotFoundError(f"{pdb} file not found")

    except Exception as error:
        raise Exception(f"An error has occured while reading {ref}:\n{error}")
    
    # Detect which chain
    ds = Path.joinpath(outdir, "data.ds")
    if chains == "all":
        chains = "*"
    
    # Writng the file
    text = f"HEADER: protein chains\n\n{pdb} {chains}"
            
    return ds, text

def extract_pocket(outdir: Path) -> Dict[str, List[int]]:
    """Extract the pocket posistions
    
    Reads the p2rank outputs to extract the positions of the best pocket that
    doesn't overlap several chains

    Args:
        outdir (pathlib.Path): Path to the output directory

    Returns:
        res_dict (dict): Dict containing as key the chain and as values the positions
    """
    
    try:
        prediction = [f for f in outdir.iterdir() if f.match("*predictions*")][0]
    except IndexError:
        raise RuntimeError(f"No predictions file after running p2rank")
    
    pred_arr = np.loadtxt(prediction, skiprows=1, delimiter=",",
                          converters={x:conv for x in range(11)},
                          dtype=str)
    
    # We go through the pockets sorted from best score to worst.
    # For each pocket, we check whether it's on a single chain.
    # If so, we retrieve it's positions else we continue.
    
    res_dict = {}
    
    try:
        for i in range(len(pred_arr)):
            res_str = pred_arr[i][9]
            res_list = [(elem.split("_")[0],
                            elem.split("_")[1]) for elem in res_str.split()]
            
            groups = groupby(res_list, key=lambda x: x[0])
            for key, g in groups:
                res_dict[key] = sorted([int(elem[1]) for elem in g])
                
            if len(res_dict) > 1:
                res_dict = {}
                continue
            else:
                break
    except Exception as error:
        raise Exception("An error has occured while reading prediction file "
                        f"from p2rank:\n{error}")
    
    return res_dict

def conv(x: str) -> str:
    
    return x.strip()

def build_pocket_text(ref: Path, res_dict: Dict[str, List[int]], outdir: Path,
                      query_chain: str) -> Tuple[Path, str]:
    """Build the pocket file

    Args:
        ref (patlib.Path): Path to reference file
        res_dict (dict): Dict containing as key the chain and as values the positions
        outdir (pathlib.Path): Path to the output directory
        query_chain (str): 'all' or value from the command line option

    Returns:
        output (pathlib.Path): Path of the pocket output file
        text (str): The text to write in the pocket output file
    """
    
    # Get the path of the pdb file used for p2rank
    pdb = Path(ref.read_text().split("\n")[0])
    # Get the file name without the extension
    pdb_id = pdb.stem
    output = Path.joinpath(outdir, "pocket.csv")
    
    try:
        chain = list(res_dict.keys())[0]
    except Exception as error:
        raise Exception("None results for p2rank, this may be due to an incorrect "
                        f"--chain value : {query_chain}")
    res_str = ''
    for elem in res_dict[chain]:
        res_str += f',{elem}'
    
    text = f"{pdb_id},{chain}{res_str}"
    
    return output, text

## ----------------------- Strcutural alignment ----------------------------- ##

def renumber_residues(ref_list: Tuple[Path, str, List[int]]) -> List[int]:
    """Renumbering reference structure

    Args:
        ref_list (tuple or list): [pathlib.Path, str(chain), list(positions)]

    Returns:
        renum (list): Renumbered positions
    """
    
    pdb = ref_list[0]
    chain = ref_list[1]
    true_pos = ref_list[2]
    renum = []
    i = 0
    resn = None
    with open(pdb, "r") as f:
        try:
            for num, line in enumerate(f):
                if line.startswith("ATOM") and line[21:22] == chain and \
                    line[17:20] != "HOH":
                    if resn is None:
                        resn = line[22:26].strip()
                    elif line[22:26].strip() != resn:
                        resn = line[22:26].strip()
                        i += 1
                    
                    if resn in true_pos and i not in renum:
                        renum.append(i)
        except IndexError:
            logging.error(f"An error has occured while reading {pdb}, line {num}"
                          f" seems to be incorrectly formatted\n:{line.strip()}")
            sys.exit(1)
    
    
    if len(renum) != len(true_pos):
        logging.error("An error has occured when renumbering the residues of "
                      "reference. This may be caused by a residue number indicated "
                      f"in the pocket file not found in the {pdb} or a "
                      "duplicated number")
        sys.exit(1)
    
    return renum

def extract_aligned_pos(id_ref: str, id_model: str,
                        ref_list: Tuple[Path, str, List[int]],
                        alignment_file: Path, keep_ref: bool) -> str:
    """Get positions aligned with the reference pocket

    Args:
        id_ref (str): Reference id
        id_model (str): Model id
        ref_list (tuple or list): [pathlib.Path, str(chain), list(positions), list(renum)]
        alignment_file (pathlib.Path): Path of the alignment file
        keep_ref (bool): Indicate whether we write the reference positions

    Returns:
        text (str): text to add to output file
    """
    
    renum_pos = ref_list[-1]
    pos_str = []
    aln = {id_ref:"", id_model:""}
    text = "" 
    
    # Read alignment from Usalign
    with open(alignment_file, "r") as f:
        ref = False
        for line in f:
            if line.startswith(">"):
                if f"{id_ref}.pdb" in line:
                    ref = True
                else:
                    ref = False
            else:
                if line.startswith("#"):
                    break
                if ref == True:
                    aln[id_ref] = line.strip()
                else:
                    aln[id_model] = line.strip()
    
    # If keep_ref is True, we write the active site of the ref
    if keep_ref == True:
        text += f">{id_ref}\n"
        # j is a counter for amino acids, incremented only when
        # an amino acid is encoutered
        j = 0
        pocket = ""
        for i, aa in enumerate(aln[id_ref]):
            if aa != "-":
                if j in renum_pos:
                    # We store i to be able to extract the aligned position
                    # in the target aligned sequence
                    pos_str.append(i)
                    pocket += aa
                j += 1
        text += pocket + "\n"
        
    else:
        j = 0
        for i, aa in enumerate(aln[id_ref]):
            if aa != "-":
                if j in renum_pos:
                    pos_str.append(i)
                j += 1
    
    # Get the positions aligned with the reference active site
    text += f">{id_model}\n"
    try:
        pocket = "".join([aln[id_model][i] for i in pos_str])
    except:
        pocket = ""
    text += pocket + "\n"
    
    return text

def build_multiple_alignment(pairwise_dir: Path, ref_file: Path,
                             pocket_file: Path) -> str:
    """Build multiple alignment

    Args:
        pairwise_dir (pathlib.Path): Path to the directory containing all
                                     pairwise alignments
        ref_file (pathlib.Path): Path to reference file
        pocket_file (pathlib.Path): Path to pocket file
    Returns:
        text (str): Text (multiple alignment) to write in the output
    """

    ref_pos = {}  # dict to store information like path, pocket chain and positions
    
    # Read the reference file and get paths
    with open(ref_file, "r") as f:
        for line in f:
            id_ref = Path(line.strip()).stem
            ref_pos[id_ref] = [line.strip()]
            
    # Reading the pocket file
    with open(pocket_file, "r") as f:
        for line in f:
            split_line = line.strip().split(",")
            ref_pos[split_line[0]].extend([split_line[1], split_line[2:]])
    
    # Renumbering references
    logging.info("Renumbering references")
    for key in ref_pos:
        renum = renumber_residues(ref_pos[key])
        ref_pos[key].append(renum)
    
    
    # Build multiple alignment
    all_pairwise = [f for f in pairwise_dir.iterdir()]
    count_ref = {ref:0 for ref in ref_pos}
    
    text = ""
    logging.info(f"Build multiple alignment")
    for i, alignment_file in enumerate(all_pairwise):
        id_model = alignment_file.stem.split("_-")[0]
        id_ref = alignment_file.stem.split("_-")[1]
        
        if count_ref[id_ref] == 0:
            text += extract_aligned_pos(id_ref, id_model, ref_pos[id_ref],
                                   alignment_file, keep_ref=True)
        else:
            text += extract_aligned_pos(id_ref, id_model, ref_pos[id_ref],
                                   alignment_file, keep_ref=False)
        
        count_ref[id_ref] += 1
    
    return text

## ----------------------- Multiple Sequence Alignment ---------------------- ##

def search_active_site_in_msa(msa: Path) -> str:
    """Search and extract active site in  MSA

    Args:
        msa (pathlib.Path): file with references, active site positions for each
        references, path to tsv file indicating the reference of each sequence
        and the path of msa

    Returns:
        str: multiple alignment of active sites to write in a file
    """
    
    ref = {}
    aln = ""
    id_file = ""
    
    # Read -M/--msa input file
    with open(msa, "r") as f:
        for i, line in enumerate(f):
            split_line = line.strip().split(",")
            if len(split_line) > 1:
                ref[split_line[0]] = {"pos":[], "aln":[], "seq":""}
                ref[split_line[0]]["pos"] = [int(x) -1 for x in split_line[1:]]
            else:
                path = Path(line.strip())
                if path.suffix in [".fasta", ".fa", ".faa"]:
                    aln = Path(line.strip())
                else:
                    id_file = Path(line.strip())
            
    if not aln.exists():
        logging.error("An error has occured while reading "
                      f"'{msa}':\n'{aln}' doesn't exists")
        sys.exit(1)
    
    if id_file != "":
        if not id_file.exists():
            logging.error("An erro has occured while reading "
                        f"'{msa}':\n'{id_file}' doesn't exists")
            sys.exit(1)
        elif id_file.exists:
            # Get the reference for each sequences
            map_target_ref = {}
            with open(id_file, "r") as f:
                for line in f:
                    split_line = line.strip().split()
                    map_target_ref[split_line[0]] = split_line[1]
    
    # Parse the multiple sequences alignment
    all_seq = {}
    with open(aln, "r") as f:
        seq_id = ""
        for line in f:
            if line.startswith(">"):
                if "tr|" in line or "sp|" in line:
                    print(line)
                    seq_id = re.search("\\w+\\|(\\w+)\\|\\w+", line).group(1)
                else:
                    seq_id = line[1:].split()[0]
                    for c in [".", ":", "|"]:
                            seq_id = seq_id.replace(c, "_")
            else:
                if seq_id in ref:
                    ref[seq_id]["seq"] += line.strip()
                else:
                    try:
                        all_seq[seq_id] += line.strip()
                    except KeyError:
                        all_seq[seq_id] = line.strip()
    
    text = ""
    
    # For each reference sequence, we extract its active site via the
    # positions indicated in the file given to the -M/--msa argument. Then we
    # extract all the characters aligned with these positions in the sequences
    # with this reference
    for ref_id in ref:
        text += f">{ref_id}\n"
        j = 0
        for i, aa in enumerate(ref[ref_id]["seq"]):
            if aa != "-":
                if j in ref[ref_id]["pos"]:
                    ref[ref_id]["aln"].append(i)
                    text += aa
                j += 1
        text += '\n'

        if len(ref) > 1:
            ref_seq = [s for s in map_target_ref if map_target_ref[s] == ref_id]
            for seq in ref_seq:
                text += f">{seq}\n"
                text += "".join(all_seq[seq][i] for i in ref[ref_id]["aln"])
                text += "\n"
        else:
            for seq in all_seq:
                text += f">{seq}\n"
                text += "".join(all_seq[seq][i] for i in ref[ref_id]["aln"])
                text += "\n"
    
    return text

## ------------------------------- Clustering ------------------------------- ##

def read_alignment(file: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Read sequences alignment in fasta format

    Args:
        file (pathlib.Path): Path of the fasta file

    Returns:
        sequences (dict): Dictionnary with the sequence id as key and
                          the corresponding sequence as value
        removed (dict) : Dictionnary containing the removed sequences
    """
    
    sequences = {}
    removed = {}
    with open(file, "r") as f:
        seq_id = ""
        for line in f:
            if line.startswith(">"):
                if "tr|" in line or "sp|" in line:
                    seq_id = re.search("\\w+\\|(\\w+)\\|\\w+", line).group(1)
                    sequences[seq_id] = ""
                else:
                    seq_id = line[1:].split()[0]
                    sequences[seq_id] = ""
                    
            else:
                n = len(re.findall("-", line))
                try:
                    if (n / len(line.strip())) > 0.1:
                        del sequences[seq_id]
                        removed[seq_id] = line.strip()
                    else:
                        sequences[seq_id] += line.strip()
                except ZeroDivisionError:
                    del removed[seq_id]
                    removed[seq_id] = "empty line"

    return sequences, removed

def read_matrix(matrix: Path) -> Dict[str, Dict[str, int]]:
    """Read the distance matrix in tsv format

    Args:
        matrix (pathlib.Path): Path of the tsv file

    Returns:
        scoring_dict (dict): Dictionnary of dictionnary for access to each distance between amino acid
    """
    
    if not matrix.exists():
        logging.error(f"{matrix} dosen't exist")
        sys.exit(1)
    
    scoring_dict  = {}
    with open(matrix, "r") as f:
        aa_order = []
        for i, line in enumerate(f):
            if i == 0:
                aa_order = line.strip().split("\t")
                if len(aa_order) == 1:
                    logging.error(f"{matrix} seems to not be a tsv file")
                    sys.exit(1)
            else:
                sl = line.strip().split("\t")
                try:
                    scoring_dict[sl[0]] = {aa_order[i]:int(sl[1:][i])
                                        for i in range(len(aa_order))}
                except:
                    logging.error("An error has occured while reading the "
                                  "distances matrix. The matrix may not be "
                                  "symetrical")
                    sys.exit(1)
                    
    return scoring_dict

def pairwise_score(scoring_dict: Dict[str, Dict[str, int]], seqA: str, seqB: str,
                   weighted_pos: List[int]) -> int:
    """Compute the score (distance) between two sequences

    Args:
        scoring_dict (dict): Dictionnary of dictionnary for access to each distance between amino acid
        seqA (str): A sequence
        seqB (str): A sequence
        weigted_pos (list): List containing positions with more weight for score calculation

    Returns:
        score (int): The score
    """
    
    score = 0
    for i, (posA, posB) in enumerate(zip(seqA, seqB)):
        if posA in ["-", "X"] or posB in ["-", "X"]:
            if i+1 in weighted_pos:
                score += 20 * 5
            else:
                score += 20
        else:
            try:
                if i+1 in weighted_pos:
                    score += scoring_dict[posA][posB] * 5
                else:
                    score += scoring_dict[posA][posB]
            except KeyError:
                logging.warning("At leats one of these two characters isn't "
                                f"in the distances matrix: {posA} {posB}, they "
                                "are given the same score as '-' and 'X'")
                
                if i+1 in weighted_pos:
                    score += 20 * 5
                else:
                    score += 20
    
    return score

def dissimilarity(sequences: Dict[str, str], scoring_dict: Dict[str, Dict[str, int]],
                  weighted_pos: List[int]) -> Tuple[List[str], np.ndarray]:
    """Build the dissimilarity/distance matrix

    Args:
        sequences (dict): Dictionnary with the sequence id as key and the corresponding sequence as value
        scoring_dict (dict): Dictionnary of dictionnary for access to each distance between amino acid
        weigted_pos (list): List containing positions with more weight for score calculation

    Returns:
        data (np.ndarray): The distance matrix
        key_list (list):  The list of sequences id
    """
    
    data = []
    key_list = list(sequences.keys())
        
    for i, key1 in enumerate(key_list):
        row = []
        for j, key2 in enumerate(key_list):
            
            if key1 == key2:
                score = 0.0
            else:
                score = pairwise_score(scoring_dict,
                                       sequences[key1],
                                       sequences[key2],
                                       weighted_pos)
                
            row.append(score)
        data.append(row)
    data = np.asarray(data)
    
    data = preprocessing.MinMaxScaler().fit_transform(
        X=data.reshape(-1,1)
    ).reshape(data.shape)

    return key_list, data

def dbscan_clustering(data: np.ndarray, threshold: float, min_samples: int,
                      threads: int) -> np.ndarray:
    """Clustering with DBSCAN

    Args:
        data (np.ndarray): The distance matrix
        threshold (float): The maximum disatnce between two samples for one to
                           be considered as in the neighborhood of the other
        min_samples (int): The number of samples in an neighborhood for a point
                           to be considered as a core point
        threads (int): Number of CPU threads

    Returns:
        labels (np.ndarray): Cluster id for each sequences
    """
    
    try:
        dbscan = DBSCAN(eps=threshold, metric="precomputed", n_jobs=threads,
                        min_samples=min_samples)
        
        labels = dbscan.fit_predict(X=data)
    except Exception as error:
        logging.error(f"An error has occured during the clustering:\n{error}")
        sys.exit(1)
    
    return labels

def formatting_output(sequences: Dict[str, str], key_list: List[str],
                      labels: np.ndarray) -> List[Tuple[str, str, int]]:
    """Format data to write output

    Args:
        sequences (dict): Dictionnary with the sequence id as key and the corresponding sequence as value
        key_list (list):  The list of sequences id present in the clustering
        labels (np.ndarray): Cluster id for each sequences

    Returns:
        G (list): Sorted data
    """
    
    try:
        G = [(key_list[i], sequences[key_list[i]], n) for i, n in enumerate(labels)]
        G = sorted(G, key=lambda x: x[2])
    except Exception as error:
        logging.error(f"An error has occured durint output formatting:\n{error}")
        sys.exit(1)
    
    return G

## ------------------------------- Weblogo ---------------------------------- ##

def build_fasta(group: List[Tuple[str, str, Optional[Any]]]) -> str:
    """Build group fasta

    Args:
        group (list): List the sequences for each sequence id in the group

    Returns:
        fasta (pathlib.Path): The path of the output file
        text (str): The text to write in the output file
    """
    
    text = ""
    for elem in group:
        text += f">{elem[0]}\n{elem[1]}\n"
    
    return text

def build_logo(lenght: int, fasta: Path, outdir: Path, n: int, prefix: str,
               out_format: str) -> None:
    """Build weblogo for a Group

    Args:
        lenght (int): Number of sequences in the group
        fasta (pathlib.Path): The path of the fasta file
        outdir (pathib.Path): Output directory
        n (int): Cluster id
        prefix (str): Prefix for the weblogo title
        out_format (str): eps or png

    """
    
    with open(fasta, "r") as fin:
        seqs = weblogo.read_seq_data(fin)
    
    try:
        data = weblogo.LogoData.from_seqs(seqs)
        options = weblogo.LogoOptions()
        options.logo_title = f"{prefix}{n}"
        options.fineprint = str(lenght)
        options.color_scheme = weblogo.chemistry
        logo_format = weblogo.LogoFormat(data, options)
        if out_format == "png":
            logo_bytes = weblogo.png_print_formatter(data, logo_format)
            output = Path.joinpath(outdir, f"{prefix}{n}.png")
        elif out_format == "eps":
            logo_bytes = weblogo.eps_formatter(data, logo_format)
            output = Path.joinpath(outdir, f"{prefix}{n}.eps")
        
        output.write_bytes(logo_bytes)
        
    except Exception as error:
        logging.error(f"An error has occured when creating the logo of" 
                      f" {prefix}{n}:\n{error}")
    
def merge_logo(outdir: Path, n: int, prefix: str, out_format: str):
    """Merge single logo files

    Args:
        outdir (pathib.Path): Output directory
        n (int): The number of groups
        prefix (str): Prefix for the weblogo title 
        out_format (str): eps or png
    """

    LOGO_PAD = 10

    IM_WIDTH = 500
    
    all_file = [f for f in outdir.iterdir() if f.match(f"{prefix}*.{out_format}")]
    
    try:
        name_to_find = outdir / f"{prefix}-1.{out_format}"
        i = all_file.index(name_to_find)
        del all_file[i]
        all_file.append(name_to_find)
    except:
        pass
    
    
    if out_format == "png":
        LOGO_WIDTH = 450
        LOGO_HEIGHT = 175
        im_height = n * LOGO_HEIGHT + ((n-1) * LOGO_PAD) + LOGO_PAD * 2
    else:
        LOGO_HEIGHT = 92*2
        im_height = n * LOGO_HEIGHT + ((n-1) * LOGO_PAD) + LOGO_PAD * 2
    
    img = Image.new(mode='RGB', size=(IM_WIDTH, im_height), color=(255,255,255))
    draw = ImageDraw.Draw(img)
    
    top_left_coord = (LOGO_PAD, LOGO_PAD)
    
    for i, f in enumerate(all_file):
        
        logo = Image.open(f)
        if out_format == "png":
            logo = logo.resize(size=(LOGO_WIDTH, LOGO_HEIGHT))
        else:
            logo.load(scale=2)
        
        if i != 0:
            top_left_coord = (top_left_coord[0],
                              top_left_coord[1] + LOGO_PAD + LOGO_HEIGHT)
            
        img.paste(logo, top_left_coord)
    
    output = Path.joinpath(outdir,
                           f"Groups_logo.{out_format}")
    
    img.save(output)
    
    for f in all_file:
        f.unlink()