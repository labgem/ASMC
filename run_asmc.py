import sys
import argparse
import datetime
import logging
import subprocess
import re
from itertools import groupby
from pathlib import Path

import yaml
import numpy as np
from sklearn import preprocessing
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    import weblogo

###############
## Functions ##
###############

def read_yaml():
    """Read the yaml file
    
    Load the content of the yaml file and check the validity of paths.
    The configuration should be place in a directory named ressources in the same
    location as run_asmc.py, e.g :

    .
    ├── ressources
    │   ├── AA_distances.tsv
    │   └── config.yml
    ├── run_asmc.py
    └── src
        ├── build_ali.py
        └── modeling.py


    Args:
        args (argparse.Namespace): the object containing all arguments

    Returns:
        yml (dict): a dictionary corresponding to the contents of the yaml file
    """
    
    main_path = Path(__file__).absolute()
    parent_path = main_path.parent
    ressources_path = Path.joinpath(parent_path, "ressources")
    
    config_path = Path.joinpath(ressources_path, "config_asmc.yml")

    if not Path.is_file(config_path):
        logging.error(f"not found the configuration file: {config_path}")
        sys.exit(1)
    
    with open(config_path, "r") as f:
        yml = yaml.safe_load(f)
        
    for key in yml:
        
        if key == "distances":
            if not Path(yml[key]).exists():
                logging.error(f"{yml[key]} doesn't exist")
        
        elif key == "usalign":
            command = f"{yml[key]} -h"
            try:
                ret = subprocess.run(command.split(), capture_output=True)
            except Exception as error:
                logging.error(f"An error has occured with USalign:\n{error} ")
                sys.exit(1)
            
            if ret.returncode != 0:
                logging.error(f"{ret.stderr.decode('utf-8')}")
                sys.exit(1)
        
        elif key == "weblogo":
            command = f"{yml[key]} --version"
            try:
                ret = subprocess.run(command.split(), capture_output=True)
            except Exception as error:
                logging.error(f"An error has occured with weblogo:\n{error}")
                sys.exit(1)
                
            if ret.returncode != 0:
                logging.error(f"{ret.stderr.decode('utf-8')}")
                sys.exit(1)
                
        elif key == "java":
            command = f"{yml[key]} --version"
            try:
                ret = subprocess.run(command.split(), capture_output=True)
            except Exception as error:
                logging.error(f"An error has occured with java:\n{error}")
                
            if ret.returncode != 0:
                logging.error(f"{ret.stderr.decode('utf-8')}")
                sys.exit(1)
        
    return yml

## -------------------------- Pocket detection ------------------------------ ##

def build_ds(ref, outdir, chains):
    """Build dataset file for p2rank

    Args:
        ref (patlib.Path): Path to reference file
        outdir (pathlib.Path): Path to the output directory
        chains (str): String indicating which chain to search

    Returns:
        ds (pathlib.Path): Path to the builded dataset file
    """
    
    # Use the 1st reference in the file
    logging.info("Using the 1st reference structure to detect pocket")
    try:
        pdb = ref.read_text().split("\n")[0]
        if not Path(pdb).exists():
            logging.error(f"Path to the 1st reference structure in reference file doesn't exist: {pdb}")
    except Exception as error:
        logging.error(f"An error has occured while reading {ref}:\n{error}")
        sys.exit(1)
    
    # Detect which chain
    ds = Path.joinpath(outdir, "data.ds")
    if chains == "all":
        chains = "*"
    
    
    # Writng the file
    text = f"HEADER: protein chains\n\n{pdb} {chains}"
    ds.write_text(text)
            
    return ds

def run_prank(yml, ds, outdir):
    """Run p2rank
    
    Execute a p2rank process with the subprocess module

    Args:
        yml (dict): a dictionary corresponding to the contents of the yaml file
        ds (pathlib.Path): Path to the dataset file
        outdir (pathlib.Path): Path to the output directory

    Returns:
        result (subprocess.CompletedProcess): The completed process
    """
    
    P2RANK = yml["p2rank"]
    
    # Execute p2rank
    command = f"{P2RANK} predict {ds} -o {outdir}"
    try:
        result = subprocess.run(command.split(), capture_output=True)
    except Exception as error:
        logging.error(f"An error has occured during the p2rank process:\n{error}")
        sys.exit(1)
    
    if result.returncode != 0:
            logging.error(f"An error has occured during the:\n{result.stderr.decode('utf-8')}")
            sys.exit(1)
    
    return result

def extract_pocket(outdir):
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
        logging.error(f"No predictions file after running p2rank")
        sys.exit(1)
    
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
        logging.error(f"An error has occured while reading prediction file from p2rank:\n{error}")
        sys.exit(1)
    
    return res_dict

def conv(x):
    
    return x.strip()

def write_pocket_file(ref, res_dict, outdir, query_chain):
    """Write the pocket file

    Args:
        ref (patlib.Path): Path to reference file
        res_dict (dict): Dict containing as key the chain and as values the positions
        outdir (pathlib.Path): Path to the output directory

    Returns:
        output (pathlib.Path): Path of the pocket output file
    """
    
    # Get the path of the pdb file used for p2rank
    pdb = Path(ref.read_text().split("\n")[0])
    # Get the file name without the extension
    pdb_id = pdb.stem
    output = Path.joinpath(outdir, "pocket.csv")
    
    try:
        chain = list(res_dict.keys())[0]
    except:
        logging.error(f"0 results for p2rank, this may be due to an incorrect --chain value : {query_chain}")
        sys.exit(1)
    
    res_str = ''
    for elem in res_dict[chain]:
        res_str += f',{elem}'
    
    text = f"{pdb_id},{chain}{res_str}"
    
    output.write_text(text)
    
    return output

## ------------------------- Homology modeling ------------------------------ ##

def run_build_ali(ref, seq, pocket, outdir, pid, log):
    """Run build_ali.py

    build_ali.py is the script used to prepare the modeling step.

    Args:
        ref (patlib.Path): Path to reference file
        seq (patlib.Path): Path to a multi fasta file
        pocket (patlib.Path): Path to the pocket file
        outdir (pathlib.Path): Path to the output directory
        pid (float): identity cutoff

    Returns:
        ret (subprocess.CompletedProcess): The completed process
    """
    
    # Absolute path of this file
    main_path = Path(__file__).absolute()
    parent_path = main_path.parent
    # Path of the build_ali.py
    src_path = Path.joinpath(parent_path, 'src', "build_ali.py")

    # Run the script    
    command = f"python3 {src_path} -o {outdir} -r {ref} -s {seq} -p {pocket} --id {pid}"
    try:
        if log is None:
            ret = subprocess.run(command.split(), check=True)
        else:
            with open(log, "a") as f_log:
                ret = subprocess.run(command.split(), check=True, stdout=f_log,
                                        stderr=subprocess.STDOUT)
    except Exception as error:
        logging.error(f"An error has occured when lauching the build_ali.py process:\n{error}")
        sys.exit(1)
        
    if ret.returncode != 0:
        logging.error(f"An error has occured during the build_ali.py process:\n"+
                      f"{ret.stderr.decode('utf-8')}")
        sys.exit(1)
    
    return ret

def run_modeling(job, outdir, threads, log):
    """Run modeling.py in parallel with GNU Parallel

    modeling.py is the script used to build model of one target sequence

    Args:
        job (pathlib.Path): Path of file containing the list of inputs
        outdir (pathlib.Path): Path to the output directory
        threads (int): Number of parallel jobs

    Returns:
        ret (subprocess.CompletedProcess): The completed process
    """
    
    # Absolute path of this file
    main_path = Path(__file__).absolute()
    parent_path = main_path.parent
    # Path of the modeling.py
    src_path = Path.joinpath(parent_path, 'src', "modeling.py")
    
    # Create the models directory which will contains the best model of
    # each target (if it pass the identity cutoff)
    model_dir = Path.joinpath(outdir, "models")
    if not model_dir.exists():
        model_dir.mkdir()
    
    # Run modeling.py in parallel
    command = f'parallel -j {threads} python3 {src_path} -o {model_dir} -a :::: {job}'
    
    try:
        if log is None:
            ret = subprocess.run(command.split())
        else:
            with open(log, "a") as f_log:
                ret = subprocess.run(command.split(), stdout=f_log,
                                        stderr=subprocess.STDOUT)
    except Exception as error:
        logging.error(f"An error as occured when launching modeling.py process:\n{error}")
        sys.exit(1)
    
    if ret.returncode == 0:
        tmp_dir = Path.joinpath(outdir, "tmp")
        ali_dir = Path.joinpath(outdir, "ali")
        rm_command = f"rm -r {tmp_dir} {ali_dir}"
        if log is None:
            subprocess.run(rm_command.split())
        else:
            with open(log, "a") as f_log:
                subprocess.run(rm_command.split(), stdout=f_log, stderr=subprocess.STDOUT)
    else:
        logging.error(f"An error has occured during the modeling.py process:\n"+
                      f"{ret.stderr.decode('utf-8')}")
   
    return ret

## ----------------------- Strcutural alignment ----------------------------- ##

def pairwise_alignment(yml, models_file, outdir, threads, log):
    """Runs USalign in parallel with GNU Parallel

    Args:
        yml (dict): a dictionary corresponding to the contents of the yaml file
        models_file (pathlib.Path): Path to the models file
        outdir (pathlib.Path): Path to the output directory
        threads (int): Number of parallel jobs

    Returns:
        pairwise_dir (pathlib.Path): Path to the directory containing all pairwise alignment
    """
    
    USALIGN = yml["usalign"]
    pairwise_dir = Path.joinpath(outdir, "pairwise")
    superposition_dir = Path.joinpath(outdir, "superposition")
    
    if not pairwise_dir.exists():
        pairwise_dir.mkdir()
        
    if not superposition_dir.exists():
        superposition_dir.mkdir()
    
    pair_file = Path.joinpath(outdir, "pair_list.txt")
    text = ""
    
    # Building the file to run gnu parallel
    with open(models_file, "r") as f:
        for i, line in enumerate(f):
            split_line = line.split()
            if not Path(split_line[0]).exists:
                logging.error(f"line {i} in {models_file}: {split_line[0]} "
                              "dosen't exist")
                sys.exit(1)
            elif not Path(split_line[1]).exists():
                logging.error(f"line {i} in {models_file}: {split_line[1]} "
                              "dosen't exist")
                sys.exit(1)
            try:
                ref = Path(split_line[1]).stem
                model = Path(split_line[0]).stem
            except IndexError:
                logging.error(f"'{models_file}' seems to not contains 2 paths on each line:\n"+
                              f"{line.strip()}")
                sys.exit(1)
            except Exception as error:
                logging.error(f"An error has occured while reading '{models_file}':\n"+
                              f"{error}")
                sys.exit(1)
                
            if ref == model:
                continue

            
            output = Path.joinpath(pairwise_dir, f"{model}_-{ref}.fasta")
            super_name = Path.joinpath(superposition_dir, f"{model}")
            text += f"{split_line[1]} {split_line[0]} -o {super_name} -outfmt 1 > {output}\n"
    
    pair_file.write_text(text)
    # Run the parallel command
    command = f"parallel -j {threads} ::: {USALIGN} :::: {pair_file}"
    if log is None:
        ret = subprocess.run(command.split())
    else:
        with open(log, "a") as f_log:
            ret = subprocess.run(command.split(), stdout=f_log, stderr=subprocess.STDOUT)
    pair_file.unlink()
    
    if ret.returncode == 0:
        all_pml = [f for f in superposition_dir.iterdir() if f.match("*.pml")]
        for pml in all_pml:
            pml.unlink()
    
    else:
        logging.error(f"An error has occured during USalign process:\n"+
                      f"{ret.stderr.decode('utf-8')}")
        sys.exit(1)
        
    return pairwise_dir


def renumber_residues(ref_list):
    """Renumbering reference structure

    Args:
        ref_list (list): [pathlib.Path, str(chain), list(positions)]

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
                if line.startswith("ATOM") and line[21:22] == chain and line[17:20] != "HOH":
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
                      f"in the pocket file not found in the {pdb}")
        sys.exit(1)
    
    return renum

def extract_aligned_pos(id_ref, id_model, ref_list, alignment_file, keep_ref):
    """Get positions aligned with the reference pocket

    Args:
        id_ref (str): Reference id
        id_model (str): Model id
        ref_list (list): [pathlib.Path, str(chain), list(positions)]
        alignment_file (pathlib.Path): Path of the alignment file
        keep_ref (bool): Indicate whether we write the reference positions

    Returns:
        text (str): text to add to output file
    """
    
    renum_pos = ref_list[-1]
    pos_str = []
    aln = {id_ref:"", id_model:""}
    text = "" 
    
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
    
    if keep_ref == True:
        text += f">{id_ref}\n"
        j = 0
        pocket = ""
        for i, aa in enumerate(aln[id_ref]):
            if aa != "-":
                if j in renum_pos:
                    pos_str.append(i)
                    pocket += aa
                j += 1
        text += pocket + "\n"
    else:
        j = 0
        pocket = ""
        for i, aa in enumerate(aln[id_ref]):
            if aa != "-":
                if j in renum_pos:
                    pos_str.append(i)
                j += 1
    
    text += f">{id_model}\n"
    try:
        pocket = "".join([aln[id_model][i] for i in pos_str])
    except:
        pocket = ""
    text += pocket + "\n"
    
    return text

def build_multiple_alignment(ref_file, pocket_file, models_file, yml, args, outdir):
    """Build multiple alignment

    Args:
        ref_file (pathlib.Path): Path to reference file
        pocket_file (pathlib.Path): Path to pocket file
        models_file (pathlib.Path): Path to the models file
        yml (dict): a dictionary corresponding to the contents of the yaml file
        args (argparse.Namespace): the object containing all arguments
        outdir (pathlib.Path): Path to the output directory

    Returns:
        taxt (str): Text (multiple alignment) to write in the output
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
    
    pair_start = datetime.datetime.now()
    logging.info(f"Start of Structural Parwise Alignment with US-align")
    pairwise_dir = pairwise_alignment(yml, models_file, outdir, args.threads, args.log)
    logging.info(f"SPA elasped time: {datetime.datetime.now() - pair_start}")
    
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

def search_active_site_in_msa(msa):
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
        logging.error(f"An error has occured while reading '{msa}':\n'{aln}' doesn't exists")
        sys.exit(1)
    
    if not id_file.exists():
        logging.error(f"An erro has occured while reading '{msa}':\n'{id_file}' doesn't exists")
        sys.exit(1)
    
    # get the reference for each sequences
    map_target_ref = {}
    with open(id_file, "r") as f:
        for line in f:
            split_line = line.strip().split()
            map_target_ref[split_line[0]] = split_line[1]
    
    # parse the multiple sequences alignment
    all_seq = {}
    with open(aln, "r") as f:
        seq_id = ""
        for line in f:
            if line.startswith(">"):
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
    
    # for each reference sequence, we extract its active site via the
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

        ref_seq = [s for s in map_target_ref if map_target_ref[s] == ref_id]
        for seq in ref_seq:
            text += f">{seq}\n"
            text += "".join(all_seq[seq][i] for i in ref[ref_id]["aln"])
            text += "\n"
    
    return text

## ------------------------------- Clustering ------------------------------- ##

def read_alignment(file, outdir):
    """Read sequences alignment in fasta format

    Args:
        file (pathlib.Path): Path of the fasta file
        outdir (pathlib.Path): Path to the output directory

    Returns:
        sequences (dict): Dictionnary with the sequence id as key and the corresponding sequence as value
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

    if len(removed) != 0:
        text = "\n".join([f"{seq_id}\t{removed[seq_id]}" for seq_id in removed])
        output = Path.joinpath(outdir, "removed_sequences.txt")
        output.write_text(text)
    
    return sequences

def read_matrix(matrix):
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

def pairwise_score(scoring_dict, seqA, seqB, weighted_pos):
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
                                "arge given the same score as '-' and 'X'")
                
                if i+1 in weighted_pos:
                    score += 20 * 5
                else:
                    score += 5
    
    return score

def dissimilarity(sequences, scoring_dict, weighted_pos):
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

def dbscan_clustering(data, threshold, min_samples, threads):
    
    try:
        dbscan = DBSCAN(eps=threshold, metric="precomputed", n_jobs=threads,
                        min_samples=min_samples)
        
        labels = dbscan.fit_predict(X=data)
    except Exception as error:
        logging.error(f"An error has occured during the clustering:\n{error}")
        sys.exit(1)
    
    return labels

def formatting_output(sequences, key_list, labels):
    
    try:
        G = [(key_list[i], sequences[key_list[i]], n) for i, n in enumerate(labels)]
        G = sorted(G, key=lambda x: x[2])
    except Exception as error:
        logging.error(f"An error has occured durint output formatting:\n{error}")
        sys.exit(1)
    
    return G

## ------------------------------- Weblogo ---------------------------------- ##

def write_fasta(group, fasta):
    
    text = ""
    for elem in group:
        text += f">{elem[0]}\n{elem[1]}\n"
    
    fasta.write_text(text)
    
    return 0

def build_logo(lenght, fasta, outdir, n, prefix, out_format):
    
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
        logging.error(f"An error has occured when creating the logo of {prefix}{n}:\n{error}")
    
    return 0

##########
## MAIN ##
##########

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, metavar="", default="./",
                        help="output directory [default: ./]")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=6,
                        help="number of cpu threads [default: 6]")
    parser.add_argument("-l", "--log", type=str, metavar="",
                        help="log file path, if it's not provied the log are display in the stdout")
    input_opt = parser.add_argument_group("References Structures options")
    input_opt.add_argument("-r","--ref", type=str, metavar="",
                        help="file containing paths to all references")
    input_opt.add_argument("-p", "--pocket", type=str, metavar="",
                        help="file indicating for each reference, the chain and"+
                        " the active site positions. If no file is provided, P2RANK "+
                        "is run to detect pockets")
    input_opt.add_argument("--chain", type=str, metavar="", default="all",
                           help="specifies chains for pocket search, separated "+
                           "by ',' only used if --pocket isn't provided [default: all]")
    targts_opt = parser.add_argument_group("Targets options",
                                        "If --seqs is given, homology modeling"+
                                        " is performed. If --models is given, "+
                                        "homology modeling is not performed "+
                                        "and if --actice-site is given just "+
                                        "the clustering is performed")
    targts_opt_ex = targts_opt.add_mutually_exclusive_group(required=True)

    targts_opt_ex.add_argument("-s","--seqs", type=str, metavar="",
                            help="multi fasta file or directory containing each single fasta file")
    targts_opt_ex.add_argument("-m","--models", type=str, metavar="",
                            help="file containing paths to all models and for each model, his reference")
    targts_opt_ex.add_argument("-M","--msa", type=str, metavar="",
                            help="file indicating active"+
                            " site positions for each references, identity_"+
                            "target_ref path and the path of an MSA")
    targts_opt_ex.add_argument("-a","--active-site", type=str, metavar="",
                               help="active site alignment in fasta format"+
                               ", can be used to create subgroup")
    targts_opt.add_argument("--id", type=float, metavar="", default=30.0,
                            help="percent identity cutoff between target and " +
                            "reference to build a model of the target, only " +
                            "used with -s, --seqs [default: 30.0]")
    dbscan_opt = parser.add_argument_group("Clustering options")
    dbscan_opt.add_argument("-e", "--eps", type=str, metavar="", default="auto",
                            help="maximum distance between two samples for them to be considered neighbors [0,1] [default: auto]")
    dbscan_opt.add_argument("--min-samples", type=str, metavar="", default="auto",
                            help="the number of samples in a neighborhood for a point to be considered as a core point [default: auto]")
    dbscan_opt.add_argument("--dbtest", type=int, choices=[0, 1], default=0,
                            help="0: use the --eps value, 1: test different values")
    dbscan_opt.add_argument('-w','--weighted-pos', type=str, metavar="", default=None,
                            help="pocket position with more weight for clustering"
                            ", positions are numbered from 1 to the total number"
                            " of positions. To give several positions, separate"
                            " them with commas, e.g: 1,6,12")
    weblogo_opt = parser.add_argument_group("Weblogo options")
    weblogo_opt.add_argument("--prefix", type=str, metavar="", default="G",
                             help="prefix for logo title before the cluster [default: G]")
    weblogo_opt.add_argument("--format", type=str, metavar="", default="png",
                             choices=["eps", "png"],
                             help="file format for output logos, 'eps' or 'png' [default: 'png']")
    
    args = parser.parse_args()
    
    start = datetime.datetime.now()

    # Configure logging
    if args.log is None:
        logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(asctime)s - %(message)s")
        
    else:
        logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s - %(asctime)s - %(message)s",
                    filename=args.log, )
    
    # Read the Config file
    yml = read_yaml()
    
    # Make output directory if doesn't exist 
    outdir = Path(args.outdir).absolute()
    if not outdir.exists():
       outdir.mkdir()
    
    if args.ref is not None:
        # check references
        ref_file = Path(args.ref).absolute()
        if not ref_file.exists():
            logging.error(f"{ref_file} doesn't exist")
            sys.exit(1)
        
        ref_list = ref_file.read_text().split()
        for r in ref_list:
            if not Path(r).exists():
                logging.error(f"An error has occured while reading {ref_file}: "
                              f"{r} doesn't exist")
            
        if not args.pocket is None:
            pocket_file = Path(args.pocket).absolute()
            if not pocket_file.exists():
                logging.error(f"{pocket_file} doesn't exist")
                sys.exit(1)
        
        else:            
            prank_output = Path.joinpath(outdir, "prank_output")
            if not prank_output.exists():
                prank_output.mkdir()
                
            ds = build_ds(ref_file, prank_output, args.chain)
            prank_results = run_prank(yml, ds, prank_output)
            pocket_dict = extract_pocket(prank_output)
            pocket_file = write_pocket_file(ref_file, pocket_dict, outdir, args.chain)
    
    elif args.seqs is not None or args.models is not None:
        logging.error(f"argument -r, --ref is required if -s, --seqs or -m, --models is used")
        sys.exit(1)
        
    if not args.seqs is None:
        seq_path = Path(args.seqs).absolute()
        if not seq_path.exists():
            logging.error(f"argument -s/--seqs '{seq_path}' doesn't exist")
            sys.exit(1)
        else:
            
            PID = args.id
            if PID < 0:
                logging.error(f"--id negative value: {PID}")
                sys.exit(1)
                
            ret_build = run_build_ali(ref_file, seq_path, pocket_file, outdir,
                                      PID, args.log)
            job_file = Path.joinpath(outdir, "job_file.txt")

            if not job_file.exists():
                logging.error(f"An error has occurend during the preparation of the homology modeling")
                sys.exit(1)
            else:
                start_model = datetime.datetime.now()
                ret_model = run_modeling(job_file, outdir, args.threads, args.log)
                logging.info(f"modeling duration: {datetime.datetime.now() - start_model}")
    
    if not args.models is None:
        models_file = Path(args.models).absolute()
        if not models_file.exists():
            logging.error(f"argument -m/--models '{models_file}' doesn't exist")
            sys.exit(1)
        
    else:
        models_file = Path.joinpath(outdir, "models.txt")
        
    if args.active_site is None:
        
        if args.msa is None:
        
            text = build_multiple_alignment(ref_file, pocket_file, models_file,
                                            yml, args, outdir)
                
            multiple_alignment = Path.joinpath(outdir, "active_site_alignment.fasta")
            multiple_alignment.write_text(text)
            
        else:
            if Path(args.msa).exists():
                text = search_active_site_in_msa(Path(args.msa))
                multiple_alignment = Path.joinpath(outdir, "active_site_alignment.fasta")
                multiple_alignment.write_text(text)
            else:
                logging.error(f"argument -M/--msa '{args.msa}' doesn't exist")

    else:
        multiple_alignment = Path(args.active_site)
    
    logging.info("Reading Multiple Alignment")
    sequences = read_alignment(multiple_alignment, outdir)
    
    logging.info("Reading Scoring Matrix")
    matrix = Path(yml["distances"])
    scoring_dict = read_matrix(matrix)
    
    if args.weighted_pos is None:
        weighted_pos = []
    else:
        try:
            weighted_pos = [int(x) for x in args.weighted_pos.split(",")]
        except ValueError:
            logging.error(f"-w/--weighted-pos accept only integers separated by ',' e.g: 1,6,12")
            sys.exit(1)
        logging.info(f"Weighted positions: {weighted_pos}")
    
    logging.info("Compute Dissimilarities")
    key_list, data = dissimilarity(sequences, scoring_dict, weighted_pos)
    perc = np.percentile(data, [25, 50, 75])
    
    logging.info(f"q1\tmed\tq3\tmean")
    logging.info(f"{perc[0]:.3f}\t{perc[1]:.3f}\t{perc[2]:.3f}\t{data.mean():.3f}")
    
    if args.dbtest == 1:
        eps_list = [0.3, 0.2, 0.1, round(perc[0], 2),
                    round(perc[0] - (perc[0] * 0.1), 2),
                    round(perc[0] - (perc[0] * 0.15), 2),
                    round(perc[0] - (perc[0] * 0.20), 2),
                    round(perc[0] - (perc[0] * 0.25), 2)]
            
    else:
        if args.eps != "auto":
            try:
                eps_list = [float(args.eps)]
            except:
                logging.error(f"argument -e, --eps invalid value : {args.eps}")
                sys.exit(1)
        else:
            eps_list = [round(perc[0] - (perc[0] * 0.1), 2)]
                
    if args.min_samples == "auto":
        if len(sequences) <= 1500:
            min_samples = 5
        else:
            min_samples = 25
    else:
        try:
            min_samples = int(args.min_samples)
        except:
            logging.error(f"argument --min-samples invalid value : {args.min_samples}")
            sys.exit(1)
                
    str_eps_list = ["q1", "q1-10p", "q1-15p", "q1-20p", "q1-25p"]
    for i, eps in enumerate(eps_list):
        if i <= 2:
            str_eps = str(eps)
            str_eps.replace(".", "_")
        else:
            str_eps = str_eps_list[0]
            del str_eps_list[0]
            
        logging.info(f"eps: {eps}\tmin_samples: {min_samples}")
        labels = dbscan_clustering(data=data, threshold=eps,
                                    min_samples=min_samples,
                                    threads=args.threads)
            
        unique, count = np.unique(labels, return_counts=True)
        logging.info(f"Number of clusters: {len(unique)}")
        logging.info({a:b for a, b in zip(unique, count)})
            
        try:
            score = silhouette_score(X=data, labels=labels, metric="precomputed")
            logging.info(f"silhouette score: {score:.3f}")
        except:
            logging.info("silhouette score: -")
                
        G = formatting_output(sequences, key_list, labels)
            
        if len(eps_list) <= 1:
            dbscan_output = Path.joinpath(outdir, f"groups_{str_eps}_min_{min_samples}.tsv")
        else:
            outdir = Path.joinpath(outdir, f"eps_{str_eps}_min_{min_samples}")
            if not outdir.exists():
                outdir.mkdir()
                
            dbscan_output = Path.joinpath(outdir, f"groups_{str_eps}_min_{min_samples}.tsv")
                
        with dbscan_output.open(mode="w") as f:
            for elem in G:
                f.write(f"{elem[0]}\t{elem[1]}\t{elem[2]}\n")
    
        for n in unique:
            group_seq = [elem for elem in G if elem[-1] == n]
            fasta = Path.joinpath(outdir, f"G{n}.fasta")
            write_fasta(group=group_seq, fasta=fasta)
            build_logo(len(group_seq), fasta, outdir, n, args.prefix, args.format)
            
        outdir = Path(args.outdir).absolute() 
        
    logging.info(f"Total Elapsed time: {datetime.datetime.now() -  start}")