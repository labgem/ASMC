import sys
import argparse
import datetime
import logging
import subprocess
from itertools import groupby
from pathlib import Path

import yaml
import numpy as np

###############
## Functions ##
###############

def read_yaml():
    """Read the yaml file
    
    Load the content of the yaml file and check the validity of paths.
    The configuration should be place in a directory named data in the same
    location as run_asmc.py, e.g :

    .
    ├── data
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
    data_path = Path.joinpath(parent_path, "data")
    
    config_path = Path.joinpath(data_path, "config_asmc.yml")

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
    
    
    prediction = [f for f in outdir.iterdir() if f.match("*predictions*")][0]
    if len(prediction) == 0:
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
    input_opt = parser.add_argument_group("References options")
    input_opt.add_argument("-r","--ref", type=str, metavar="",
                        help="file containing paths to all references")
    input_opt.add_argument("-p", "--pocket", type=str, metavar="",
                        help="file indicating for each reference, the chain and"+
                        " the pocket positions. If no file is provided, P2RANK "+
                        "is run to detect pockets")
    input_opt.add_argument("--chain", type=str, metavar="", default="all",
                           help="Specifies chains for pocket search, separated "+
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
    targts_opt_ex.add_argument("-a","--active-site", type=str, metavar="",
                               help="active site alignment in fasta format")
    targts_opt.add_argument("--id", type=float, metavar="", default=30.0,
                            help="percent identity cutoff between target and " +
                            "reference to build a model of the target, only " +
                            "used with -s, --seqs [default: 30.0]")
    
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
        ref_file = Path(args.ref).absolute()
        if not ref_file.exists():
            logging.error(f"{ref_file} doesn't exist")
            sys.exit(f"{ref_file} doesn't exist")
            
            
        if not args.pocket is None:
            pocket_file = Path(args.pocket).absolute()
            if not pocket_file.exists():
                logging.error(f"{pocket_file} doesn't file")
                sys.exit(f"{pocket_file} doesn't file")
        
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
                sys.exit(f"--id negative value: {PID}")
                
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
     
    logging.info(f"Total Elapsed time: {datetime.datetime.now() -  start}")