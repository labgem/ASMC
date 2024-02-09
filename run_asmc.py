import sys
import argparse
import datetime
import logging
import subprocess
import yaml
from pathlib import Path

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
            
    logging.info(f"Total Elapsed time: {datetime.datetime.now() -  start}")