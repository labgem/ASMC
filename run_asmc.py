import asmc
import sys
import yaml
import argparse
import subprocess
import logging
import datetime
import numpy as np
from pathlib import Path
from sklearn.metrics import silhouette_score

###############
## Functions ##
###############

def read_yaml(args):
    """Read the yaml file
    
    Load the content of the yaml file and check the validity of paths.
    The configuration should be place in a directory named ressources in the same
    location as run_asmc.py, e.g :

    .
    ├── ressources
    │   ├── AA_distances.tsv
    │   └── config_asmc.yml
    ├── run_asmc.py
    └── ascm
        ├── __inti__.py
        ├── asmc.py
        └── ...


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
        
        elif key == "usalign" and (args.msa is None and args.active_site is None):
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
                
        elif key == "java" and args.ref is not None and args.pocket is None:
            command = f"{yml[key]} --version"
            try:
                ret = subprocess.run(command.split(), capture_output=True)
            except Exception as error:
                logging.error(f"An error has occured with java:\n{error}")
                sys.exit(1)
                
            if ret.returncode != 0:
                logging.error(f"{ret.stderr.decode('utf-8')}")
                sys.exit(1)
    
    return yml

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
            logging.error("An error has occured during the p2rank process:\n"
                          f"{result.stderr.decode('utf-8')}")
            sys.exit(1)
    
    return result

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
    src_path = Path.joinpath(parent_path, 'asmc', "build_ali.py")

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
        logging.error(f"An error has occured when lauching the build_ali.py "
                      f"process:\n{error}")
        sys.exit(1)
        
    if ret.returncode != 0:
        logging.error(f"An error has occured during the build_ali.py process:\n"
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
    src_path = Path.joinpath(parent_path, 'asmc', "modeling.py")
    
    # Create the models directory which will contains the best model of
    # each target (if it pass the identity cutoff)
    model_dir = Path.joinpath(outdir, "models")
    if not model_dir.exists():
        model_dir.mkdir()
    
    # Run modeling.py with gnu parallel
    command = f'parallel -j {threads} python3 {src_path} -o {model_dir} -a :::: {job}'
    
    try:
        if log is None:
            ret = subprocess.run(command.split())
        else:
            with open(log, "a") as f_log:
                ret = subprocess.run(command.split(), stdout=f_log,
                                        stderr=subprocess.STDOUT)
    except Exception as error:
        logging.error(f"An error as occured when launching modeling.py "
                      f"process:\n{error}")
        sys.exit(1)
    
    if ret.returncode == 0:
        tmp_dir = Path.joinpath(outdir, "tmp")
        ali_dir = Path.joinpath(outdir, "ali")
        rm_command = f"rm -r {tmp_dir} {ali_dir}"
        if log is None:
            subprocess.run(rm_command.split())
        else:
            with open(log, "a") as f_log:
                subprocess.run(rm_command.split(), stdout=f_log,
                               stderr=subprocess.STDOUT)
    else:
        logging.error(f"An error has occured during the modeling.py process:\n"
                      f"{ret.stderr.decode('utf-8')}")
   
    return ret

def pairwise_alignment(yml, models_file, outdir, threads, log):
    """Runs USalign in parallel with GNU Parallel

    Args:
        yml (dict): a dictionary corresponding to the contents of the yaml file
        models_file (pathlib.Path): Path to the models file
        outdir (pathlib.Path): Path to the output directory
        threads (int): Number of parallel jobs
        log (pathlib.Path): Path to the log file

    Returns:
        pairwise_dir (pathlib.Path): Path to the directory containing all pairwise alignment
    """
    
    # Name or path of Usalign binaries
    USALIGN = yml["usalign"]
    # Usalign output directories
    pairwise_dir = Path.joinpath(outdir, "pairwise")
    superposition_dir = Path.joinpath(outdir, "superposition")
    
    if not pairwise_dir.exists():
        pairwise_dir.mkdir()
        
    if not superposition_dir.exists():
        superposition_dir.mkdir()
    
    pair_file = Path.joinpath(outdir, "pair_list.txt")
    text = ""
    
    # Building the file to be used with gnu parallel
    with open(models_file, "r") as f:
        for i, line in enumerate(f):
            split_line = line.split()
            if not Path(split_line[0]).exists:
                logging.error(f"line {i} in {models_file}: {split_line[0]} "
                              "doesn't exist")
                sys.exit(1)
            elif not Path(split_line[1]).exists():
                logging.error(f"line {i} in {models_file}: {split_line[1]} "
                              "doesn't exist")
                sys.exit(1)
            try:
                ref = Path(split_line[1]).stem
                model = Path(split_line[0]).stem
            except IndexError:
                logging.error(f"'{models_file}' seems to not contains 2 paths" 
                              f"on each line:\n {line.strip()}")
                sys.exit(1)
            except Exception as error:
                logging.error(f"An error has occured while reading"
                              f"'{models_file}':\n{error}")
                sys.exit(1)
                
            if ref == model:
                continue

            output = Path.joinpath(pairwise_dir, f"{model}_-{ref}.fasta")
            super_name = Path.joinpath(superposition_dir, f"{model}")
            text += f"{split_line[1]} {split_line[0]} -o {super_name} -outfmt 1"
            text += f" > {output}\n"
    
    pair_file.write_text(text)
    
    # Run the parallel command
    command = f"parallel -j {threads} ::: {USALIGN} :::: {pair_file}"
    if log is None:
        ret = subprocess.run(command.split())
    else:
        with open(log, "a") as f_log:
            ret = subprocess.run(command.split(), stdout=f_log,
                                 stderr=subprocess.STDOUT)
    pair_file.unlink()
    
    # Remove pymol scripts
    if ret.returncode == 0 or ret.returncode == 101:
        all_pml = [f for f in superposition_dir.iterdir() if f.match("*.pml")]
        for pml in all_pml:
            pml.unlink()
    
    else:
        logging.error(f"{ret.returncode}")
        logging.error("An error has occured during USalign process:\n"
                      f"{ret.stderr.decode('utf-8')}")
        sys.exit(1)
        
    return pairwise_dir

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
                        help="log file path, if it's not provied the log are"
                        " display in the stdout")
    input_opt = parser.add_argument_group("References Structures options")
    input_opt.add_argument("-r","--ref", type=str, metavar="",
                        help="file containing paths to all references")
    input_opt.add_argument("-p", "--pocket", type=str, metavar="",
                        help="file indicating for each reference, the chain and"
                        " the active site positions. If no file is provided, "
                        "P2RANK is run to detect pockets")
    input_opt.add_argument("--chain", type=str, metavar="", default="all",
                           help="specifies chains for pocket search, separated "
                           "by ',' only used if --pocket isn't provided "
                           "[default: all]")
    targts_opt = parser.add_argument_group("Targets options",
                                        "If --seqs is given, homology modeling"
                                        " is performed. If --models is given, "
                                        "homology modeling is not performed "
                                        "and if --actice-site is given just "
                                        "the clustering is performed")
    targts_opt_ex = targts_opt.add_mutually_exclusive_group(required=True)

    targts_opt_ex.add_argument("-s","--seqs", type=str, metavar="",
                            help="multi fasta file or directory containing each"
                            " single fasta file")
    targts_opt_ex.add_argument("-m","--models", type=str, metavar="",
                            help="file containing paths to all models and for "
                            "each model, his reference")
    targts_opt_ex.add_argument("-M","--msa", type=str, metavar="",
                            help="file indicating active"+
                            " site positions for each references, identity_"
                            "target_ref path and the path of an MSA")
    targts_opt_ex.add_argument("-a","--active-site", type=str, metavar="",
                               help="active site alignment in fasta format"
                               ", can be used to create subgroup")
    targts_opt.add_argument("--id", type=float, metavar="", default=30.0,
                            help="percent identity cutoff between target and " 
                            "reference to build a model of the target, only " 
                            "used with -s, --seqs [default: 30.0]")
    dbscan_opt = parser.add_argument_group("Clustering options")
    dbscan_opt.add_argument("-e", "--eps", type=str, metavar="", default="auto",
                            help="maximum distance between two samples for them"
                            " to be considered neighbors [0,1] [default: auto]")
    dbscan_opt.add_argument("--min-samples", type=str, metavar="", default="auto",
                            help="the number of samples in a neighborhood for "
                            "a point to be considered as a core point "
                            "[default: auto]")
    dbscan_opt.add_argument("--test", type=int, choices=[0, 1], default=0,
                            help="0: use the --eps value, 1: test different "
                            "values")
    dbscan_opt.add_argument('-w','--weighted-pos', type=str, metavar="",
                            default=None,
                            help="pocket position with more weight for clustering"
                            ", positions are numbered from 1 to the total number"
                            " of positions. To give several positions, separate"
                            " them with commas, e.g: 1,6,12")
    weblogo_opt = parser.add_argument_group("Weblogo options")
    weblogo_opt.add_argument("--prefix", type=str, metavar="", default="G",
                             help="prefix for logo title before the cluster"
                             " [default: G]")
    weblogo_opt.add_argument("--format", type=str, metavar="", default="png",
                             choices=["eps", "png"],
                             help="file format for output logos, 'eps' or 'png'"
                             " [default: 'png']")
    
    args = parser.parse_args()
    
    start = datetime.datetime.now()

    # Configure logging
    if args.log is None:
        logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(asctime)s - %(message)s")
        
    else:
        logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s - %(asctime)s - %(message)s",
                    filename=args.log)
        
    # Read the Config file
    yml = read_yaml(args)
    
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
                sys.exit(1)
            
        if not args.pocket is None:
            pocket_file = Path(args.pocket).absolute()
            if not pocket_file.exists():
                logging.error(f"{pocket_file} doesn't exist")
                sys.exit(1)
        
        else:            
            prank_output = Path.joinpath(outdir, "prank_output")
            if not prank_output.exists():
                prank_output.mkdir()
                
            ds, ds_text = asmc.build_ds(ref_file, prank_output, args.chain)
            ds.write_text(ds_text)
            prank_results = run_prank(yml, ds, prank_output)
            pocket_dict = asmc.extract_pocket(prank_output)
            pocket_file, pocket_text = asmc.build_pocket_text(ref_file,
                                                              pocket_dict,
                                                              outdir,
                                                              args.chain)
            pocket_file.write_text(pocket_text)
    
    elif args.seqs is not None or args.models is not None:
        logging.error(f"argument -r, --ref is required if -s, --seqs or -m, "
                      "--models is used")
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
                logging.error(f"An error has occurend during the preparation "
                              "of the homology modeling")
                sys.exit(1)
            else:
                start_model = datetime.datetime.now()
                ret_model = run_modeling(job_file, outdir, args.threads, args.log)
                logging.info("modeling duration: "
                             f"{datetime.datetime.now() - start_model}")
    
    if not args.models is None:
        models_file = Path(args.models).absolute()
        if not models_file.exists():
            logging.error(f"argument -m/--models '{models_file}' doesn't exist")
            sys.exit(1)
        
    else:
        models_file = Path.joinpath(outdir, "models.txt")
        
    if args.active_site is None:
        
        if args.msa is None:
            
            pair_start = datetime.datetime.now()
            logging.info(f"Start of Structural Parwise Alignment with US-align")
            pairwise_dir = pairwise_alignment(yml=yml,
                                              models_file=models_file,
                                              outdir=outdir,
                                              threads=args.threads,
                                              log=args.log)
            logging.info(f"SPA elasped time: {datetime.datetime.now() - pair_start}")
            
            text = asmc.build_multiple_alignment(pairwise_dir, ref_file,
                                                 pocket_file)
                
            multiple_alignment = Path.joinpath(outdir,
                                               "active_site_alignment.fasta")
            multiple_alignment.write_text(text)
            
        else:
            if Path(args.msa).exists():
                text = asmc.search_active_site_in_msa(Path(args.msa))
                multiple_alignment = Path.joinpath(outdir,
                                                   "active_site_alignment.fasta")
                multiple_alignment.write_text(text)
            else:
                logging.error(f"argument -M/--msa '{args.msa}' doesn't exist")

    else:
        multiple_alignment = Path(args.active_site)
    
    logging.info("Reading Multiple Alignment")
    sequences, removed = asmc.read_alignment(multiple_alignment)
    if len(removed) != 0:
        text = "\n".join([f"{seq_id}\t{removed[seq_id]}" for seq_id in removed])
        output = Path.joinpath(outdir, "removed_sequences.txt")
        output.write_text(text)
    
    logging.info("Reading Scoring Matrix")
    matrix = Path(yml["distances"])
    scoring_dict = asmc.read_matrix(matrix)
    
    if args.weighted_pos is None:
        weighted_pos = []
    else:
        try:
            weighted_pos = [int(x) for x in args.weighted_pos.split(",")]
        except ValueError:
            logging.error(f"-w/--weighted-pos accept only integers separated "
                          "by ',' e.g: 1,6,12")
            sys.exit(1)
        logging.info(f"Weighted positions: {weighted_pos}")
    
    logging.info("Compute Dissimilarities")
    key_list, data = asmc.dissimilarity(sequences, scoring_dict, weighted_pos)
    perc = np.percentile(data, [25, 50, 75])
    
    logging.info(f"q1\tmed\tq3\tmean")
    logging.info(f"{perc[0]:.3f}\t{perc[1]:.3f}\t{perc[2]:.3f}\t{data.mean():.3f}")
    
    if args.test == 1:
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
            logging.error("argument --min-samples invalid value : "
                          f"{args.min_samples}")
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
        labels = asmc.dbscan_clustering(data=data, threshold=eps,
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
                
        G = asmc.formatting_output(sequences, key_list, labels)
            
        if len(eps_list) <= 1:
            dbscan_output = Path.joinpath(outdir,
                                          f"groups_{str_eps}_min_{min_samples}.tsv")
        else:
            outdir = Path.joinpath(outdir, f"eps_{str_eps}_min_{min_samples}")
            if not outdir.exists():
                outdir.mkdir()
                
            dbscan_output = Path.joinpath(outdir,
                                          f"groups_{str_eps}_min_{min_samples}.tsv")
                
        with dbscan_output.open(mode="w") as f:
            for elem in G:
                f.write(f"{elem[0]}\t{elem[1]}\t{elem[2]}\n")
    
        for n in unique:
            group_seq = [elem for elem in G if elem[-1] == n]
            fasta = Path.joinpath(outdir, f"G{n}.fasta")
            fasta_text = asmc.build_fasta(group=group_seq)
            fasta.write_text(fasta_text)
            asmc.build_logo(len(group_seq), fasta, outdir, n, args.prefix,
                            args.format)
        
        asmc.merge_png(outdir, len(unique), args.prefix)  
        outdir = Path(args.outdir).absolute() 
        
    logging.info(f"Total Elapsed time: {datetime.datetime.now() -  start}")