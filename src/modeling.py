import os
import argparse
from pathlib import Path

import modeller
from modeller import automodel

###############
## Functions ##
###############

def modeling(aln, outdir):
    """Build model

    Args:
        aln (str): String containing the path to the alignmentand the reference name
        outdir (pathlib.Path): Path to the output directory 

    Returns:
        (int): 0
    """

    tmp_dir = Path.joinpath(outdir.parent, 'tmp')
    os.chdir(str(tmp_dir))
    
    env = modeller.Environ()
    modeller.log.level(output=0, notes=0, warnings=0, errors=1, memory=0)
    aln_file = Path(aln.split("+")[0])
    seq = aln_file.stem
    ref = aln.split("+")[1:]
    
    a = automodel.AutoModel(env, alnfile=str(aln_file),
                            knowns=ref, sequence=seq,
                            assess_methods=(automodel.assess.DOPE,
                                            automodel.assess.GA341))
    
    a.starting_model = 1
    a.ending_model = 2
    a.make()
    
    model_name = ""
    dope = 9999
    for model in a.outputs:
        if model["DOPE score"] < dope:
            dope = model["DOPE score"]
            model_name = model["name"]
    
    mv = Path.joinpath(outdir, f"{model_name.split('.')[0]}.pdb")
    Path.joinpath(tmp_dir, model_name).rename(mv)
    
    return 0

##########
## Main ##
##########

if __name__ == "__main__":
    

    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--ali", type=str, metavar="",
                        help="alignment between references and target")
    parser.add_argument("-o", "--outdir", type=str, metavar="", default="./",
                        help="output directory [default: ./]")
    
    args = parser.parse_args()
        
    modeling(args.ali, Path(args.outdir))