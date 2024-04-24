import argparse
import re
from pathlib import Path

import modeller

###############
## Functions ##
###############

def multi_to_single(file, outdir):
    """Transform multi fasta file into single fasta files

    Replaces potential '|', '.' and ':' in the sequence id to avoid file
    and command parsing problems.
    
    For uniprot, we keep only the unique identifier between the pipes.

    Args:
        file (pathlib.Path): Path of the multi fasta file
        outdir (pathlib.Path): Path of the output directory

    Returns:
        tmp_dir (pathlib.Path): Path of the temporary directory where new files are stored
    """
    
    sequences = {}

    with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                
                l = line
                
                if "tr|" in line or "sp|" in l:
                    seq_id = re.search("\\|(\\w+)\\|", l).group(1)
                    sequences[seq_id] = f">{seq_id} {''.join(line.split()[1:])}\n"
                else:
                    seq_id = l.split()[0][1:]
                    for c in [".", ":", "|"]:
                        seq_id = seq_id.replace(c, "_")
                    sequences[seq_id] = f">{seq_id} {''.join(line.split()[1:])}\n"
            else:
                sequences[seq_id] += line
                       
    tmp_dir = Path.joinpath(outdir, "tmp")
    if not tmp_dir.exists():
        tmp_dir.mkdir()
        
    for key in sequences:
        single_fasta = Path.joinpath(tmp_dir,
                                    f"{key}.fasta")
        
        single_fasta.write_text(sequences[key])
    
    return tmp_dir

def Build_ref_target_ali(fasta, ref, pocket, outdir, PID):
    """Build alignment for the modeling

    Args:
        fasta (pathlib.Path): Path of the single fasta file
        ref (pathlib.Path): Path of the file containing the references path
        pocket (pathlib.Path): Path of the file containing the pocket chain and position
        outdir (pathlib.Path): Path of the output directory
        PID (int): cutoff of %id

    Returns:
        (tuple(int, float, str)): (return code, %id, path of the ref with the best %id)
    """
    
    # get the chain to use for modeling for each reference
    ref_chain = {}
    with open(pocket, "r") as f:
        for line in f:
            split_line = line.split(',')
            ref_chain[split_line[0]] = split_line[1]
    
    # create modeller objects
    env = modeller.Environ()
    modeller.log.level(output=0, notes=0, warnings=0, errors=1, memory=0)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    aln = modeller.Alignment(env)
    template = modeller.Model(env)
    for r in ref:
        # case of references are in the set of sequences
        if fasta.stem == Path(r).stem:
            return (1, 0.0, "")
        elif r == "":
            continue
        else:
            # use the chain specified in the pocket file
            segment = (f"FIRST:{ref_chain[Path(r).stem]}",
                       f"LAST:{ref_chain[Path(r).stem]}")
            # read the structure and add it in the alignment
            template.read(file=r,model_format="PDB", model_segment=segment)
            aln.append_model(template, align_codes=Path(r).stem, atom_files=r)

    # add the sequence in the alignment
    aln.append(file=str(fasta), alignment_format="FASTA", align_codes=fasta.stem)

    # One or more references
    if len(ref_chain) > 1:
        aln.salign(gap_function=True,
                feature_weights=(1.0,0,0,0,0,0),
                max_gap_length=20, overhang=0, similarity_flag=True,
                gap_penalties_1d=(-450, 0),
                gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.)
                )
    else:
        aln.salign(gap_function=True, alignment_type='PAIRWISE',
                feature_weights=(1.0,0,0,0,0,0),
                max_gap_length=20, overhang=0, similarity_flag=True,
                gap_penalties_1d=(-450, 0),
                gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.)
                )
    
    
    # Identify the reference with the best %id with the target sequence
    pid = 0
    best_ref =""
    for seq in aln[:-1]:
        x = seq.get_sequence_identity(aln[-1])
        if x > pid:
            pid = x
            best_ref = seq._Sequence__get_atom_file()
    
    # condition to build model for the target
    if pid >= PID:
        models_text = f"{Path.joinpath(outdir.parent, 'models', f'{fasta.stem}.pdb')} {best_ref}\n"
        models_file = Path.joinpath(outdir.parent, "models.txt")
        
        # create or complete the model file
        if not models_file.exists():
            models_file.write_text(models_text)
        else:
            with models_file.open(mode="a") as f:
                f.write(models_text)
        
        # write the alignment file
        output = Path.joinpath(outdir, f"{fasta.stem}.ali")
        aln.write(file=str(output), alignment_format="PIR")
    
        # create or complete job file
        job_file = Path.joinpath(outdir.parent, "job_file.txt")
        ref_text = "+".join([code for code in ref_chain])
        text = f"{output}+{ref_text}\n"
        if job_file.exists():
            with job_file.open(mode="a") as f:
                f.write(text)
        else:
            job_file.write_text(text)
    
        return (0, pid, best_ref)
    
    else:
        return (2, pid, best_ref)

##########
## Main ##
##########

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--ref", type=str, metavar="",
                        help="file containing paths to all references and the chain for modeling")
    parser.add_argument("-s", "--seq", type=str, metavar="",
                        help="multi fasta file or directory containing all single fasta")
    parser.add_argument("-p", "--pocket", type=str, metavar="",
                        help="file indicating for each reference, the chain and the pocket positions")
    parser.add_argument("-o", "--outdir", type=str, metavar="", default="./",
                        help="output directory [default: ./]")
    parser.add_argument("--id", type=float, metavar="",
                        help=r"%id cutoff between target and reference to build a model of the target")
    
    args = parser.parse_args()
    
    PID = args.id
    
    outdir = Path(args.outdir).absolute()
    if not outdir.exists():
       outdir.mkdir()
       
    ref_file = Path(args.ref)
    pocekt_file = Path(args.pocket)
    seq_path = Path(args.seq)
    if seq_path.is_dir():
        all_fasta = [f for f in seq_path.iterdir()]
    else:
        tmp_dir = multi_to_single(file=seq_path, outdir=outdir)
        all_fasta = [f for f in tmp_dir.iterdir()]

    all_ref = ref_file.read_text().split("\n")
    
    ali_dir = Path.joinpath(outdir, "ali")
    #ali_dir = outdir
    if not ali_dir.exists():
        ali_dir.mkdir()
        
    for i, fasta in enumerate(all_fasta):

        aln = Build_ref_target_ali(fasta, all_ref, pocekt_file, ali_dir,
                                   PID)
        
        identity_file = Path.joinpath(outdir, "identity_targets_refs.tsv")
        if not identity_file.exists():
            identity_file.write_text(f"{fasta.stem}\t{Path(aln[2]).stem}\t{aln[1]:.2f}\t{PID:.2f}\n")
        else:
            with identity_file.open(mode="a") as f:
                f.write(f"{fasta.stem}\t{Path(aln[2]).stem}\t{aln[1]:.2f}\t{PID:.2f}\n")