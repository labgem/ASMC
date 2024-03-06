from pymol import cmd
from pathlib import Path

def find_ref(target_name, identity_file):
    
    with open(identity_file, "r") as f:
        for line in f:
            if line.startswith(str(target_name)):
                return line.split("\t")[1].strip()
            
    return 1
    
def get_active_site_pos(ref_name, active_site_file):
    
    with open(active_site_file, "r") as f:
        for line in f:
            if line.startswith(ref_name):
                split_line = line.strip().split(',')
                chain, pos = split_line[1], split_line[2:]
                return chain, pos
            
    return 1

def target(target):
    
    working_directory = Path.cwd().absolute()
    
    # Set the path of identity_target_ref.csv
    id_target_ref_path = Path.joinpath(working_directory,
                                       "identity_target_ref.tsv")
    
    # Get the reference name for the target
    ref_name = find_ref(target, id_target_ref_path)
    # Set the reference structure path
    ref_path = Path.joinpath(working_directory,
                             "superposition",
                             target+".pdb")
    
    # Set the target structure path
    target_path = Path.joinpath(working_directory,
                                "models",
                                target+".pdb")
    
    # Load the structures
    cmd.load(ref_path)
    ref_name = f"{target}-{ref_name}"
    cmd.set_name(target, ref_name)
    
    cmd.load(target_path)
    
def active_site(active_site_file):
    
    active_site_path = Path(active_site_file).absolute()
    
    # Get the object list
    object_list = cmd.get_object_list('(all)')
    # Get the name of teh reference structure object
    object_ref_name = [elem for elem in object_list if "-" in elem][0]
    ref_name = object_ref_name.split('-')[-1]
    target_name = object_ref_name.split("-")[0]
    
    # Read the pocket/active_site_file
    chain, pos = get_active_site_pos(ref_name, active_site_path)
    
    # Create reference active site ofbject
    string_selection = f"model {object_ref_name} and chain {chain} and resi {'+'.join(pos)}"
    ref_active_site_name = f"{ref_name}_active_site"
    cmd.create(ref_active_site_name, string_selection)
    cmd.show_as("sticks", ref_active_site_name)

    # Get the residues list of the reference active site
    space = {"ref":[]}
    cmd.iterate(f"{ref_active_site_name} and name CA", "ref.append((resi, resn))",
                space=space)
    
    target_active_site_name = f"{target_name}_active_site"
    print("Ref - Target")
    
    # For each residue in the reference active site
    for i, (resi, resn) in enumerate(space["ref"]):
        # Search for the target residue aligned with it
        string_selection = f"{target_name} and br. all and name CA near_to 1.0 of {ref_active_site_name} and resi {resi} and name CA"
        inner_space = {"res":[]}
        if i == 0:
            cmd.create(target_active_site_name, string_selection)
            cmd.iterate(f"{target_active_site_name} and name CA", "res.append((resi, resn))",
                        space=inner_space)
            try:
                print(f"{resi} {resn} - {inner_space['res'][0][0]}"
                    f" {inner_space['res'][0][1]}")
            except IndexError:
                print(f"{resi} {resn} - Gap")
                continue
        else:
            cmd.select("tmp", string_selection)
            cmd.iterate(f"tmp and name CA", "res.append((resi, resn))",
                        space=inner_space)
            try:
                print(f"{resi} {resn} - {inner_space['res'][0][0]}"
                    f" {inner_space['res'][0][1]}")
            except IndexError:
                print(f"{resi} {resn} - Gap")
                continue
            cmd.copy_to(target_active_site_name, "tmp")
            cmd.delete("tmp")
         
    cmd.show_as("sticks", target_active_site_name)
    
    cmd.center(ref_active_site_name)
    cmd.orient(ref_active_site_name)
  
    # set cartoon transparency
    cmd.set("cartoon_transparency", "0.6", "(all)")
    
cmd.extend(target)
cmd.extend(active_site)