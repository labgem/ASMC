<!-- ![ASMC logo](docs/asmc.png) -->
<!-- center and resize the image with html syntax -->
<p align="center">
  <img src="docs/asmc.png" alt="ASMC logo" />
</p>

# ASMC

ASMC combines (i) homology modeling of family members (MODELLER), (ii) ligand-binding pocket search (P2RANK), (iii) structural alignment of modeled active sites (USalign) and (iv) density-based spatial clustering of obtained alignments (DBSCAN) in a single command line. Clustering step can be carried out on either structural or sequence alignment.

<!-- ![ASMC workflow](docs/ASMC_workflow.svg) -->
<!-- center and resize the image with html syntax -->
<p align="center">
  <img src="docs/ASMC_workflow.svg" alt="ASMC workflow" />
</p>

## Installation

### Download

Download the latest GitHub release to obtain the code: [https://github.com/labgem/ASMC/releases](https://github.com/labgem/ASMC/releases)

### Python requirements

- Python ≥ 3.8
- biopython ≥ 1.81
- numpy
- scikit-learn
- pyyaml
- pillow
- weblogo

You can install the python dependencies with `pip`, `conda` or `mamba` with the following commands (in the ASMC directory) and the files given in the releases:

**pip**

```
pip install ./
```

**conda** or **mamba**
```
conda env create -n env_name -f env.yml
```

The command `pip install ./`, create the `asmc` command and adds it the User PATH environment variable. We therefore recommend that you also use this command after installation via conda or mamba. 

It's recommended to use the `pip install ./` command after the installation via conda or mamba. Indeed, this creates `asmc` command and adds it to the User PATH envrionment variable.

Installation via conda and mamba includes the MODELLER installation, but you still need to request the license key.

### External Software dependencies

- P2RANK - for ligand-binding pocket detection ([https://github.com/rdk/p2rank](https://github.com/rdk/p2rank))
- MODELLER - for homology modeling ([https://salilab.org/modeller/](https://salilab.org/modeller/))
- USalign - for structural alignment ([https://zhanggroup.org/US-align/](https://zhanggroup.org/US-align/))

 Please follow the links above and the instructions given by their authors.

## Configuration

In `ASMC/resources`, add a file exactly named `config_asmc.yml`. This file should contain 3 information:

- the path to the `ASMC/resources/AA_distances.tsv`
- the path of P2RANK executable (or it's name if the location is in the PATH)
- the path of USalign executable (or it's name if the location is in the PATH)

Example:
```yaml
distances: "/home/User/ASMC/resources/AA_distances.tsv"
usalign: "USalign"  # location in the PATH e.g: /usr/bin ...
p2rank: "/home/User/p2rank_2.4.1/prank"  # another location
```

The keys should be identical to this example.

## Quick Usage

Run ASMC in a blind way (unknown active site) using a multi fasta file that should contain at least 100 sequences for clustering to be sufficiently relevant.

```
asmc run --log run_asmc.log --threads 6 -r reference_file -s sequences.fasta
```

`reference_file` should contains the path to the reference(s) structure(s), e.g:
```
<path>/RefA.pdb
<path>/RefB.pdb
```

NB: For more details, see the [wiki](https://github.com/labgem/ASMC/wiki/Options-and-Usages)
