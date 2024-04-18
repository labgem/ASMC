# ASMC

ASMC combines (i) homology modeling of family members (MODELLER), (ii) ligand-binding pocket search (P2RANK), (iii) structural alignment of modeled active sites (USalign) and (iv) density-based spatial clustering of obtained alignments (DBSCAN) in a single command line. Clustering step can be carried out on either structural or sequence alignment.

## Requirements

Requires Python>=3.8

Python dependencies:
- numpy
- scikit-learn
- pyyaml
- pillow
- biopython>=1.81
- weblogo

External dependencies:
- p2rank - for ligand-binding pocket detection ([https://github.com/rdk/p2rank](https://github.com/rdk/p2rank))
- modeller - for homology modeling ([https://salilab.org/modeller/](https://salilab.org/modeller/))
- USalign - for structural alignment ([https://zhanggroup.org/US-align/](https://zhanggroup.org/US-align/))

## Installation

### Download

Download the latest GitHub release to obtained the code: [https://github.com/labgem/ASMC/releases](https://github.com/labgem/ASMC/releases)

### Python dependencies

You can install the python dependencies with `pip`, `conda` or `mamba` with theses command and the files given in the releases:

**pip**
```
pip install -r requirements.txt
```

**conda**
```
conda env create -n env_name -f env.yml
```

**mamba**
```
mamba env create -n env_name -f env.yml
```

Installation via conda and mamba includes the modeller installation, but you still need to request the licence key.

### External dependencies

To install and configure the external dependencies, please follow the links in the requirements sections and the instructions given by their authors.

### Configuration

Before you can run the pipeline, you need to add a file named `config_asmc.yml` in `ASMC/ressources`. This file should contain the path to the `ASMC/ressources/AA_distances.tsv` and the path (or alias or binary name if it's in your PATH) to the p2rank and USalign executables, e.g:

```
distances: "<path>/ASMC/ressources/AA_distances.tsv"
usalign: "<path>/USalign"
p2rank: "<path>/prank"
```

## Quick Usage

```
python ASMC/run_asmc.py --threads 6 -r reference_file -s multi_fasta_file
```

the reference_file should contains the path to the reference(s) structure(s), e.g:
```
<path>/RefA.pdb
<path>/RefB.pdb
```

## Recommanded Usage

It is advisable to define the positions of interest manually, based on the literature or your own expertise.

It is possible to run only the pocket detection and then refine the selection yourself:

```
python ASMC/run_asmc.py --threads 6 -r reference_file -s multi_fasta_file --end pocket
```

Then you should modified (or created) the *pocket.csv* file returned.

Finally you can run the pipeline with this command :

```
python ASMC/run_asmc.py --threads 6 -r reference_file -p pocket.csv -s multi_fasta_file
```

For more details, see the [wiki](https://github.com/labgem/ASMC/wiki)