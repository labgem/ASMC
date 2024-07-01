# ASMC

ASMC combines (i) homology modeling of family members (MODELLER), (ii) ligand-binding pocket search (P2RANK), (iii) structural alignment of modeled active sites (USalign) and (iv) density-based spatial clustering of obtained alignments (DBSCAN) in a single command line. Clustering step can be carried out on either structural or sequence alignment.

<!-- ![ASMC workflow](images/ASMC_workflow.svg) -->
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

You can install the python dependencies with `pip`, `conda` or `mamba` with the following commands and the files given in the releases:

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

Installation via conda and mamba includes the MODELLER installation, but you still need to request the license key.

### External Software dependencies

- P2RANK - for ligand-binding pocket detection ([https://github.com/rdk/p2rank](https://github.com/rdk/p2rank))
- MODELLER - for homology modeling ([https://salilab.org/modeller/](https://salilab.org/modeller/))
- USalign - for structural alignment ([https://zhanggroup.org/US-align/](https://zhanggroup.org/US-align/))

 Please follow the links above and the instructions given by their authors.

## Configuration

In `ASMC/resources`, add a file exactly named `config_asmc.yml`. This file should contain 3 information:

- the path to the `ASMC/resources/AA_distances.tsv`
- the path of P2RANK executable
- the path of USalign executable

Example:
```yaml
distances: "/home/User/ASMC/resources/AA_distances.tsv"
usalign: "/home/User/bin/USALIGN/USalign"
p2rank: "/home/User/bin/p2rank_2.4.1/prank"
```

The keys should be identical to this example.

## Quick Usage

Run ASMC in a blind way (unknown active site) using a multi fasta file that should contain at least 100 sequences for clustering to be sufficiently relevant.

```
python ASMC/run_asmc.py --log run_asmc.log --threads 6 -r reference_file -s sequences.fasta
```

`reference_file` should contains the path to the reference(s) structure(s), e.g:
```
<path>/RefA.pdb
<path>/RefB.pdb
```

NB: For more details, see the [wiki](https://github.com/labgem/ASMC/wiki/Options-and-Usages)
