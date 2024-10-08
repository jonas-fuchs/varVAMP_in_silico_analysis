## [varVAMP](https://github.com/jonas-fuchs/varVAMP) in silico analysis

This repository contains the complete code and input files for the in silico analysis of [varVAMP schemes](https://github.com/jonas-fuchs/ViralPrimerSchemes) and the code and data to reproduce the phylogenetic tree of HEV in the publication:

**"varVAMP: automated pan-specific primer design for tiled full genome sequencing and qPCR of highly diverse viral pathogens."**

### Requirements

For the phylogenetic tree, you will need ``R >= 4.2`` and for the remaining analysis ``python3 >= 3.9`` installed.

### Run the code

First clone this repository:

```bash
git clone github.com/jonas-fuchs/varVAMP_in_silico_analysis
cd varVAMP_in_silico_analysis
```

To reproduce the figures run:

```bash
pip install -r requirements.txt
python3 in_silico_eval.py
```

To reproduce the phylogenetic tree visualization run (might take a few mins):

```bash
Rscript tree_vis.R -t tree/HepE.aln.treefile -c tree/clus.tabular -o tree.pdf
```

To reproduce the input preparation for olivar, run:

```bash
python3 snp-sites_to_olivar.py
```


