# The chemosynthetic microbiome of deep-sea hydrothermal vents across space and time
[![forthebadge](https://forthebadge.com/images/badges/cc-by-nd.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/powered-by-coffee.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)

[![deepseamicrobiologylab](https://img.shields.io/badge/BY-DeepseaMicrobiologyLab-blue)](https://marine.rutgers.edu/deep-seamicrobiology/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13899446.svg)](https://doi.org/10.5281/zenodo.13899446)
[![made-with-Markdown](https://img.shields.io/badge/Coded%20in-R-red.svg)](https://www.r-project.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)

Data and code relative to the Grosche et al. 2025 manuscript about the chemosynthetic microbiome of deepsea hydrothermal vents in the East Pacific Rise. Sequences are available through the NCBI SRA archive with accession numbers PRJNA1092253 for the amplicon data and through the MG-RAST archive with accession numbers MGM4773053.3 and MGM4773182.3 for the metatranscriptomic data.

The repository contains the following files:

- README.md (this readme file)
- 16S_rRNA folder with 2 subfolders: dada2 and phyloseq
  - grosche_et_al_dada2.r (the r code to reproduce the dada2 analysis)
  - grosche_et_al_phyloseq.r (the r code to reproduce the figures)
  - a series of .csv files to replicate the figures and analyses reported in the paper
- Metatranscriptomic folder
  - Grosche_et_al_metatranscriptomic.r (the r code to reproduce the figures)
  - metaxa_metatranscriptomic.r (the r code to reproduce the figures)
  - a series of .csv files to replicate the figures and analyses reported in the paper

Please cite as:
Ashley Grosche, Matteo Selci, Francesco Smedile, Donato Giovannelli, Sara Borin, Nadine Le Bris, and Costantino Vetriani. 2024. The chemosynthetic microbiome of deep-sea hydrothermal vents across space and time. -Submitted-
