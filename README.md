PCA+CCA
=======

This repository contains the analysis scripts necessary for reproducing the results in "Expression
 reflects population structure." It contains three main files.

- Snakefile: A snakemake file that can be used to reproduce the entire analysis from start to
 finish, incuding downloading the data from public locations.
- pca_cca/util.py: A utility module containing the methods used in the analysis - implementations
 of coupled PCA and CCA as well as the leave-one-out projections, regression comparison, and
 permutation tests.
- make_config.py: A utility for constructing the Snakemake config file, useful for adjusting
 the number of components used, and selecting the methods and populations to use in the analysis

Running `python make_config.py` will produce the Snakemake config file equivalent to the
standard analysis, the only change being the number of permutations is reduced from 10M used
 in the main analysis to 10K because this is by far the most time consuming analysis step.
 To replicate the analysis in the paper, simply clone and `cd` into the repository, then type

```
$: python make_config.py
$: Snakemake results
```

Other targets include `projection_figures` which will make the PCA+CCA projection and
leave-one-out cross-validation plots, `preprocess_data` which will do all steps up
to and including the quantification of corrected transcript levels and genotype PCs, and `get_data`
 which will simply download the necessary data from public locations.

# Dependencies

This was developed using Python 3.6.2 with snakemake 3.13.3, numpy 1.13.3, pandas 0.20.3
 and sklearn 0.19.1.
The simplest way to meet these requirements is to use anaconda. For example:

```
:$ conda create -n PCA_CCA python=3.6.2
:$ conda install -n PCA_CCA -c bioconda snakemake
:$ source activate PCA_CCA
```

# Use as a library

This repository is currently indended for replication of the results in the manuscript. Advanced
users may be interested in using the utility methods in their own analysis, which is possible
by importing the utility module provided. At this time, the module is not documented for use
as a library in other analyses and such use is not officially supported. We intend to provide
a simple-to-use, well-documented implementation of our methods for use by others in their own
analyses at a later date.