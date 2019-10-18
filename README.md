# Sexually antagonistic variation and the evolution of dimorphic sexual systems 

## Overview

This is the GitHub repository for the development of a population genetics theory project about the evolution of dioecy and incipient sex chromosomes. The published article is now [available through the publisher](https://www.journals.uchicago.edu/doi/full/10.1086/702847). 

* Shortly after publication, the authors noticed a mistake in the analytic results pertaining to the invasion of sexually antagonistic alleles (hereafter, “SA alleles”) linked to a unisexual sterility allele. The correction is available [here](https://www.journals.uchicago.edu/doi/full/10.1086/705014).

You can also contact me directly if you would like a reprint. Here you can find all of the necessary code to reproduce the analyses and the pre-publication manuscript. The *Mathematica* `.nb` files where the models are developed are not tracked by git, but can be downloaded as an online appendix from the publisher website above (or contact me and I will be happy to provide if you don't have access to the appendices). Aside from the `.nb` files, all necessary code for creating the figures can be found in the `./R/functions-*.R` files. 


## Citing information

Olito, C. and Connallon, T. 2019. Sexually antagonistic variation and the evolution of dimorphic sexual systems. *American Naturalist* 193: 688--701. doi:10.1086/702847

* Olito, C. and Connallon, T. 2019. *Correction*. *American Naturalist* 194:741--742.

## Reproducing the manuscript

The easiest way to reproduce the manuscript is to simply clone the repo, run `createFigs.R`, and then compile the manuscript file `doc/dioecySA.tex` using whatever default LaTeX editor/engine you have. To reproduce the simulations etc., you will need to create the  run the R scripts `recSim-*-RunSims.R`


## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the dioecySA github [issues page](https://github.com/colin-olito/dioecySA/issues). If you would like to report a bug/issue, and do not have a github account (and don't want to get one), please send a brief email detailing the problem you encountered to colin.olito at monash dot edu.



