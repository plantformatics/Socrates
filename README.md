# Socrates

![](https://github.com/plantformatics/Socrates/blob/master/tutorial_figures/Socrates_home_image.png)

`Socrates` is an R package for analyzing Single-Cell Assay for Transposase Accessible Chromatin Sequencing (scATAC-seq) data. `Socrates` takes as input **1.** barcode meta data and **2.** a cell x feature sparse matrix in triplet format (see example data in inst/extdata). The main contribution of `Socrates` compared to previously established methods is a regularized quasi-binomial logistic regression for single-cell chromatin accessibility for normalizing accessibility profiles across peaks and cells. 

Coming updates to `Socrates` will provides tool for several processing and analytical steps: 

1. Clustering
2. Batch effect removal 
3. Cell-type annotation 
4. Co-accessibility 
5. Motif analysis 
6. Gene accessibility 
7. Pseudotime 
8. scRNA-seq integration

and much more!

All users need to begin is a counts matrix (binarized) in triplet format. 

---

If you use `Socates` in your own study, please consider citing the following article:

> [**A cis-regulatory atlas in maize at single-cell resolution**](https://www.biorxiv.org/content/10.1101/2020.09.27.315499v3)
> *Alexandre P. Marand, Zongliang Chen, Andrea Gallavotti, Robert J. Schmitz*
> bioRxiv 2020.09.27.315499; doi: [https://doi.org/10.1101/2020.09.27.315499](https://www.biorxiv.org/content/10.1101/2020.09.27.315499v3)

---

Current release: 10/27/20 BETA v0.0.9

## Installation

`Socates` requires R v3.6.3 or greater. 

```
# download the devtools package if not currently installed
install.packages("devtools")
library(devtools)

# install
devtools::install_github("plantformatics/Socrates", ref="main")
```

---

## Tutorials
[1. Pre-processing, dimensionality reduction, and clustering](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Socrates_simple_clustering_tutorial.html)

[2. Comparison of normalization methods](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Socrates_clustering_comparison_tutorial.html)

[3. Identify co-accessible ACRs](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Socrates_coACRs_tutorial.html)

...  more coming soon
