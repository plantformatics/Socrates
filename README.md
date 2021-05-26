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

All users need to begin is a counts matrix (binarized) in triplet format, or a BED file (Columns 1-3: locations of Tn5 insertions, Column 4: barcode, and Column 5 strand) with genome annotation data (gff3/gtf and fai/chromosome.sizes). 

---

If you use `Socates` in your own study, please consider citing the following article:

> *Alexandre P. Marand, Zongliang Chen, Andrea Gallavotti, Robert J. Schmitz*. (2021). 
> [**A cis-regulatory atlas in maize at single-cell resolution**](https://doi.org/10.1016/j.cell.2021.04.014). 
> [**Cell**](https://doi.org/10.1016/j.cell.2021.04.014),  doi:10.1016/j.cell.2021.04.014

---

Current release: 03/16/21 BETA v0.0.9

## Installation

`Socates` requires R v4.0.0 or greater. 

```
# download the devtools package if not currently installed
install.packages("devtools")
library(devtools)

# install
devtools::install_github("plantformatics/Socrates", ref="main")
```

---

## Tutorials
[1. Loading data, quality control, identifying cells, and creating a Socrates object](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Load_and_QC_data.html)

[2. Loading pre-processed data, dimensionality reduction, and clustering](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Socrates_simple_clustering_tutorial.html)

[3. Comparison of normalization methods](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Socrates_clustering_comparison_tutorial_html.html)

[4. Identify co-accessible ACRs](https://htmlpreview.github.io/?https://github.com/plantformatics/Socrates/blob/main/vignettes/Socrates_coACRs_tutorial.html)

...  more coming soon
