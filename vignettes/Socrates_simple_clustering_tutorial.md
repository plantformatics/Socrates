RRscca Tutorial: Simple Clustering Workflow
================
Alexandre P. Marand
10/28/2020

## Preprocessing and modeling chromatin accessibility with regularized quasibinomial logistic regression

scATAC-seq data is highly sparse and essentially binary for diploid cells, providing a significant challenge for downstream analyses. To mitigate these technical effects, we developed a model-based approach that we term Regularized quasibinomial logistic Regression for Single-Cell Chromatin Accessibility (RRscca). Inspired by innotation in scRNA-seq methods, namely the [SCTransform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) function implemented in [Seurat](https://satijalab.org/seurat/), RRscca removes technical variation confounding the accessibility signal from each cell and peak by fitting a generalized linear model between binary counts of chromatin accessibility and per-cell log10 read depths for each peak independently (`y~x`, where `y` represents a binary numeric vector of accessibility states across cells at a given peak, and `x` represents a vector of log10\[sum of accessible peaks\] per cell). `RRscca` explicitly avoids overfitting by sampling representative peaks with kernel regression, learning global quasibinomial model parameters (including a term for over-dispersion), and projecting the learned parameters to all peaks (regularization). The models are then used to extract Pearson's residuals that represent read depth-normalized chromatin accessibility profiles for each cell and peak. **Technical Note:** You will need access to 8G memory to run this example.

------------------------------------------------------------------------

## Load RRscca and raw data

Here, we demonstrate how to load the following inputs to create an `RRscca` object:

1.  **binary sparse peak x cell matrix** (a gzipped text file in triplet tsv format).
2.  **meta-data** saved as a tsv document containing various per-cell metrics.

Meta-data is not necessarily required, but useful for evaluating technical effects in clustering, as well as other downstream steps. An example of input data formats for 1,500 PBMC cells and how to construct an `RRscca` object from scratch is shown below.

``` r
# load library
suppressWarnings(library("RRscca"))

# specify paths to raw data in RRscca package
input <- system.file("extdata", "pbmc_atac_10x.1.5K_cells.sparse.gz", package="RRscca")
meta <- system.file("extdata", "pbmc_atac_10x.1.5K_cells.metadata.txt", package="RRscca")

# load raw data for viewing
input.format <- read.table(input)
meta.format <- read.table(meta)

# view
head(input.format)
```

    ##                     V1                 V2 V3
    ## 1   chr1_713551_714783 CACTGAACAGTTGCAT-1  1
    ## 2   chr1_762249_763345 CACTGAACAGTTGCAT-1  1
    ## 3   chr1_839587_841090 CACTGAACAGTTGCAT-1  1
    ## 4   chr1_911144_911956 CACTGAACAGTTGCAT-1  1
    ## 5 chr1_1003706_1006207 CACTGAACAGTTGCAT-1  1
    ## 6 chr1_1057212_1058134 CACTGAACAGTTGCAT-1  1

``` r
head(meta.format)
```

    ##                                cellID total duplicate chimeric unmapped lowmapq
    ## AGGCCCAAGTGAAGGA-1 AGGCCCAAGTGAAGGA-1 15947      7372      108       91     900
    ## GCTTGCTAGGAATGGA-1 GCTTGCTAGGAATGGA-1 31714     15398      353      162    1582
    ## ATGTCTTGTTTGATCG-1 ATGTCTTGTTTGATCG-1 23718     10135      228      137    1302
    ## GGTAGGATCTTATCAC-1 GGTAGGATCTTATCAC-1 24968     12068      306      139    1449
    ## CTACAGATCTGAACGT-1 CTACAGATCTGAACGT-1 18681      9362      219      116    1037
    ## TCTATTGAGAGACTCG-1 TCTATTGAGAGACTCG-1 28407     11432      368      171    1602
    ##                    mitochondrial passed_filters    cell_id is__cell_barcode
    ## AGGCCCAAGTGAAGGA-1             0           7476 _cell_1520                1
    ## GCTTGCTAGGAATGGA-1             0          14219 _cell_5602                1
    ## ATGTCTTGTTTGATCG-1             5          11911 _cell_1987                1
    ## GGTAGGATCTTATCAC-1             0          11006 _cell_6099                1
    ## CTACAGATCTGAACGT-1            68           7879 _cell_3691                1
    ## TCTATTGAGAGACTCG-1             0          14834 _cell_8134                1
    ##                    TSS_fragments DNase_sensitive_region_fragments
    ## AGGCCCAAGTGAAGGA-1          3744                             6524
    ## GCTTGCTAGGAATGGA-1          6354                            12101
    ## ATGTCTTGTTTGATCG-1          4662                             9791
    ## GGTAGGATCTTATCAC-1          7051                             9290
    ## CTACAGATCTGAACGT-1          4443                             6604
    ## TCTATTGAGAGACTCG-1          6619                            12212
    ##                    enhancer_region_fragments promoter_region_fragments
    ## AGGCCCAAGTGAAGGA-1                       976                      3839
    ## GCTTGCTAGGAATGGA-1                      2066                      6511
    ## ATGTCTTGTTTGATCG-1                      1886                      4688
    ## GGTAGGATCTTATCAC-1                       764                      7376
    ## CTACAGATCTGAACGT-1                       705                      4589
    ## TCTATTGAGAGACTCG-1                      2003                      6728
    ##                    on_target_fragments blacklist_region_fragments
    ## AGGCCCAAGTGAAGGA-1                6821                          4
    ## GCTTGCTAGGAATGGA-1               12753                         11
    ## ATGTCTTGTTTGATCG-1               10371                          6
    ## GGTAGGATCTTATCAC-1               10122                          4
    ## CTACAGATCTGAACGT-1                7039                          7
    ## TCTATTGAGAGACTCG-1               13051                          9
    ##                    peak_region_fragments peak_region_cutsites
    ## AGGCCCAAGTGAAGGA-1                  5401                10513
    ## GCTTGCTAGGAATGGA-1                 10004                19385
    ## ATGTCTTGTTTGATCG-1                  7466                14300
    ## GGTAGGATCTTATCAC-1                  8811                17346
    ## CTACAGATCTGAACGT-1                  5675                11134
    ## TCTATTGAGAGACTCG-1                  9959                19264

``` r
# load data into RRscca object
rrscca.object <- loadSparseData(input=input, meta=meta, verbose=T)
```

    ##  - loading gzipped sparse matrix ...

    ##  - loading meta data ...

------------------------------------------------------------------------

## Filter peaks and cells

For the remainder of this tutorial, we will work from a 5K PBMC data set publically available from the 10X Genomics [website](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_pbmc_10k_nextgem). A precompiled `RRscca` binary object containing this data is automatically available after loading the `RRscca` package.

``` r
str(obj)
```

    ## List of 2
    ##  $ counts:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   .. ..@ i       : int [1:35591942] 0 14 17 20 27 30 32 51 70 79 ...
    ##   .. ..@ p       : int [1:5001] 0 11371 23721 33022 37976 44327 52028 54701 63390 70036 ...
    ##   .. ..@ Dim     : int [1:2] 90327 5000
    ##   .. ..@ Dimnames:List of 2
    ##   .. .. ..$ : chr [1:90327] "chr1_10002383_10004137" "chr1_10010207_10011578" "chr1_100151234_100151687" "chr1_100165510_100166120" ...
    ##   .. .. ..$ : chr [1:5000] "CCTGCTAAGTGAAACT-1" "ATTACTCAGCTACGCC-1" "TATCGAGGTCACAGTT-1" "CAGTATGAGTTCTCCC-1" ...
    ##   .. ..@ x       : num [1:35591942] 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. ..@ factors : list()
    ##  $ meta  :'data.frame':  5000 obs. of  20 variables:
    ##   ..$ cellID                          : Factor w/ 9668 levels "AAACGAAAGAGCTGTG-1",..: 3216 2033 7499 2615 7518 9359 1986 7988 7508 7115 ...
    ##   ..$ total                           : int [1:5000] 54346 65983 49121 16175 22659 34150 9523 33538 33572 41436 ...
    ##   ..$ duplicate                       : int [1:5000] 27690 31177 25083 6531 9363 16138 4271 16523 13221 19528 ...
    ##   ..$ chimeric                        : int [1:5000] 666 919 649 175 176 486 97 336 301 484 ...
    ##   ..$ unmapped                        : int [1:5000] 250 329 296 79 114 206 67 173 148 219 ...
    ##   ..$ lowmapq                         : int [1:5000] 2905 3831 3016 899 1248 1876 593 1554 2031 2422 ...
    ##   ..$ mitochondrial                   : int [1:5000] 188 21 34 0 1 9 74 58 584 4 ...
    ##   ..$ passed_filters                  : int [1:5000] 22647 29706 20043 8491 11757 15435 4421 14894 17287 18779 ...
    ##   ..$ cell_id                         : Factor w/ 9668 levels "_cell_0","_cell_1",..: 2465 1151 7222 1797 7244 9289 1097 7765 7233 6797 ...
    ##   ..$ is__cell_barcode                : int [1:5000] 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ TSS_fragments                   : int [1:5000] 11850 16354 12167 3581 4650 9028 2236 6533 4863 10587 ...
    ##   ..$ DNase_sensitive_region_fragments: int [1:5000] 19538 24968 16643 7147 9718 12837 3750 13067 12062 16229 ...
    ##   ..$ enhancer_region_fragments       : int [1:5000] 2680 2843 1625 1273 1851 1257 514 2164 2688 1747 ...
    ##   ..$ promoter_region_fragments       : int [1:5000] 12255 16868 12599 3681 4680 9335 2316 6893 4857 10983 ...
    ##   ..$ on_target_fragments             : int [1:5000] 20661 26662 18148 7520 10312 13941 3929 13703 12918 17244 ...
    ##   ..$ blacklist_region_fragments      : int [1:5000] 17 24 9 3 5 9 1 3 43 10 ...
    ##   ..$ peak_region_fragments           : int [1:5000] 16741 21493 15077 5739 7492 11453 3016 11263 8002 13996 ...
    ##   ..$ peak_region_cutsites            : int [1:5000] 32599 42140 29604 11082 14365 22495 5897 21877 15419 27306 ...
    ##   ..$ nSites                          : num [1:5000] 11371 12350 9301 4954 6351 ...
    ##   ..$ log10nSites                     : num [1:5000] 4.06 4.09 3.97 3.69 3.8 ...

To reduce the effects of outlier peaks and cells on clustering, it is often helpful to remove cells with few accessible peaks, and peaks with extreme accessibility profiles (i.e. peaks that are accessible in all, or very few cells). Specifically, the cell x peak matrix can be filtered by adjusting heurstic frequency thresholds after visual inspection of cell and peak accessibility distributions. Let's investigate these distributions further to determine reasonable thresholds for this particular data set.

``` r
# estimate log10 number of accessible regions per cell
cell.counts <- Matrix::colSums(obj$counts)

# estimate peak accessibility frequency across cells
site.freq <- Matrix::rowMeans(obj$counts)

# plot distributions
layout(matrix(c(1:2), ncol=2))
par(mar=c(3,3,1,1))
plot(density(cell.counts), main="log10 cell counts", log="x")
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 10 x values <= 0 omitted from
    ## logarithmic plot

``` r
abline(v=1000, col="red")
plot(density(site.freq), main="peak accessibility frequency", log="x")
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 7 x values <= 0 omitted from
    ## logarithmic plot

![](https://github.com/plantformatics/RRscca/blob/master/tutorial_figures/view_distributions-1.png)

It appears that most cells have a median over around 7,000 accessible peaks. The cells on the lower tail of the distribution may reflect broken nuclei, so we'll remove cells with less than 1,000 open peaks from the analysis. The distribution of average peak accessibilities doesnt show any clear (lower-tail) cut-offs, therefore, we will use the default thresholds (`min.t=0.001, max.t=0.005`). The arguments `min.t` and `max.t` set the minimum peak accessibility frequency and upper (99.5%) quantile cut-offs, respectively.

``` r
# filter matrix 
obj <- cleanData(obj, min.c=1000, min.t=0.001, max.t=0.005, verbose=T)
```

    ##    * Input: cells = 5000 | peaks = 90327

    ##    * Filtered: cells = 4984 | peaks = 88953

#### Filtering peaks enriched in pre-computed clusters

**NOTE** If users have generated crude clusters, such as *in silico* sorting described by [Cusanovich et al. 2018](https://linkinghub.elsevier.com/retrieve/pii/S0092-8674(18)30855-9), the parameters `min.p` and `preclusterID` allow users to set minimum peak accessibility frequencies for pre-specied groups. Below, we constrain peaks to be accessible in at least 5% of cells in at least one `crude_cluster`. See `?cleanData` for more details.

``` r
# simulate 10 random clusters, 
obj$meta$crude_clusters <- sample(factor(seq(1:10)), nrow(obj$meta), replace=T)
obj <- cleanData(obj, min.p=0.05, preclusterID="crude_clusters")
```

------------------------------------------------------------------------

## Normalization

With a filtered cell x peak matrix in hand, we can now calculate normalized accessibility scores (Pearson's residuals) across all cells and peaks using regularized quasibinomial logistic regression. Note that the function `regModel` can be parallelized by setting `nthreads` to a number greater than 1. Parallel implementations depend on the `doSNOW` library. In the example below, we set the number of threads to 4 to speed-up the analysis. We will run multiple normalization schemes to compare their performance.

``` r
# run regularized regression
obj <- regModel(obj, verbose=T, nthreads=4)
```

    ##  - regularizing logistic model parameters ...

    ##    * formula: ~zlog10nSites

    ##    * estimating geometric mean ...

    ##    * density sampling on peak space ...

    ##    * fitting parameters to 5000 peaks ...

    ##    * finding outliers ...

    ##    * found 84 outliers ...

    ##    * regularizing all coefficients ...

    ##    * estimating residuals on the full data set ...

      |======================================================================| 100%

    ##    * residual range: -3.00377470612868 - 135.858192652319

------------------------------------------------------------------------

## Reducing dimensions

### Singular Value Decomposition (SVD)

After normalizing peak x cell chromatin accessibility profiles using one of the aforementioned methods, we can reduce the dimensionality of residual matrix to remove noise and better model cell-cell relationships in a reduced space. Dimensionality reduction is implemented via Singular Value Decomposition (SVD) from the `irlba` package. We will reduce the dimensions of all the normalization methods by specifically selecting the slots from each approach one-by-one.

``` r
# reduce dimensionality of the residual matrix
obj <- reduceDims(obj, n.pcs=50, cor.max=0.7, verbose=T)
```

    ##  - reduce dimensions with SVD ...

    ##  - standardizing reduced dimensions per cell ...

The above commands estimate the first 50 singular values, and removes singular values correlated with technical variation (read depth) above a Spearman Correlation Coefficient of 0.7.

### Projecting into a reduced dimensionality with UMAP

Similarity between cells is most easily visualized on two dimensions. Uniform Manifold Approximation Projection (UMAP) has gained popularity in single-cell approaches owing to its scalability and capacity for maintaining both global and local data structure, greatly aiding data interpretations. To generate UMAP embeddings, we can run the `projectUMAP` function which relies on `uwot::umap`:

``` r
# run projectUMAP
obj <- projectUMAP(obj, verbose=T)
```

    ##  - non-linear dimensionality reduction with UMAP ...

    ## 11:59:30 UMAP embedding parameters a = 1.577 b = 0.8951

    ## 11:59:30 Read 4984 rows and found 50 numeric columns

    ## 11:59:30 Using Annoy for neighbor search, n_neighbors = 15

    ## 11:59:30 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 11:59:31 Writing NN index file to temp file /var/folders/_6/vtf34rcn6bd6xfg1_c3p34_m0000gp/T//RtmpNP8NNb/filefb8a657471d8
    ## 11:59:31 Searching Annoy index using 2 threads, search_k = 1500
    ## 11:59:31 Annoy recall = 100%
    ## 11:59:32 Commencing smooth kNN distance calibration using 2 threads
    ## 11:59:33 Initializing from normalized Laplacian + noise
    ## 11:59:33 Commencing optimization for 500 epochs, with 109596 positive edges
    ## 11:59:38 Optimization finished

We can quickly visualize the reduced embedding. As you can see below, they are quite similar. In cases where speed and memory usage are central factors, it may be more advisable to use TF-IDF normalization in place of RRscca. One benefit of using model-based approaches is that the normalized values can be readily interpretted for down-stream analyses.

``` r
# plot UMAP results
plot(obj$UMAP, pch=16, cex=0.2)
```

![](https://github.com/plantformatics/RRscca/blob/master/tutorial_figures/plotUMAPraw-1.png)

------------------------------------------------------------------------

## Graph-based clustering

Visualization of the UMAP embeddings suggests several groups of cells with distinct identities. We can cluster cells into groups using graph-based clustering by leaning on Louvain and Leiden clustering approaches provided by the popular `Seurat` package. Below is a wrapper for running graph-based clustering in the SVD space via Seurat. Cluster resolution can be adjusted from the default (0.4) using the `res` argument. Reasonable `res` values range from 0 to 2 (high values yield higher cluster numbers). Cluster membership is appended to the meta data.frame under the column 'LouvainClusters' by default. For a list of tuneable parameters, run `?callClusters` in the R console.

``` r
# run clustering
obj <- callClusters(obj, res=0.6, verbose=T)
```

    ##  - filtering outliers in UMAP manifold (z-score e.thresh = 3) ...

    ##  - creating seurat object for graph-based clustering ...

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## Warning: The following arguments are not used: reduction

    ## Warning: The following arguments are not used: reduction

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4926
    ## Number of edges: 326242
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 100 random starts: 0.8871
    ## Number of communities: 13
    ## Elapsed time: 6 seconds

    ## 1 singletons identified. 12 final clusters.

    ##    * removing low quality clusters ...

    ##    * filtering per-cluster outliers (z-score filtDistClst2 = 5) ...

    ##    * total number of cells surviving subcluster filtering = 4890

## Plotting results

Having run the clustering algorithm, we can now visualize the different groupings on the UMAP embedding.

``` r
plotUMAP(obj)
```

![](https://github.com/plantformatics/RRscca/blob/master/tutorial_figures/plotClusters-1.png)

## Saving results

The `RRscca` object is updated iteratively after each processing step. To save results for sharing or future exploration, you can use the following command to save a snapshot of current stage of analysis:

``` r
saveRDS(obj, file="RRscca_object.rds")
```

## Accessing results

Results are appended to the `RRscca` object after each function, as described above. Below describes the location of different data sets up to this point.

-   Raw counts
    -   `obj$counts`
    -   generated by `loadSparseData`
-   Raw meta-data
    -   `obj$meta`
    -   generated by `loadSparseData`
-   SVD/PCA
    -   `obj$PCA` (Default, can be customized)
    -   generated by `reduceDims`
-   Initial UMAP
    -   `obj$UMAP` (Default, can be customized)
    -   generated by `projectUMAP`
-   Cluster + meta-data (filtered cells)
    -   `obj$Clusters` (Default, can be customized)
    -   generated by `callClusters`

------------------------------------------------------------------------

### Session Information

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Catalina 10.15.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] RRscca_0.0.0.9000
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Seurat_3.2.2          Rtsne_0.15            colorspace_1.4-1     
    ##   [4] deldir_0.1-29         ellipsis_0.3.1        ggridges_0.5.2       
    ##   [7] spatstat.data_1.4-3   leiden_0.3.3          listenv_0.8.0        
    ##  [10] ggrepel_0.8.2         RSpectra_0.16-0       codetools_0.2-16     
    ##  [13] splines_3.6.3         knitr_1.30            itertools_0.1-3      
    ##  [16] polyclip_1.10-0       jsonlite_1.7.1        ica_1.0-2            
    ##  [19] cluster_2.1.0         png_0.1-7             uwot_0.1.8           
    ##  [22] shiny_1.5.0           sctransform_0.3.1     compiler_3.6.3       
    ##  [25] httr_1.4.2            Matrix_1.2-18         fastmap_1.0.1        
    ##  [28] lazyeval_0.2.2        later_1.1.0.1         htmltools_0.5.0      
    ##  [31] tools_3.6.3           rsvd_1.0.3            igraph_1.2.6         
    ##  [34] gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
    ##  [37] reshape2_1.4.4        dplyr_1.0.2           Rcpp_1.0.5           
    ##  [40] spatstat_1.64-1       vctrs_0.3.4           nlme_3.1-150         
    ##  [43] iterators_1.0.13      lmtest_0.9-38         xfun_0.18            
    ##  [46] stringr_1.4.0         globals_0.13.1        mime_0.9             
    ##  [49] miniUI_0.1.1.1        lifecycle_0.2.0       irlba_2.3.3          
    ##  [52] goftest_1.2-2         future_1.19.1         MASS_7.3-53          
    ##  [55] zoo_1.8-8             scales_1.1.1          doSNOW_1.0.19        
    ##  [58] promises_1.1.1        spatstat.utils_1.17-0 parallel_3.6.3       
    ##  [61] RColorBrewer_1.1-2    yaml_2.2.1            reticulate_1.18      
    ##  [64] pbapply_1.4-3         gridExtra_2.3         ggplot2_3.3.2        
    ##  [67] rpart_4.1-15          stringi_1.5.3         S4Vectors_0.24.4     
    ##  [70] foreach_1.5.1         BiocGenerics_0.32.0   BiocParallel_1.20.1  
    ##  [73] shape_1.4.5           rlang_0.4.8           pkgconfig_2.0.3      
    ##  [76] matrixStats_0.57.0    evaluate_0.14         lattice_0.20-41      
    ##  [79] ROCR_1.0-11           purrr_0.3.4           tensor_1.5           
    ##  [82] patchwork_1.0.1       htmlwidgets_1.5.2     cowplot_1.1.0        
    ##  [85] tidyselect_1.1.0      RcppAnnoy_0.0.16      plyr_1.8.6           
    ##  [88] magrittr_1.5          R6_2.4.1              IRanges_2.20.2       
    ##  [91] snow_0.4-3            generics_0.0.2        DelayedArray_0.12.3  
    ##  [94] pillar_1.4.6          mgcv_1.8-33           fitdistrplus_1.1-1   
    ##  [97] survival_3.2-7        abind_1.4-5           tibble_3.0.4         
    ## [100] future.apply_1.6.0    crayon_1.3.4          KernSmooth_2.23-17   
    ## [103] plotly_4.9.2.1        rmarkdown_2.5         viridis_0.5.1        
    ## [106] grid_3.6.3            data.table_1.13.2     FNN_1.1.3            
    ## [109] digest_0.6.27         xtable_1.8-4          tidyr_1.1.2          
    ## [112] httpuv_1.5.4          glmnet_4.0-2          stats4_3.6.3         
    ## [115] munsell_0.5.0         viridisLite_0.3.0
