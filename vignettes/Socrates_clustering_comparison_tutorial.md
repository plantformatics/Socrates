RRscca Tutorial: Comparing Normalization Methods
================
Alexandre P. Marand
10/27/2020

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
obj <- loadSparseData(input=input, meta=meta, verbose=T)
```

    ##  - loading gzipped sparse matrix ...

    ##  - loading meta data ...

------------------------------------------------------------------------

## Filter peaks and cells

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

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 13 x values <= 0 omitted from
    ## logarithmic plot

``` r
abline(v=500, col="red")
plot(density(site.freq), main="peak accessibility frequency", log="x")
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 7 x values <= 0 omitted from
    ## logarithmic plot

![](https://github.com/plantformatics/RRscca/blob/master/tutorial_figures/view_distributions-2.png)

It appears that most cells have a median over around 7,000 accessible peaks. The cells on the lower tail of the distribution may reflect broken nuclei, so we'll remove cells with less than 1,000 open peaks from the analysis. The distribution of average peak accessibilities doesnt show any clear (lower-tail) cut-offs, therefore, we will use the default thresholds (`min.t=0.001, max.t=0.005`). The arguments `min.t` and `max.t` set the minimum peak accessibility frequency and upper (99.5%) quantile cut-offs, respectively.

``` r
# filter matrix 
obj <- cleanData(obj, min.c=1000, min.t=0.001, max.t=0.005, verbose=T)
```

    ##    * Input: cells = 1500 | peaks = 89905

    ##    * Filtered: cells = 1498 | peaks = 88907

#### Filtering peaks enriched in pre-computed clusters

**NOTE** If users have generated crude clusters, such as *in silico* sorting described by [Cusanovich et al. 2018](https://linkinghub.elsevier.com/retrieve/pii/S0092-8674(18)30855-9), the parameters `min.p` and `preclusterID` allow users to set minimum peak accessibility frequencies for pre-specied groups. Below, we constrain peaks to be accessible in at least 5% of cells in at least one `crude_cluster`. See `?cleanData` for more details.

``` r
# simulate 10 random clusters, 
obj$meta$crude_clusters <- sample(factor(seq(1:10)), nrow(obj$meta), replace=T)
obj <- cleanData(obj, min.p=0.05, preclusterID="crude_clusters")
```

------------------------------------------------------------------------

## Normalization

With a filtered cell x peak matrix in hand, we can now calculate normalized accessibility scores (Pearson's residuals) across all cells and peaks using regularized quasibinomial logistic regression. Note that the function `regModel` can be parallelized by setting `nthreads` to a number greater than 1. Parallel implementations depend on the `doSNOW` library. In the example below, we set the number of threads to 4 to speed-up the analysis. In addition to regularized modeling, we provide additional functions to normalize chromatin accessibility profiles. Currently included in this release (supplement to `regModel`) are `tfidf`, `logisticModel`, and `regModel2`. In all cases, the output is saved in the 'residuals' slot of the output object (`obj$residuals`). We will run multiple normalization schemes to compare their performance. Each scheme is saved to a different slot.

#### Regularized logistic regression

The recommended normalization procedure.

``` r
# run regularized regression
obj <- regModel(obj, verbose=T, nthreads=4, slotName="NORM1")
```

    ##  - regularizing logistic model parameters ...

    ##    * formula: ~zlog10nSites

    ##    * estimating geometric mean ...

    ##    * density sampling on peak space ...

    ##    * fitting parameters to 5000 peaks ...

    ##    * finding outliers ...

    ##    * found 73 outliers ...

    ##    * regularizing all coefficients ...

    ##    * estimating residuals on the full data set ...

      |======================================================================| 100%

    ##    * residual range: -2.98779001451022 - 171.579569686042

#### TF-IDF

`tfidf` has the benefit of keeping the normalized data in a sparse format that allows users to conserve and reduce memory usage. We will compare TF-IDF normalization with quasibinomial later on in this tutorial. The TF-IDF function was adapted from Andrew Hill, the original implementation can be found [here](http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/)

``` r
# run TF-IDF normalization
obj <- tfidf(obj, slotName="NORM2")
```

#### Logistic regression without regularization

`logisticModel` does not regularize parameters by explicitly learning a model each peak independently. To run logistic regression with 4 threads without any regularization, run the following line of code. See `?logisticModel` for more details.

``` r
# run logisticModel
obj <- logisticModel(obj, verbose=T, nthreads=4, slotName="NORM3")
```

    ##  - running logistic regression ...

    ##    * formula: ~zlog10nSites

#### Regularized logistic regression with cell and peak sub-sampling.

As the number of cells can quickly become prohibitive (we have tested up to 60K cells across 160K peaks, requiring upwards of 50G memory), we extended the peak sub-sampling procedure for sampling cells. The function `regModel2` samples subsets of cells uniforming for factors specified in a column from the meta data, such as sampling 1,000 cells from different biological replicates or tissues. Sampling down to around 1,000 cells dramatically speeds up the computation with little effects on clustering. Here, we will sample 750 cells (half of the total 1,500) for building the models. See `?regModel2` for more details.

``` r
# run regularized model with cell sampling
obj$meta$library <- factor(1)
obj <- regModel2(obj, subcells=750, verbose=T, slotName="NORM4")
```

    ##  - regularizing logistic model parameters ...

    ##    * formula: ~zlog10nSites

    ##    * sampling cells ...

    ##    * sampled a total of 750 cells ...

    ##    * estimating geometric mean ...

    ##    * density sampling on peak space ...

    ##    * fitting parameters to 5000 peaks ...

    ##    * finding outliers ...

    ##    * found 69 outliers ...

    ##    * regularizing all coefficients ...

    ##    * estimating residuals on the full data set ...

    ##    * residual range: -2.91461115065793 - 141.815634021284

------------------------------------------------------------------------

## Reducing dimensions

### Singular Value Decomposition (SVD)

After normalizing peak x cell chromatin accessibility profiles using one of the aforementioned methods, we can reduce the dimensionality of residual matrix to remove noise and better model cell-cell relationships in a reduced space. Dimensionality reduction is implemented via Singular Value Decomposition (SVD) from the `irlba` package. We will reduce the dimensions of all the normalization methods by specifically selecting the slots from each approach one-by-one.

``` r
# reduce dimensionality of the residual matrix
obj <- reduceDims(obj, n.pcs=50, cor.max=0.7, verbose=T, residuals_slotName="NORM1", svd_slotName="SVD1")
```

    ##  - reduce dimensions with SVD ...

    ##  - standardizing reduced dimensions per cell ...

``` r
obj <- reduceDims(obj, n.pcs=50, cor.max=0.7, verbose=T, residuals_slotName="NORM2", svd_slotName="SVD2")
```

    ##  - reduce dimensions with SVD ... 
    ##  - standardizing reduced dimensions per cell ...

``` r
obj <- reduceDims(obj, n.pcs=50, cor.max=0.7, verbose=T, residuals_slotName="NORM3", svd_slotName="SVD3")
```

    ##  - reduce dimensions with SVD ... 
    ##  - standardizing reduced dimensions per cell ...

``` r
obj <- reduceDims(obj, n.pcs=50, cor.max=0.7, verbose=T, residuals_slotName="NORM4", svd_slotName="SVD4")
```

    ##  - reduce dimensions with SVD ... 
    ##  - standardizing reduced dimensions per cell ...

The above commands estimate the first 50 singular values, and removes singular values correlated with technical variation (read depth) above a Spearman Correlation Coefficient of 0.7.

### Projecting into a reduced dimensionality with UMAP

Similarity between cells is most easily visualized on two dimensions. Uniform Manifold Approximation Projection (UMAP) has gained popularity in single-cell approaches owing to its scalability and capacity for maintaining both global and local data structure, greatly aiding data interpretations. To generate UMAP embeddings, we can run the `projectUMAP` function which relies on `uwot::umap`:

``` r
# run projectUMAP
obj <- projectUMAP(obj, svd_slotName="SVD1", umap_slotName="UMAP1",verbose=T)
```

    ##  - non-linear dimensionality reduction with UMAP ...

    ## 12:43:44 UMAP embedding parameters a = 1.577 b = 0.8951

    ## 12:43:44 Read 1498 rows and found 50 numeric columns

    ## 12:43:44 Using Annoy for neighbor search, n_neighbors = 15

    ## 12:43:44 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 12:43:44 Writing NN index file to temp file /var/folders/_6/vtf34rcn6bd6xfg1_c3p34_m0000gp/T//RtmpJbEByC/filefe0f15a94d6d
    ## 12:43:44 Searching Annoy index using 2 threads, search_k = 1500
    ## 12:43:44 Annoy recall = 100%
    ## 12:43:45 Commencing smooth kNN distance calibration using 2 threads
    ## 12:43:45 Initializing from normalized Laplacian + noise
    ## 12:43:45 Commencing optimization for 500 epochs, with 32716 positive edges
    ## 12:43:47 Optimization finished

``` r
obj <- projectUMAP(obj, svd_slotName="SVD2", umap_slotName="UMAP2",verbose=T)
```

    ##  - non-linear dimensionality reduction with UMAP ...
    ## 12:43:47 UMAP embedding parameters a = 1.577 b = 0.8951
    ## 12:43:47 Read 1498 rows and found 49 numeric columns
    ## 12:43:47 Using Annoy for neighbor search, n_neighbors = 15
    ## 12:43:47 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 12:43:47 Writing NN index file to temp file /var/folders/_6/vtf34rcn6bd6xfg1_c3p34_m0000gp/T//RtmpJbEByC/filefe0f1e1a3a37
    ## 12:43:47 Searching Annoy index using 2 threads, search_k = 1500
    ## 12:43:47 Annoy recall = 100%
    ## 12:43:48 Commencing smooth kNN distance calibration using 2 threads
    ## 12:43:48 Initializing from normalized Laplacian + noise
    ## 12:43:48 Commencing optimization for 500 epochs, with 33866 positive edges
    ## 12:43:50 Optimization finished

``` r
obj <- projectUMAP(obj, svd_slotName="SVD3", umap_slotName="UMAP3",verbose=T)
```

    ##  - non-linear dimensionality reduction with UMAP ...
    ## 12:43:50 UMAP embedding parameters a = 1.577 b = 0.8951
    ## 12:43:50 Read 1498 rows and found 50 numeric columns
    ## 12:43:50 Using Annoy for neighbor search, n_neighbors = 15
    ## 12:43:50 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 12:43:50 Writing NN index file to temp file /var/folders/_6/vtf34rcn6bd6xfg1_c3p34_m0000gp/T//RtmpJbEByC/filefe0f4be40051
    ## 12:43:50 Searching Annoy index using 2 threads, search_k = 1500
    ## 12:43:50 Annoy recall = 100%
    ## 12:43:51 Commencing smooth kNN distance calibration using 2 threads
    ## 12:43:51 Initializing from normalized Laplacian + noise
    ## 12:43:51 Commencing optimization for 500 epochs, with 32330 positive edges
    ## 12:43:53 Optimization finished

``` r
obj <- projectUMAP(obj, svd_slotName="SVD4", umap_slotName="UMAP4",verbose=T)
```

    ##  - non-linear dimensionality reduction with UMAP ...
    ## 12:43:53 UMAP embedding parameters a = 1.577 b = 0.8951
    ## 12:43:53 Read 1498 rows and found 50 numeric columns
    ## 12:43:53 Using Annoy for neighbor search, n_neighbors = 15
    ## 12:43:53 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 12:43:53 Writing NN index file to temp file /var/folders/_6/vtf34rcn6bd6xfg1_c3p34_m0000gp/T//RtmpJbEByC/filefe0f65d0ecc3
    ## 12:43:53 Searching Annoy index using 2 threads, search_k = 1500
    ## 12:43:53 Annoy recall = 100%
    ## 12:43:54 Commencing smooth kNN distance calibration using 2 threads
    ## 12:43:54 Initializing from normalized Laplacian + noise
    ## 12:43:54 Commencing optimization for 500 epochs, with 32848 positive edges
    ## 12:43:56 Optimization finished

We can quickly visualize the reduced embedding. As you can see below, they are quite similar. In cases where speed and memory usage are central factors, it may be more advisable to use TF-IDF normalization in place of RRscca. One benefit of using model-based approaches is that the normalized values can be readily interpretted for down-stream analyses.

``` r
# plot UMAP results
layout(matrix(c(1:4), ncol=2, byrow=T))
par(mar=c(4,4,1,1))
plot(obj$UMAP1, pch=16, cex=0.2, main="regModel")
plot(obj$UMAP2, pch=16, cex=0.2, main="tfidf")
plot(obj$UMAP3, pch=16, cex=0.2, main="LogisticModel")
plot(obj$UMAP4, pch=16, cex=0.2, main="regModel2")
```

![](https://github.com/plantformatics/RRscca/blob/master/tutorial_figures/plotUMAPraw-2.png)

------------------------------------------------------------------------

## Graph-based clustering

Visualization of the UMAP embeddings suggests several groups of cells with distinct identities. We can cluster cells into groups using graph-based clustering by leaning on Louvain and Leiden clustering approaches provided by the popular `Seurat` package. Below is a wrapper for running graph-based clustering in the SVD space via Seurat. Cluster membership is appended to the meta data.frame under the column 'LouvainClusters' by default. For a list of tuneable parameters, run `?callClusters` in the R console.

``` r
# run clustering
obj <- callClusters(obj, svd_slotName="SVD1", umap_slotName="UMAP1", cluster_slotName="Clusters1", verbose=T)
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
    ## Number of nodes: 1476
    ## Number of edges: 80057
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 100 random starts: 0.9024
    ## Number of communities: 9
    ## Elapsed time: 0 seconds

    ##    * removing low quality clusters ...

    ##    * filtering per-cluster outliers (z-score filtDistClst2 = 5) ...

    ##    * total number of cells surviving subcluster filtering = 1467

``` r
obj <- callClusters(obj, svd_slotName="SVD2", umap_slotName="UMAP2", cluster_slotName="Clusters2", verbose=T)
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
    ## Number of nodes: 1478
    ## Number of edges: 98021
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 100 random starts: 0.8824
    ## Number of communities: 8
    ## Elapsed time: 0 seconds

    ##    * removing low quality clusters ...

    ##    * filtering per-cluster outliers (z-score filtDistClst2 = 5) ...

    ##    * total number of cells surviving subcluster filtering = 1460

``` r
obj <- callClusters(obj, svd_slotName="SVD3", umap_slotName="UMAP3", cluster_slotName="Clusters3", verbose=T)
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
    ## Number of nodes: 1486
    ## Number of edges: 81146
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 100 random starts: 0.8953
    ## Number of communities: 8
    ## Elapsed time: 0 seconds

    ##    * removing low quality clusters ...

    ##    * filtering per-cluster outliers (z-score filtDistClst2 = 5) ...

    ##    * total number of cells surviving subcluster filtering = 1477

``` r
obj <- callClusters(obj, svd_slotName="SVD4", umap_slotName="UMAP4", cluster_slotName="Clusters4", verbose=T)
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
    ## Number of nodes: 1481
    ## Number of edges: 82556
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 100 random starts: 0.9005
    ## Number of communities: 8
    ## Elapsed time: 0 seconds

    ##    * removing low quality clusters ...

    ##    * filtering per-cluster outliers (z-score filtDistClst2 = 5) ...

    ##    * total number of cells surviving subcluster filtering = 1475

## Plotting results

Having run the clustering algorithm, we can now visualize the different groupings on the UMAP embedding.

``` r
layout(matrix(c(1:4), ncol=2, byrow=T))
par(mar=c(4,4,1,1))
plotUMAP(obj, cluster_slotName="Clusters1", main="regModel")
plotUMAP(obj, cluster_slotName="Clusters2", main="tfidf")
plotUMAP(obj, cluster_slotName="Clusters3", main="logisticModel")
plotUMAP(obj, cluster_slotName="Clusters4", main="regModel2")
```

![](https://github.com/plantformatics/RRscca/blob/master/tutorial_figures/plotClusters-2.png)

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
