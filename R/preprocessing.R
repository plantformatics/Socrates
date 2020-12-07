###################################################################################################
###################################################################################################
###################################################################################################
#' loadSparseData
#'
#' This function loads the sparse matrix and meta data and returns a list - counts & meta
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#'
#' @param input Path to sparse cell x peak text file in triplet format. Input data can be gzipped.
#' Required.
#' @param meta Path to meta data for cells in sparse cell x peak data set in tsv format. Headers
#' and barcode/cellID rownames are required. Meta data can be gzipped. Providing meta data is
#' highly recommended, but not required.
#' @param verbose logical. Defaults to FALSE.
#'
#' @rdname loadSparseData
#' @export
#'
loadSparseData <- function(input=NULL,
                           meta=NULL,
                           verbose=F){

    # check inputs
    if(is.null(input)){
        stop(" - Error: input is required (path/to/triplet/format/sparse/cell/by/peak/matrix.sparse)")
    }

    # check file types
    .filetype <- function(path){
        f = file(path)
        ext = summary(f)$class
        close.connection(f)
        ext
    }
    input.type <- .filetype(input)
    meta.type <- .filetype(meta)

    # load sparse
    if(input.type == "gzfile"){

        # verbose
        if(verbose){message(" - loading gzipped sparse matrix ... ")}
        a <- read.table(gzfile(input))
        a$V1 <- factor(a$V1)
        a$V2 <- factor(a$V2)

    }else{

        # verbose
        if(verbose){message(" - loading sparse matrix ... ")}
        a <- read.table(input)
        a$V1 <- factor(a$V1)
        a$V2 <- factor(a$V2)

    }

    # load meta
    if(!is.null(meta)){
        if(meta.type == "gzfile"){

            # verbose
            if(verbose){message(" - loading gzipped meta data ... ")}
            b <- read.table(gzfile(meta))

        }else{

            # verbose
            if(verbose){message(" - loading meta data ... ")}
            b <- read.table(meta)

        }
    }

    # convert to sparseMatrix format
    a <- Matrix::sparseMatrix(i=as.numeric(a$V1),
                              j=as.numeric(a$V2),
                              x=as.numeric(a$V3),
                              dimnames=list(levels(a$V1),levels(a$V2)))

    # order barcodes
    if(!is.null(meta)){
        both <- intersect(rownames(b), colnames(a))
        a <- a[,both]
        b <- b[both,]
    }else{
        b <- data.frame(cellIDs=colnames(a), row.names=colnames(a))
    }

    # make sure nSites is calculated
    b$nSites   <- Matrix::colSums(a)
    b$log10nSites <- log10(b$nSites)

    # return
    return(list(counts=a,meta=b))
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Filter sparse matrix
#'
#' This function filters the cell x peak matrix based on pre-clustering results, a minimum peak
#' usage frequency, and number of accessible peaks per cell
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#'
#' @param obj list (counts and meta) output from loadSparseData
#' @param preclusterID character column name in meta data for LSI or other clustering results to
#' condition ACR peak frequency. Sets the minimum frequency of accessibility across cells within a
#' single cluster, where sites that pass this threshold in at least one cluster are retained.
#' Defaults to NULL.
#' @param min.p numeric, minimum frequency of accessible cells per cluster (see preclusterID).
#' Sites above this threshold in at least one cluster are retained. Acceptable values range from 0
#' to 1. Defaults to 0.01.
#' @param min.t filters ACRs below the specified frequency of accessible cells across all cells.
#' Acceptable values range from 0 to 1.
#' @param min.c minimum number of accessible ACRs per cell. Defaults to 100.
#' @param max.t remove top X% of ACRs by frequency across cells (remove constitutive accessible
#' regions). Defaults to 0.005
#' @param min.cells.precluster minimum number cells from precluster (preclusterID) to include
#' cluster for thresholds set by min.p. Defaults to 10
#' @param verbose Defaults to FALSE
#'
#' @rdname cleanData
#' @export
#'
cleanData <- function(obj,
                      preclusterID=NULL,
                      min.p=0.01,
                      min.t=0.001,
                      min.c=100,
                      max.t=0.005,
                      min.cells.precluster=10,
                      verbose=F){

    # check object
    if(is.null(obj$counts)){
        stop(" - 'counts' slot missing ... terminating")
    }
    if(is.null(obj$meta)){
        stop(" - 'meta' slot missing ... terminating")
    }

    # verbose
    if(verbose){message("   * Input: cells = ", ncol(obj$counts), " | peaks = ", nrow(obj$counts))}

    # set vars
    x <- obj$counts
    y <- obj$meta

    # min peak count from dist
    y <- y[colnames(x),]

    # split cells by cluster, ie LSI_clusters
    if(!is.null(preclusterID)){
        num.clusts <- table(y[,preclusterID])
        num.clusts <- num.clusts[num.clusts > min.cells.precluster]
        clusts <- names(num.clusts)
        keep <- lapply(clusts, function(i){
            cells <- rownames(subset(y, y[,preclusterID]==i))
            sub.cells <- x[,colnames(x) %in% cells]
            if(length(cells)>1){
                peak.props <- Matrix::rowMeans(sub.cells)
                return(names(peak.props)[which(peak.props >= min.p)])
            }else{
                return(names(sub.cells)[which(sub.cells > 0)])
            }

        })
        keep <- do.call(c, keep)
        keep <- unname(keep, force=T)
        keep <- unique(keep)
        x <- x[rownames(x) %in% keep,]
    }

    # filter peaks by over all frequency
    x <- x[Matrix::rowSums(x)>(ncol(x)*min.t),]
    x <- x[Matrix::rowSums(x)<(quantile(Matrix::rowSums(x), c(1-max.t))),]

    # final clean
    if(min.c < 1){
        x <- x[,Matrix::colSums(x) > quantile(Matrix::colSums(x), min.c)]
    }else{
        x <- x[,Matrix::colSums(x)>min.c]
    }

    # last
    x <- x[Matrix::rowSums(x)>0,]
    x <- x[,Matrix::colSums(x)>0]

    # verbose
    if(verbose){message("   * Filtered: cells = ", ncol(x), " | peaks = ", nrow(x))}

    # return
    obj$counts <- x
    obj$meta <- obj$meta[colnames(obj$counts),]
    return(obj)

}


###################################################################################################
###################################################################################################
###################################################################################################
#' standardize or center residuals
#'
#' iterate over residual matrix, scaling or centering values across cells for each ACR.
#'
#' @param obj list object containing slot 'residuals' output from logisticModel or regModel functions.
#' @param center logical, whether or not to center residuals. Has no effect if scale is TRUE.
#' Defaults to TRUE.
#' @param scale logical, whether or not to center and standardize residuals. Defaults to TRUE.
#' @param bins numeric, number of bins to split ACRs into. Defaults to 1024.
#' @param verbose logical. Defaults to FALSE.
#'
#' @rdname standardizeResiduals
#' @export
#'
standardizeResiduals <- function(obj,
                                 center=T,
                                 scale=F,
                                 bins=1024,
                                 verbose=F){

    # verbose
    if(verbose){message(" - standardizing deviance scores ...")}

    # set up bins
    bin_ind <- ceiling(x = 1:nrow(obj$residuals) / bins)
    max_bin <- max(bin_ind)
    ids <- rownames(obj$residuals)

    # run in bins
    dev <- lapply(seq(1:max_bin), function(x){
        peaks_bin <- rownames(obj$residuals)[bin_ind == x]
        if(center==T & scale==F){
            t(apply(obj$residuals[peaks_bin,], 1, function(z){z-mean(z, na.rm=T)}))
        }else{
            t(apply(obj$residuals[peaks_bin,], 1, function(z){(z-mean(z, na.rm=T))/sd(z, na.rm=T)}))
        }
    })
    rm(obj$residuals)

    # merge and return
    dev <- do.call(rbind, dev)
    dev <- dev[ids,]

    # return scaled
    obj$residuals <- dev
    return(obj)
}
