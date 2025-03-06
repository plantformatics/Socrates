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
#' convertSparseData
#'
#' This function loads the sparse matrix and meta data and returns a list - counts & meta
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#'
#' @param obj Object output by generateMatrix.
#' @param verbose logical. Defaults to FALSE.
#'
#' @rdname convertSparseData
#' @export
#'
convertSparseData <- function(obj,
                              verbose=F){

    # specify data sets
    a <- obj$counts

    # make sure bins/cells are factors
    if(verbose){message(" - converting triplet format to sparseMatrix")}
    a$V1 <- factor(a$V1)
    a$V2 <- factor(a$V2)

    # convert to sparseMatrix format
    a <- Matrix::sparseMatrix(i=as.numeric(a$V1),
                              j=as.numeric(a$V2),
                              x=as.numeric(a$V3),
                              dimnames=list(levels(a$V1),levels(a$V2)))

    # order barcodes
    if(!is.null(obj$meta.v3)){
        meta.use <- obj$meta.v3
    }else{
        if(!is.null(obj$meta.v2)){
            meta.use <- obj$meta.v2
        }else{
            if(!is.null(obj$meta.v1)){
                meta.use <- obj$meta.v1
            }else{
                meta.use <- obj$meta
            }
        }
    }
    b <- meta.use

    # align barcodes
    both <- intersect(rownames(b), colnames(a))
    a <- a[,both]
    b <- b[both,]

    # make sure nSites is calculated
    b$nSites   <- Matrix::colSums(a)
    b$log10nSites <- log10(b$nSites)
    obj$counts <- a
    obj$meta <- b
  
    # return
    return(obj)
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
                      min.t=0.01,
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


estGeneActivity <- function(obj, FeatureName="gene", pRange=50000, dRange=100, con=10000, per.cell.scale=10000){

    ## Load Tn5 file with GRanges
    Tn5_Grange <-  GRanges(seqnames = obj$bed$V1,
                           ranges = IRanges(start = obj$bed$V2,
                                            end = obj$bed$V3,
                                            names = obj$bed$V4))

    ## Select Features
    if(FeatureName == "gene"){
        Ann <- genes(obj$gff)
    }else if(FeatureName == "transcript"){
        Ann <- transcripts(obj$gff)
    }else{
        message("ERROR: Feature name should be 'gene' or 'transcript")
    }

    ## Find overlap in genes
    message(" - building sparse gene body matrix ...")
    hits_Within <- suppressWarnings(findOverlaps(Tn5_Grange,  Ann, minoverlap=1,
                                type=c("within"), select="all", ignore.strand = TRUE))
    Intersect <- paste(names(Ann)[hits_Within@to],
                       names(Tn5_Grange)[hits_Within@from],sep="/_Com_/")
    Intersect <- table(Intersect)
    Intersect <- data.frame(geneID = as.factor(as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                    split="/_Com_/"), "[", 1))),
                           cellID = as.factor(as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                 split="/_Com_/"), "[", 2))),
                           activity = as.numeric(as.character(Intersect)))
    obj$gba <- sparseMatrix(i=as.numeric(Intersect$geneID),
                            j=as.numeric(Intersect$cellID),
                            x=Intersect$activity,
                            dimnames=list(levels(Intersect$geneID),
                                          levels(Intersect$cellID)))
    rm(Intersect)
    rm(hits_Within)
    all.genes <- rownames(obj$gba)

    # find overlaps in peaks
    message(" - building sparse ACR matrix ...")
    acrs <- GRanges(seqnames=obj$acr$V1,
                    ranges=IRanges(start=obj$acr$V2,
                                   end=obj$acr$V3,
                                   names=paste(obj$acr$V1,obj$acr$V2,obj$acr$V3, sep="_")))
    peak_hits <- suppressWarnings(findOverlaps(Tn5_Grange, acrs,
                                               minoverlap=1,
                                               type=c("within"),
                                               select="all",
                                               ignore.strand=T))

    peak_int <- paste(names(acrs)[peak_hits@to],
                       names(Tn5_Grange)[peak_hits@from],sep="/_Com_/")
    peak_int <- table(peak_int)
    peak_int <- data.frame(acrID = as.factor(as.character(lapply(strsplit(as.character(rownames(peak_int)),
                                                                 split="/_Com_/"), "[", 1))),
                           cellID = as.factor(as.character(lapply(strsplit(as.character(rownames(peak_int)),
                                                                 split="/_Com_/"), "[", 2))),
                           activity = as.numeric(as.character(peak_int)))
    obj$pa <- sparseMatrix(i=as.numeric(peak_int$acrID),
                           j=as.numeric(peak_int$cellID),
                           x=peak_int$activity,
                           dimnames=list(levels(peak_int$acrID),
                                         levels(peak_int$cellID)))
    #message("   peak matrix = ", nrow(obj$pa), " | ", ncol(obj$pa))
    all.peaks <- rownames(obj$pa)

    # find peaks within the gene range
    message(" - building ACR x gene body weight matrix ...")
    upstream <- promoters(Ann, upstream=pRange, downstream=0)
    upstream_peaks <- suppressWarnings(findOverlaps(acrs, upstream, minoverlap=1,
                                   type=c("within"),
                                   select="all",
                                   ignore.strand=T))
    peak_genes <- data.frame(geneID=names(Ann)[upstream_peaks@to],
                             acrID=names(acrs)[upstream_peaks@from])
    all.genes <- unique(c(all.genes), as.character(peak_genes$geneID))
    all.peaks <- unique(c(all.peaks), as.character(peak_genes$acrID))
    peak_genes$distance <- distance(Ann[peak_genes$geneID,],
                                    acrs[peak_genes$acrID,])
    peak_genes$geneID <- factor(as.character(peak_genes$geneID), levels=all.genes)
    peak_genes$acrID <- factor(as.character(peak_genes$acrID), levels=all.peaks)
    peak_genes$weight <- 2^-(peak_genes$distance/con)
    peak_genes <- peak_genes[complete.cases(peak_genes),]
    obj$pgw <- sparseMatrix(i=as.numeric(peak_genes$geneID),
                            j=as.numeric(peak_genes$acrID),
                            x=peak_genes$weight,
                            dims=c(length(levels(peak_genes$geneID)),
                                   length(levels(peak_genes$acrID))),
                            dimnames=list(levels(peak_genes$geneID),
                                          levels(peak_genes$acrID)))

    shared <- intersect(rownames(obj$pa), colnames(obj$pgw))
    #message("   peak gene weighted matrix = ", nrow(obj$pgw), " | ", ncol(obj$pgw))
    #message("   number of shared peaks = ", length(shared))
    obj$peak_gene_score <- obj$pgw[,shared] %*% obj$pa[shared,]
    #message("   regulatory activity = ", nrow(obj$peak_gene_score), " | ", ncol(obj$peak_gene_score))
    #message(" - gene body accessibility")
    #print(head(obj$gba[,1:5]))
    #message(" - peak accessibility")
    #print(head(obj$pa[,1:5]))
    #message(" - regulatory activity")
    #print(head(obj$peak_gene_score[,1:5]))

    # aggregate regulatory score with gene body score
    shared.cells <- intersect(colnames(obj$peak_gene_score),
                              colnames(obj$gba))
    shared.genes <- intersect(rownames(obj$peak_gene_score),
                              rownames(obj$gba))
    obj$gene_activity <- obj$peak_gene_score[shared.genes,shared.cells] + obj$gba[shared.genes,shared.cells]

    # scale gene activity
    obj$scaled_gene_activity <- obj$gene_activity %*% Diagonal(x=per.cell.scale/Matrix::colSums(obj$gene_activity))
    
    # return
    return(obj)
}
