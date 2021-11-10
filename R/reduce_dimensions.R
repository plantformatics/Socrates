###################################################################################################
###################################################################################################
###################################################################################################
#' Reduce dimensions of residual matrix
#'
#' Run SVD on Pearson residuals using IRLBA.
#'
#' @importFrom irlba irlba
#' @importFrom RcppML nmf
#'
#' @param obj list, object containing "pearson_residuals" output by the function regModel.
#' @param method character, string denoting dimension reduction method. Can be one of "SVD" or
#' "NMF". Defaults to SVD.
#' @param n.pcs numeric, number of singular values to calculate.
#' @param scaleVars logical, whether or not to scale PCs by variance explained (or to scale
#' NMF components by scale factors). Default to TRUE.
#' @param num.vars number of highly variable ACRs/bins to use for dimensionality reduction.
#' Defaults to 5000. Set to NULL to use all ACRs/bins.
#' @param regNum number of peaks/bins to use for regularization. regNum must be equal or
#' less than num.vars. Defaults to 5000. 
#' @param cor.max float, maximum spearman correlation between log10nSites (log10 number of
#' accessible peaks) and singular value to keep. Singular values with correlations greater than
#' cor.max are removed. Ranges from 0 to 1. Default set to 0.75.
#' @param doL2 logical, whether or not to L2 normalize barcodes.
#' @param doL1 logical, whether or not to L1 normalize barcodes
#' @param stdLSI logical, whether or not to standardize barcodes.
#' @param refit_residuals logical, whether or not to use quasibinomial logistic regression residuals
#' for num.vars features. Only applicable when num.vars < nrow(obj$counts) and when the normalization
#' was performed with tfidf. Defaults to FALSE.
#' @param residuals_slotName character, character string of the desired residual slotName. Defaults
#' to "residuals".
#' @param svd_slotName character, character string for naming the SVD output in the returned
#' object. Defaults to "PCA".
#' @param verbose logical. Defaults to FALSE.
#' @param ... Additional arguments to be passed to RcppML::nmf
#'
#' @rdname reduceDims
#' @export
reduceDims <- function(obj,
                       method="SVD",
                       n.pcs=50,
                       scaleVar=T,
                       num.var=5000,
		               regNum=5000,
                       cor.max=0.75,
                       doL2=F,
                       doL1=F,
                       doSTD=T,
                       refit_residuals=F,
                       residuals_slotName="residuals",
                       svd_slotName="PCA",
                       verbose=FALSE,
                       ...){
    
    # sub functions
    l2norm <- function(x){x/sqrt(sum(x^2))}
    l1norm <- function(x){x/(sum(x))}
    RowVar <- function(x){
        spm <- t(x)
        stopifnot( methods::is( spm, "dgCMatrix" ) )
        ans <- sapply( base::seq.int(spm@Dim[2]), function(j) {
            if( spm@p[j+1] == spm@p[j] ) { return(0) } # all entries are 0: var is 0
            mean <- base::sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
            sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
                mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
        names(ans) <- spm@Dimnames[[2]]
        ans
    }
    
    # check if slotName exists
    if(is.null(obj[[residuals_slotName]])){
        message(" ERROR: useSlot - ", residuals_slotName, " - is missing ...")
        stop("exiting")
    }
    
    # if use subset
    if(!is.null(num.var)){
        if(!is.null(obj$meta$library) & length(unique(obj$meta$library)) > 1){
            df <- lapply(unique(obj$meta$library), function(z){
                ids <- rownames(subset(obj$meta, obj$meta$library==z))
                frac <- length(ids)/nrow(obj$meta)
                row.var <- RowVar(obj[[residuals_slotName]][,ids])
                row.means <- Matrix::rowMeans(obj[[residuals_slotName]][,ids])
                vals <- (loess(row.var~row.means)$residuals*frac)
                names(vals) <- names(row.means)
                return(vals)
            }
            df <- as.matrix(do.call(cbind, df))
            adj.row.var <- rowSums(df, na.rm=T)
            adj.row.var <- adj.row.var[order(adj.row.var, decreasing=T)]
            topSites <- names(head(adj.row.var, n=num.var))            
            
        }else{
            row.var <- RowVar(obj[[residuals_slotName]])
            row.means <- Matrix::rowMeans(obj[[residuals_slotName]])
            adj.row.var <- loess(row.var~row.means)$residuals
            names(adj.row.var) <- names(row.var)
            adj.row.var <- adj.row.var[order(adj.row.var, decreasing=T)]
            topSites <- names(head(adj.row.var, n=num.var))
        }
        
        # input matrix
        if(refit_residuals & obj$norm_method =="tfidf"){
            if(regNum > length(topSites)){
                regNum <- length(topSites)
            }
            test.dat <- list(counts=obj$counts[topSites,], meta=obj$meta)
            M <- regModel(test.dat, subpeaks=regNum, verbose=verbose)$residuals
            M <- Matrix(t(apply(M, 1, function(x){x - min(x, na.rm=T)})), sparse=T)
            M <- M[Matrix::rowSums(M) > 0,]
        }else{
            M <- obj$residuals[topSites,]
            if(obj$norm_method != "tfidf"){
                M <- t(apply(M, 1, function(x){
                    x - min(x, na.rm=T)
                }))
            }
            M <- M[Matrix::rowSums(M) > 0,]
        }
    }else{
        M <- obj[[residuals_slotName]]
    }
    
    # choose method
    if(method=="SVD"){
        
        # set method  used
        obj$rdMethod <- "SVD"
        
        # verbose
        if(verbose){message(" - reduce dimensions with SVD ... ")}
        
        # run
        pcs <- irlba(t(M), nv=n.pcs)
        pc <- pcs$u
        
        # scale by % variance
        if(scaleVar){
            pc <- pc %*% diag(pcs$d)
        }
        pc[is.na(pc)] <- 0
        
    }else if(method=="NMF"){
        
        # set method used
        obj$rdMethod <- "NMF"
        
        # if use subset
        if(verbose){message(" - running NMF...")}
        
        # run NMF
        pcs <- RcppML::nmf(M, n.pcs, verbose=verbose, ...)
        pcs$u <- t(pcs$h)
        pcs$v <- pcs$w
        if(scaleVar){
            pc <- t(pcs$h) %*% Diagonal(x=1/pcs$d)
        }else{
            pc <- t(pcs$h)
        }
        pc[is.na(pc)] <- 0
        
    }
    
    # add colnames
    rownames(pc) <- colnames(obj[[residuals_slotName]])
    colnames(pc) <- paste0("PC_", seq(1:ncol(pc)))
    
    # remove PCs with correlation to read-depth
    if(verbose){message(" - removing components correlated to read depth...")}
    if(cor.max < 1){
        depth <- Matrix::colSums(obj$counts[,rownames(pc)])
        cors <- apply(pc, 2, function(u) cor(u,depth,method="spearman"))
        cors <- abs(cors)
        idx.keep <- cors < cor.max
        pc <- pc[,idx.keep]
    }else{
        idx.keep <- seq(1:ncol(pc))
    }
    
    # standardize, L2, or L1 reduced dimensions per cell
    if(verbose){message(" - normalizing reduced dimensions...")}
    if(doL2){
        pc <- t(apply(pc, 1, l2norm))
    }else if(doL1){
        pc <- t(apply(pc, 1, l1norm))
    }else if(doSTD){
        pc <- t(apply(pc, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
    }
    pc[is.na(pc)] <- 0

    # return
    obj[[svd_slotName]] <- pc
    model_slot <- paste0(svd_slotName, "_model")
    obj[[model_slot]] <- list()
    obj[[model_slot]]$d <- pcs$d
    obj[[model_slot]]$v <- pcs$v
    obj[[model_slot]]$u <- pcs$u
    obj[[model_slot]]$keep_pcs <- idx.keep
    obj$hv_sites <- rownames(M)
    return(obj)
    
    
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Project cells into a reduced embedding with UMAP
#'
#' @importFrom uwot umap
#'
#' @param obj list, object containing 'PCA' for projecting into a reduced embedding with UMAP.
#' @param m.dist numeric, m_dist parameter for uwot::umap. Defaults to 0.1
#' @param k.near numeric, k-nearest neighbors used by the uwot::umap algorithm. Defaults to 15.
#' @param metric character, distance metric used by uwot::umap. Defaults to cosine.
#' @param svd_slotName character, name of desired svd_slotName for running UMAP. Defaults to "PCA".
#' @param umap_slotName character, name of desired umap_slotName for return UMAP results. Defaults
#' to "UMAP".
#' @param verbose, logical. Defaults to FALSE.
#'
#' @rdname projectUMAP
#' @export
projectUMAP <- function(obj,
                        m.dist=0.01,
                        k.near=40,
                        metric="cosine",
                        svd_slotName="PCA",
                        umap_slotName="UMAP",
                        verbose=FALSE,
                        seed=1){

    # checks
    if(is.null(obj[[svd_slotName]])){
        message(" - ERROR: ", svd_slotName, " slot is not present in object")
        stop()
    }

    # verbose
    if(verbose){message(" - non-linear dimensionality reduction with UMAP ...")}

    # run UMAP (uwot implementation)
    set.seed(seed)
    umap.res <- uwot::umap(obj[[svd_slotName]],
                           verbose=verbose,
                           min_dist=m.dist,
                           n_neighbors=k.near,
                           metric=metric,
                           ret_model=T)

    # convert to data frame
    umap.out <- as.data.frame(umap.res$embedding)
    rownames(umap.out) <- rownames(obj[[svd_slotName]])
    colnames(umap.out) <- c("umap1", "umap2")

    # return object
    obj$meta$umap1 <- umap.out[rownames(obj$meta),]$umap1
    obj$meta$umap2 <- umap.out[rownames(obj$meta),]$umap2
    obj[[umap_slotName]] <- umap.out
    obj$UMAP_model <- umap.res
    return(obj)
}
