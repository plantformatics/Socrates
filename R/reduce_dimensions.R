###################################################################################################
###################################################################################################
###################################################################################################
#' Reduce dimensions of residual matrix
#'
#' Run SVD on Pearson residuals using IRLBA.
#'
#' @importFrom irlba irlba
#'
#' @param obj list, object containing "pearson_residuals" output by the function regModel.
#' @param n.pcs numeric, number of singular values to calculate.
#' @param scaleVars logical, whether or not to scale PCs by variance explained. Default to TRUE.
#' @param cor.max float, maximum spearman correlation between log10nSites (log10 number of
#' accessible peaks) and singular value to keep. Singular values with correlations greater than
#' cor.max are removed. Ranges from 0 to 1. Default set to 0.75.
#' @param doL2 numeric, whether to L2 normalize barcodes (doL2=1), L2 normalize PCs (doL2=2),
#' or no L2 normalization (NULL or FALSE).
#' @param stdLSI numeric, whether to standarize barcodes (stdLSI=1), standardize PCs (stdLSI=2),
#' or no SVD standarization (NULL or FALSE). Defaults to 1 (standardize components for each
#' barcode).
#' @param residuals_slotName character, character string of the desired residual slotName. Defaults
#' to "residuals".
#' @param svd_slotName character, character string for naming the SVD output in the returned
#' object. Defaults to "PCA".
#' @param verbose logical. Defaults to FALSE.
#'
#' @rdname reduceDims
#' @export
reduceDims <- function(obj,
                       n.pcs=50,
                       scaleVar=T,
                       cor.max=0.75,
                       doL2=NULL,
                       stdLSI=1,
                       residuals_slotName="residuals",
                       svd_slotName="PCA",
                       verbose=FALSE){

    # verbose
    if(verbose){message(" - reduce dimensions with SVD ... ")}

    # check if slotName exists
    if(is.null(obj[[residuals_slotName]])){
        message(" ERROR: useSlot - ", residuals_slotName, " - is missing ...")
        stop("exiting")
    }

    # run
    pcs <- irlba(t(obj[[residuals_slotName]]), nv=n.pcs)
    pc <- pcs$u

    # scale by % variance
    if(scaleVar){
        pc <- pc %*% diag(pcs$d)
    }

    # add colnames
    rownames(pc) <- colnames(obj[[residuals_slotName]])
    colnames(pc) <- paste0("PC_", seq(1:ncol(pc)))

    # remove PCs with correlation to read-depth
    if(cor.max < 1){
        depth <- Matrix::colSums(obj$counts[,rownames(pc)])
        cors <- apply(pc, 2, function(u) cor(u,depth,method="spearman"))
        idx.keep <- abs(cors) < cor.max
        pc <- pc[,idx.keep]
    }

    # convert l2 norm
    if(!is.null(doL2)){
        if(doL2 != FALSE){

            # verbose
            if(verbose){message(" - L2 norm of reduced dimensions ... ")}

            if(doL2==1){
                pc <- t(apply(pc, 1, function(x) x/(sqrt(sum(x^2)))))
            }else if(doL2==2){
                pc <- apply(pc, 2, function(x) x/(sqrt(sum(x^2))))
            }

            if(cor.max < 1){
                depth <- Matrix::colSums(obj$counts[,rownames(pc)])
                cors <- apply(pc, 2, function(u) cor(u,depth,method="spearman"))
                pc <- pc[,abs(cors) < cor.max]
            }
        }
    }

    # standarize reduced dimensions per cell
    if(!is.null(stdLSI)){
        if(stdLSI != FALSE){

            # verbose
            if(stdLSI==1){
                if(verbose){message(" - standardizing reduced dimensions per cell ... ")}
                pc <- t(apply(pc, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
            }else if(stdLSI==2){
                if(verbose){message(" - standardizing reduced dimensions by components ... ")}
                pc <- apply(pc, 2, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)})
            }
        }
    }

    # return
    obj[[svd_slotName]] <- pc
    obj$SVD_model <- list()
    obj$SVD_model$d <- pcs$d
    obj$SVD_model$v <- pcs$v
    obj$SVD_model$u <- pcs$u
    obj$SVD_model$keep_pcs <- idx.keep
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
