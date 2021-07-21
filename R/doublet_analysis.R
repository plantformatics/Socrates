###################################################################################################
###################################################################################################
###################################################################################################
#' detectDoublets
#'
#' This function estimates doublet likelihood via synthetic doublet creation and enrichment 
#' analysis. The underlying code is heavily adapted from ArchR described by Granja and Corces et al. 
#' If you use detectDoublets in your analysis, please cite the original ArchR paper. 
#'
#' @importFrom nabor knn
#' @importFrom parallel mclapply
#'
#' @param obj Socrates object. For best results, use downstream of the function cleanCells prior to 
#' normalization. 
#' @param normModel Character string. normalization model for embedding synthetic doublets. Options include
#' 'regModel', 'regModel2', 'logisticModel', and 'tfidf'. TF-IDF is the default due to its quicker run time.
#' But the other models are also available if the user wish to maintain a consistent normalization scheme 
#' as the intended analysis. 
#' @param nTrials Numeric. Number of times nSample cells are simulated. 
#' @param nSample Numeric. Number of synthetic doublets to create per trial.
#' @param k Numeric. Number of nearest neighbors to search for when estimating doublet enrichment via knn 
#' from the narbor package. 
#' @param n.pcs Numeric. Number of PCS/SVD components to retain when reducing dimensions. 
#' @param seed Numeric. 
#' @param threads Numeric. Number of threads to use for mclapply from the parallel package.
#'  
#'
#' @rdname detectDoublets
#' @export
#'
detectDoublets <- function(obj=NULL,
                           nTrials=5, 
                           nSample=1000, 
                           k=10, 
                           n.pcs=50,
                           seed=1, 
                           threads=1,
                           svd_slot="PCA"){
    
    # pre-checks
    if(is.null(obj)){
        stop("! Socrates object is required ...")
    }
    if(is.null(obj[[svd_slot]])){
        stop("! Socrates object", svd_slot, " is empty. Please run reduceDims first or check that the specified slot name is correct ...")
    }
    
    # hidden functions
    .sampleSparseMat <- function(mat = NULL, sampleRatio = 0.5){
        total <- length(mat@x)
        sampleTo <- floor(total * (1-sampleRatio))
        mat@x[sample(seq_len(total), sampleTo)] <- 0
        mat <- drop0(mat)
        mat
    }
    .computeKNN <- function(data=NULL, 
                            query=NULL, 
                            k=50, 
                            includeSelf=FALSE, 
                            ...){
        if(is.null(query)){
            query <- data
            searchSelf <- TRUE
        }else{
            searchSelf <- FALSE
        }
        if(searchSelf & !includeSelf){
            knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
            knnIdx <- knnIdx[,-1,drop=FALSE]
        }else{
            knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
        }
        knnIdx
    }
    .projectSVD <- function(sobj,
                            u=NULL,
                            v=NULL,
                            d=NULL,
                            n.pcs=NULL,
                            idx.keep=NULL,
                            normModel="regModel"){
        
        # run normalization
        if(normModel == 'regModel'){
            sobj <- regModel(sobj, nthreads=threads, make.sparse=T)
        }else if(normModel == 'regModel2'){
            sobj <- regModel2(sobj, nthreads=threads, make.sparse=T)
        }else if(normModel == 'tfidf'){
            sobj <- tfidf(sobj)
        }else if(normModel == 'logisticModel'){
            sobj <- logisticModel(sobj, nthreads=threads)
        }
        mat <- sobj$residuals
        
        # get V matrix
        V <- Matrix::t(mat) %*% v %*% Matrix::diag(1/d)
        
        # diagonal
        if(n.pcs > length(d)){
            n.pcs <- length(d)
        }
        svdDiag <- matrix(0, nrow=n.pcs, ncol=n.pcs)
        diag(svdDiag) <- d
        matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
        matSVD <- as.matrix(matSVD)
        rownames(matSVD) <- colnames(mat)
        colnames(matSVD) <- paste0("PC_",seq_len(ncol(matSVD)))
        matSVD <- matSVD[,idx.keep]
        return(matSVD)
        
    }
    
    # specify objects
    sampleRatio1 <- c(0.5)    
    sampleRatio2 <- c(0.5)
    mat <- obj$counts
    
    # simulated matrix
    set.seed(seed)
    message(" - Creating synthetic doublets ...")
    simMat <- mclapply(seq_len(nTrials), function(y){
            if(y %% 5 == 0){
                message("   ~ iterated over ", y, " of ", nTrials, " trials...")
                gc()
            }
            outs <- lapply(seq_along(sampleRatio1), function(x){
                idx1 <- sample(seq_len(ncol(mat)), nSample, replace = TRUE)
                idx2 <- sample(seq_len(ncol(mat)), nSample, replace = TRUE)
                simulatedMat <- .sampleSparseMat(mat = mat[,idx1], sampleRatio = sampleRatio1[x]) + 
                    .sampleSparseMat(mat = mat[,idx2], sampleRatio = sampleRatio2[x])
                b <- data.frame(cellIDs=paste0("sim.",y,'.',seq(1:ncol(simulatedMat))), 
                                row.names=paste0("sim.",y,'.',seq(1:ncol(simulatedMat))))
                b$nSites   <- Matrix::colSums(simulatedMat)
                b$log10nSites <- log10(b$nSites)
                colnames(simulatedMat) <- rownames(b)
                rownames(simulatedMat) <- rownames(mat)
                simobj <- list(counts=simulatedMat, meta=b)
                simSVD <- .projectSVD(simobj, 
                                      u=obj$SVD_model$u, 
                                      v=obj$SVD_model$v, 
                                      d=obj$SVD_model$d, 
                                      n.pcs=n.pcs,
                                      idx.keep=obj$SVD_model$keep_pcs,
                                      normModel=obj$norm_method)
                return(simSVD)
            })
            outs <- do.call(rbind, outs)
        }, mc.cores = threads)
    simMat <- do.call(rbind, simMat)
    simMat <- as.matrix(simMat)
    simMat <- t(apply(simMat, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
    message(" - Created ", nrow(simMat), " synthetic doublets ...")
    
    # merge
    allMat <- rbind(simMat, obj[[svd_slot]])
    num.cells <- nrow(obj[[svd_slot]])
    
    # run UMAP model
    proj.umap <- uwot::umap_transform(X = as.matrix(allMat), 
                                      model = obj$UMAP_model)
    rownames(proj.umap) <- rownames(allMat)
    colnames(proj.umap) <- c("umap1", "umap2")
    out <- list()
    
    ##############################################################################
    # Compute Doublet Scores from SVD Embedding
    ##############################################################################
    message(" - Computing KNN doublets (SVD)...")
    knnDoub <- .computeKNN(allMat[-seq_len(nrow(simMat)),], allMat[seq_len(nrow(simMat)),], k)
        
    #Compile KNN Sums
    countKnn <- rep(0, num.cells)
    names(countKnn) <- rownames(allMat[-seq_len(nrow(simMat)),])
    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <- countKnn[as.integer(names(tabDoub))] + tabDoub
    scaleTo <- 10000
    scaleBy <- scaleTo / num.cells
    
    #P-Values
    pvalBinomDoub <- unlist(lapply(seq_along(countKnn), function(x){
        countKnnx <- round(countKnn[x] * scaleBy)
        sumKnnx <- round(sum(countKnn) * scaleBy)
        pbinom(countKnnx - 1, sumKnnx, 1 / scaleTo, lower.tail = FALSE)
    }))
    
    #Adjust
    padjBinomDoub <- p.adjust(pvalBinomDoub, method = "bonferroni")
    
    #Convert To Scores
    doubletScore <- -log10(pmax(padjBinomDoub, 4.940656e-324))
    doubletEnrich <- (countKnn / sum(countKnn)) / (1 / num.cells)
    doubletEnrich <- 10000 * doubletEnrich / length(countKnn) #Enrichment Per 10000 Cells in Data Set  
    
    #Store Results
    out$doubletEnrichSVD <- doubletEnrich
    out$doubletScoreSVD <- doubletScore
    
    ##############################################################################
    # Compute Doublet Scores from UMAP Embedding
    ##############################################################################
    message(" - Computing KNN doublets (UMAP)...")
    knnDoub <- .computeKNN(proj.umap[-seq_len(nrow(simMat)),], proj.umap[seq_len(nrow(simMat)),], k)
    
    #Compile KNN Sums
    countKnn <- rep(0, num.cells)
    names(countKnn) <- rownames(allMat[-seq_len(nrow(simMat)),])
    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <- countKnn[as.integer(names(tabDoub))] + tabDoub
    scaleTo <- 10000
    scaleBy <- scaleTo / num.cells
    
    #P-Values
    pvalBinomDoub <- unlist(lapply(seq_along(countKnn), function(x){
        countKnnx <- round(countKnn[x] * scaleBy)
        sumKnnx <- round(sum(countKnn) * scaleBy)
        pbinom(countKnnx - 1, sumKnnx, 1 / scaleTo, lower.tail = FALSE)
    }))
    
    #Adjust
    padjBinomDoub <- p.adjust(pvalBinomDoub, method = "bonferroni")
    
    #Convert To Scores
    doubletScore <- -log10(pmax(padjBinomDoub, 4.940656e-324))
    doubletEnrich <- (countKnn / sum(countKnn)) / (1 / num.cells)
    doubletEnrich <- 10000 * doubletEnrich / length(countKnn) #Enrichment Per 10000 Cells in Data Set  
    
    #Store Results
    out$doubletEnrichUMAP <- doubletEnrich
    out$doubletScoreUMAP <- doubletScore
    out$dSVD <- allMat
    out$dUMAP <- proj.umap
    obj$doublets <- out
    
    return(obj)
        
}



###################################################################################################
###################################################################################################
###################################################################################################
#' filterDoublets
#'
#' This function estimates doublet likelihood via synthetic doublet creation and enrichment 
#' analysis. The underlying code is heavily adapted from ArchR described by Granja and Corces et al. 
#' If you use detectDoublets in your analysis, please cite the original ArchR paper. 
#'
#'
#' @param obj Socrates object. For best results, use downstream of the function cleanCells prior to 
#' normalization. 
#' @param filterRatio Numeric.
#' @param embedding Character.
#' @param verbose Boolean.
#'
#' @rdname filterDoublets
#' @export
#'
filterDoublets <- function(obj=NULL, filterRatio=1.5, embedding="UMAP", verbose=T){
    
    # number of cells
    num.cells <- nrow(obj$meta)
    
    # estimate # cells to remove
    rm.counts <- floor(filterRatio * ((num.cells^2) / 100000))
    keep.count <- num.cells - rm.counts
    
    # remove 
    if(embedding == "UMAP"){
        obj$meta$doubletscore <- obj$doublets$doubletEnrichUMAP[rownames(obj$meta)]
    }else{
        obj$meta$doubletscore <- obj$doublets$doubletEnrichSVD[rownames(obj$meta)]
    }
    if(verbose){message("   * Doublet filtering * Input: cells = ", ncol(obj$counts), " | peaks = ", nrow(obj$counts))}
    obj$meta <- obj$meta[order(obj$meta$doubletscore, decreasing=F),]
    obj$meta <- head(obj$meta, n=keep.count)
    obj$counts <- obj$counts[,rownames(obj$meta)]
    obj$counts <- obj$count[Matrix::rowSums(obj$counts)>0,]
    obj$counts <- obj$count[,Matrix::colSums(obj$counts)>0]
    new.num.cells <- ncol(obj$counts)
    ratio <- signif(((1-new.num.cells/num.cells)*100), digits=3)
    if(verbose){message("   * Doublet filtering * Filtered (", ratio, "%) : cells = ", ncol(obj$counts), " | peaks = ", nrow(obj$counts))}
    return(obj)
}


