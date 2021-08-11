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
#' @param nTrials Numeric. Number of times nSample cells are simulated.
#' @param nSample Numeric. Number of synthetic doublets to create per trial.
#' @param k Numeric. Number of nearest neighbors to search for when estimating doublet enrichment via knn
#' from the narbor package.
#' @param n.pcs Numeric. Number of PCS/SVD components to retain when reducing dimensions.
#' @param threads Numeric. Number of threads to use for mclapply from the parallel package.
#'
#' @rdname detectDoublets
#' @export
#'
detectDoublets <- function(obj=NULL,
                           nTrials=5,
                           nSample=1000,
                           k=10,
                           n.pcs=50,
                           threads=1){
    
    # pre-checks
    if(is.null(obj)){
        stop("! Socrates object is required ...")
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
            sobj <- regModel(sobj, nthreads=threads)
        }else if(normModel == 'regModel2'){
            sobj <- regModel2(sobj, nthreads=threads)
        }else if(normModel == 'tfidf'){
            sobj <- tfidf(sobj)
        }else if(normModel == 'logisticModel'){
            sobj <- logisticModel(sobj, nthreads=threads)
        }
        
        # get V matrix
        V <- Matrix::t(sobj$residuals) %*% v %*% Matrix::diag(1/d)
        
        # diagonal
        if(n.pcs > ncol(u)){
            n.pcs <- ncol(u)
        }
        svdDiag <- matrix(0, nrow=n.pcs, ncol=n.pcs)
        diag(svdDiag) <- d
        matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
        matSVD <- as.matrix(matSVD)
        rownames(matSVD) <- colnames(sobj$residuals)
        colnames(matSVD) <- paste0("PC_",seq_len(ncol(matSVD)))
        matSVD <- matSVD[,idx.keep]
        return(matSVD)
        
    }
    
    
    ##############################################################################
    # Simulate doublets
    ##############################################################################
    
    # specify objects
    sampleRatio1 <- c(0.5)
    sampleRatio2 <- c(0.5)
    mat <- obj$counts
    
    # scale nTrials by number of cells
    nTrials <- nTrials * max( floor(nrow(obj$meta) / nSample), 1 )
    
    # simulated matrix
    set.seed(1)
    message(" - Creating synthetic doublets ...")
    simMat <- mclapply(seq_len(nTrials), function(y){
        outs <- lapply(seq_along(sampleRatio1), function(x){
            idx1 <- sample(seq_len(ncol(mat)), nSample, replace = TRUE)
            idx2 <- sample(seq_len(ncol(mat)), nSample, replace = TRUE)
            simulatedMat <- .sampleSparseMat(mat = mat[,idx1], sampleRatio = sampleRatio1[x]) +
                .sampleSparseMat(mat = mat[,idx2], sampleRatio = sampleRatio2[x])
            simulatedMat@x[simulatedMat@x > 0] <- 1
            b <- data.frame(cellIDs=paste0("sim.",y,'.',seq(1:ncol(simulatedMat))),
                            row.names=paste0("sim.",y,'.',seq(1:ncol(simulatedMat))))
            b$nSites <- Matrix::colSums(simulatedMat)
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
    message(" - Created ", nrow(simMat), " synthetic doublets ...")
    
    # project original
    ogMat <- .projectSVD(obj,
                         u=obj$SVD_model$u,
                         v=obj$SVD_model$v,
                         d=obj$SVD_model$d,
                         n.pcs=n.pcs,
                         idx.keep=obj$SVD_model$keep_pcs,
                         normModel=obj$norm_method)
    
    # merge
    allMat <- rbind(simMat, ogMat)
    #allMat <- t(apply(allMat, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
    num.cells <- nrow(ogMat)
    rm(ogMat)
    
    # run UMAP model
    set.seed(1)
    proj.umap <- uwot::umap_transform(X = as.matrix(allMat),
                                      model = obj$UMAP_model)
    rownames(proj.umap) <- rownames(allMat)
    colnames(proj.umap) <- c("umap1", "umap2")
    out <- list()
    
    
    ##############################################################################
    # Compute Doublet Scores from SVD Embedding
    ##############################################################################
    message(" - Computing KNN doublets (SVD)...")
    knnDoub <- .computeKNN(allMat[-seq_len(nrow(simMat)),], allMat[seq_len(nrow(simMat)),], k=k)
    knnDoub.svd <- knnDoub
    
    #Compile KNN Sums
    countKnn <- rep(0, num.cells)
    names(countKnn) <- rownames(allMat[-seq_len(nrow(simMat)),])
    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <- countKnn[as.integer(names(tabDoub))] + tabDoub
    scaleTo <- 10000
    scaleBy <- scaleTo / nrow(num.cells)
    countknn.svd <- countKnn
    
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
    knnDoub <- .computeKNN(proj.umap[-seq_len(nrow(simMat)),], proj.umap[seq_len(nrow(simMat)),], k=k)
    knnDoub.umap <- knnDoub
    
    #Compile KNN Sums
    countKnn <- rep(0, num.cells)
    names(countKnn) <- rownames(allMat[-seq_len(nrow(simMat)),])
    tabDoub <- table(as.vector(knnDoub))
    countKnn[as.integer(names(tabDoub))] <- countKnn[as.integer(names(tabDoub))] + tabDoub
    scaleTo <- 10000
    scaleBy <- scaleTo / num.cells
    countknn.umap <- countKnn
    
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
    out$knnDoublet.umap <- knnDoub.umap
    out$knnDoublet.svd <- knnDoub.svd
    out$countknn.umap <- countknn.umap
    out$countknn.svd <- countknn.svd
    obj$doublets <- out
    obj$meta$doubletscore <- obj$doublets$doubletEnrichUMAP[rownames(obj$meta)]
    return(obj)
    
}



###################################################################################################
###################################################################################################
###################################################################################################
#' filterDoublets
#'
#' This function identifies potential doublets. A new column in the meta slot (d.type) denotes
#' the doublet status. Doublets can also be directly filtered from the data set by setting
#' removeDoublets to TRUE.
#'
#'
#' @param obj Socrates object. For best results, use downstream of the function cleanCells prior to
#' normalization.
#' @param filterRatio Numeric. Defaults to 1.5.
#' @param embedding Character string. Which embedding to use for doublet estimation filtering. Defaults
#' to "UMAP".
#' @param libraryVar Character string. Meta data column containing library information. Defaults to
#' NULL, which is specifically for single library/replicate objects. Users with multiple libraries/replicates
#' in the object must set the libraryVar parameter to the column name specifying each cells library ID. 
#' The default column name used by mergeSocratesRDS is "sampleID". 
#' @param verbose Boolean. Defaults to TRUE.
#' @param umap_slotname Character string. Defaults to "UMAP".
#' @param svd_slotname Character string. Defaults to "PCA".
#' @param removeDoublets Boolean. Defaults to FALSE.
#'
#' @rdname filterDoublets
#' @export
#'
filterDoublets <- function(obj=NULL, filterRatio=1.5, embedding="UMAP", libraryVar="sampleID",
                           verbose=T, umap_slotname="UMAP", svd_slotname="PCA", removeDoublets=F){
    
    # split by library/replicate
    o.ids <- rownames(obj$meta)
    if(is.null(libraryVar)){
        libraryVar <- "temp"
        obj$meta[,libraryVar] <- 1
    }else if(! libraryVar %in% colnames(obj$meta) & !is.null(libraryVar)){
        message(" !! The parameter `libraryVar`: ", libraryVar, ", does not exist as a column in the meta data slot ")
        stop(" - please make sure that the parameter `libraryVar` is correctly set or is set to NULL for single-library objects ...")
    }
    obj$meta[,libraryVar] <- as.character(obj$meta[,libraryVar])
    libs <- unique(as.character(obj$meta[,libraryVar]))
    
    # iterate over each library/replicate
    outs <- lapply(libs, function(z){
        
        # subset by library/replicate
        obj.lib <- subset(obj$meta, obj$meta[,libraryVar]==z)
        num.cells <- nrow(obj.lib)
        ids <- rownames(obj.lib)
        
        # estimate # cells to remove
        rm.counts <- floor(filterRatio * ((num.cells^2) / 100000))
        keep.count <- num.cells - rm.counts
        
        # select embedding
        if(embedding == "UMAP"){
            obj.lib$doubletscore <- obj$doublets$doubletEnrichUMAP[rownames(obj.lib)]
        }else{
            obj.lib$doubletscore <- obj$doublets$doubletEnrichSVD[rownames(obj.lib)]
        }
        
        # reorder and call doublets
        obj.lib <- obj.lib[order(obj.lib$doubletscore, decreasing=T),]
        obj.lib$d.type <- c(rep("doublet", rm.counts), rep("singlet", keep.count))
        obj.lib <- obj.lib[ids,]
        return(obj.lib)
    })
    outs <- do.call(rbind, outs)
    
    # filter doublets
    if(removeDoublets){
        if(verbose){message("   * Doublet filtering * Input: cells = ", ncol(obj$counts), " | peaks = ", nrow(obj$counts))}
        obj$meta <- subset(outs, outs$d.type == "singlet")
        obj$counts <- obj$counts[,rownames(obj$meta)]
        obj$counts <- obj$count[Matrix::rowSums(obj$counts)>0,]
        obj$counts <- obj$count[,Matrix::colSums(obj$counts)>0]
        obj[[umap_slotname]] <- obj[[umap_slotname]][colnames(obj$counts),]
        obj[[svd_slotname]] <- obj[[svd_slotname]][colnames(obj$counts),]
        new.num.cells <- ncol(obj$counts)
        ratio <- signif(((1-new.num.cells/num.cells)*100), digits=3)
        if(verbose){message("   * Doublet filtering * Filtered (", ratio, "%) : cells = ", ncol(obj$counts), " | peaks = ", nrow(obj$counts))}
        
    }else{
        obj$meta <- outs[o.ids,]
    }
    
    # remove temp column if necessary
    if(is.null(libraryVar)){
        obj$meta$temp <- NULL
    }
    
    # return object
    return(obj)
}


