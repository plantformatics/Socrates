###################################################################################################
###################################################################################################
###################################################################################################
#' coAccess
#'
#' This function identifies co-accessible ACRs within clusters or across the entire data set
#'
#' @import Matrix
#' @import cicero
#' @import parallel
#'
#' @param obj Socrates object. Required.
#' @param genome path to genome file. See 'inst/extdata/genome' for an example. Required.
#' @param byGroup logical. Whether to run co-accessibility tests for each factor in 'groupID'.
#' Defaults to False.
#' @param groupID character. Column name in obj$meta for partitioning during co-accessibility
#' analysis. Defaults to "LouvainClusters".
#' @param eFDR logical. Whether or not to use empirical FDR to filter co-accessible ACRs. Defaults
#' to FALSE. Setting eFDR to 'TRUE' doubles memory usage.
#' @param fdr_thresh float. fdr threshold when using empirical FDR. Defaults to 0.05.
#' @param win_size numeric. Genomic distance to constrain co-accessible ACRs. Defaults to 500000.
#' @param sample_num numeric. Number of random windows to sample for determining distance_parameter
#'  during 'run_cicero'. See cicero for more details.
#' @param k numeric. Number of nearest neighbors for aggregating cells. Defaults to 30.
#' @param conn_slotName character. Slot name for depositing co-accessible ACRs. Defaults to coACRs.
#' @param nthreads numeric. Number of threads to use for estimating co-ACRs for each cluster or
#' celltype. Defaults to 1.
#' @param verbose logical. Defaults to FALSE.
#'
#' @rdname coAccess
#' @export
#'
coAccess <- function(obj,
                     genome=NULL,
                     byGroup=FALSE,
                     groupID="LouvainClusters",
                     eFDR=FALSE,
                     fdr_thresh=0.05,
                     win_size=500000,
                     sample_num=100,
                     k=30,
                     conn_slotName="coACRs",
                     nthreads=1,
                     verbose=FALSE){

    # hidden functions
    .loadMeta <- function(cds, meta, clusterID){

        # svd and raw
        ids <- colnames(cds)
        meta <- meta[ids,]
        meta$UMAP1 <- meta$umap1
        meta$UMAP2 <- meta$umap2
        ids.2 <- rownames(meta)
        cds <- cds[,colnames(cds) %in% ids.2]
        umap.d <- t(meta[,c("UMAP1","UMAP2")])

        # UMAP output
        cds@reducedDimA <- umap.d
        cds@reducedDimS <- umap.d
        cds@dim_reduce_type <- "UMAP"
        pData(cds)$Cluster <- as.factor(meta[,clusterID])

        # return data
        return(cds)
    }

    # check assumptions
    if(is.null(genome)){
        stop("! argument 'genome' is required, exiting...")
    }
    if(byGroup==T){
        if(is.null(groupID)){
            stop(" ! argument groupID is required when 'byGroup' is set to TRUE, exiting...")
        }
        if(!groupID %in% colnames(obj$Clusters)){
            stop(" ! argument groupID must be a column name in the 'meta' slot of obj, exiting...")
        }
    }

    # load genome
    genome <- read.table(genome)

    # create cicero objects and process
    obs <- as.data.frame(summary(obj$counts))
    obs$i <- as.factor(rownames(obj$counts)[obs$i])
    obs$j <- as.factor(colnames(obj$counts)[obs$j])
    cds <- make_atac_cds(obs, binarize=T)
    cds <- cds[,colnames(cds) %in% rownames(obj$Clusters)]
    pData(cds) <- obj$Clusters[colnames(exprs(cds)),]
    cds <- cds[Matrix::rowSums(exprs(cds))>0,]
    cds <- cds[,Matrix::colSums(exprs(cds))>0]
    cds <- detectGenes(cds)
    cds <- estimateSizeFactors(cds)
    cds <- .loadMeta(cds, obj$Clusters, groupID)

    # run per cluster or not
    if(eFDR==FALSE){

        if(byGroup==TRUE){

            # get cluster ids
            clusts <- unique(obj$Clusters[,groupID])

            # run in parallel
            outs <- mclapply(clusts, function(z){
                cluster.ids <- rownames(subset(obj$Clusters, obj$Clusters[,groupID]==z))
                sub.cds <- cds[,colnames(exprs(cds)) %in% cluster.ids]
                sub.cds <- sub.cds[Matrix::rowSums(exprs(sub.cds))>0,]
                sub.cds <- sub.cds[,Matrix::colSums(exprs(sub.cds))>0]
                sub.conns <- runCicero(sub.cds, genome=genome, k=k, win=win_size, sample_num=sample_num)
                if(verbose){message(" - found ", nrow(sub.conns), " potential co-accessible ACRs in group ",z, "...")}
                sub.conns$group <- z
                return(sub.conns)
            }, mc.cores=nthreads)

            # join into a single df
            outs <- as.data.frame(do.call(rbind, outs))
            obj[[conn_slotName]] <- outs

        }else{
            obj[[conn_slotName]] <- runCicero(cds, genome=genome, k=k, win=win_size, sample_num=sample_num)
            if(verbose){message(" - found ", nrow(obj[[conn_slotName]]), " potential co-accessible ACRs ...")}
        }

    }else{

        # shuffle
        shuf <- obs
        shuf$i <- shuf$i[sample(length(shuf$i))]
        shuf$j <- shuf$j[sample(length(shuf$j))]
        shufcds <- make_atac_cds(shuf, binarize=T)

        # add metadata and filter
        shufcds <- shufcds[,colnames(shufcds) %in% rownames(obj$Clusters)]
        pData(shufcds) <- obj$Clusters[colnames(exprs(shufcds)),]
        shufcds <- shufcds[Matrix::rowSums(exprs(shufcds))>0,]
        shufcds <- shufcds[,Matrix::colSums(exprs(shufcds))>0]

        # process
        shufcds <- detectGenes(shufcds)
        shufcds <- estimateSizeFactors(shufcds)

        # update cluster/UMAP info
        shufcds <- .loadMeta(shufcds, obj$Clusters, groupID)

        # run cicero by cluster or no
        if(byGroup==TRUE){

            # get cluster ids
            clusts <- unique(obj$Clusters[,groupID])

            # run in parallel
            outs <- mclapply(clusts, function(z){
                cluster.ids <- rownames(subset(obj$Clusters, obj$Clusters[,groupID]==z))
                sub.cds <- cds[,colnames(exprs(cds)) %in% cluster.ids]
                sub.cds <- sub.cds[Matrix::rowSums(exprs(sub.cds))>0,]
                sub.cds <- sub.cds[,Matrix::colSums(exprs(sub.cds))>0]
                sub.shufcds <- shufcds[,colnames(exprs(shufcds)) %in% cluster.ids]
                sub.shufcds <- sub.shufcds[Matrix::rowSums(exprs(sub.shufcds))>0,]
                sub.shufcds <- sub.shufcds[,Matrix::colSums(exprs(sub.shufcds))>0]

                # run
                sub.conns <- runCicero(sub.cds, genome=genome, k=k, win=win_size, sample_num=sample_num)
                if(verbose){message(" - found ", nrow(sub.conns), " potential co-accessible ACRs in original group ",z, "...")}
                sub.shufconns <- runCicero(sub.shufcds, genome=genome, k=k, win=win_size, sample_num=sample_num)
                if(verbose){message(" - found ", nrow(sub.conns), " potential co-accessible ACRs in shuffled group ",z, "...")}

                # get FDR
                fdr.conns <- getFDR(sub.conns, sub.shufconns, fdr=fdr_thresh, verbose=verbose)
                fdr.conns$group <- z

                # return
                return(fdr.conns)

            }, mc.cores=nthreads)

            # join into a single df
            outs <- as.data.frame(do.call(rbind, outs))
            obj[[conn_slotName]] <- outs

        }else{
            conns1 <- runCicero(cds, genome=genome, k=k, win=win_size, sample_num=sample_num)
            if(verbose){message(" - found ", nrow(conns1), " potential co-accessible ACRs in observed data ...")}
            conns2 <- runCicero(shufcds, genome=genome, k=k, win=win_size, sample_num=sample_num)
            if(verbose){message(" - found ", nrow(conns2), " potential co-accessible ACRs in expected data ...")}
            fdr.conns <- getFDR(conns1, conns2, fdr=fdr_thresh, verbose=verbose)
            obj[[conn_slotName]] <- fdr.conns
        }

    }

    # return
    return(obj)

}


###################################################################################################
###################################################################################################
###################################################################################################
#' runCicero
#'
#' This function is a wrapper for several cicero commands
#'
#' @import cicero
#'
#' @param cds cicero object. Required.
#' @param genome path to genome file. See 'inst/extdata/genome' for an example. Required.
#' @param k numeric. Number of nearest neighbors for binning cells. See cicero for more details.
#' Defaults to 30.
#' @param win numeric. Window size for constraining co-accessible ACRs. Defaults to 500000.
#' @param sample_num numeric. Number of random windows to sample for determining distance_parameter
#'  during 'run_cicero'. See cicero for more details.
#' @param silent logical. Defaults to TRUE.
#'
#' @rdname runCicero
#'
runCicero <- function(cds,
                      genome=NULL,
                      k=30,
                      win=500000,
                      sample_num=100,
                      silent=TRUE){

    # get UMAP coordinates
    umap_coords <- t(reducedDimA(cds))
    umap_coords <- umap_coords[colnames(cds),]
    rownames(umap_coords) <- colnames(exprs(cds))
    cicero_cds <- make_cicero_cds(cds, reduced_coordinates=umap_coords, k=k)

    # run cicero (connections)
    conns <- run_cicero(cicero_cds, genome, window=win, sample_num=sample_num, silent=silent)

    # return
    return(conns)

}


###################################################################################################
###################################################################################################
###################################################################################################
#' getFDR
#'
#' This function filters coACRs by empirical FDR
#'
#' @param observed observed coACRs.
#' @param expected shuffled coACRs.
#' @param fdr float. Defaults to 0.05.
#' @param grid grid size for determining thresholds. Defaults to 100.
#' @param verbose logical. Defaults to TRUE.
#'
#' @rdname getFDR
#'
getFDR <- function(obs,
                   exp,
                   fdr=0.05,
                   grid=100,
                   verbose=F){

    # split into +/-
    n.exp <- subset(exp, exp$x < 0)
    p.exp <- subset(exp, exp$x > 0)

    # get counts
    pos.nexp <- nrow(p.exp)
    neg.nexp <- nrow(n.exp)
    pos.nobs <- nrow(subset(obs, obs$x > 0))
    neg.nobs <- nrow(subset(obs, obs$x < 0))
    if(verbose){message(" - number expected links = (+) ",pos.nexp, " | (-) ",neg.nexp)}
    if(verbose){message(" - number observed links = (+) ",pos.nobs, " | (-) ",neg.nobs)}

    # generate range of thresholds
    p.vals <- seq(from=0, to=1, length.out=grid)
    n.vals <- seq(from=0, to= -1, length.out=grid)

    # iterate over grid
    if(verbose){message(" - scanning positive thresholds ...")}
    p.thresh <- c()
    for(i in p.vals){
        num.exp <- sum(p.exp$x > as.numeric(i))
        c.fdr <- num.exp/(pos.nexp)
        if(is.na(c.fdr)){
            c.fdr <- 0
        }
        p.thresh <- c(p.thresh, c.fdr)
        if(verbose){message(" - (+) correlation threshold = ", i, " | FDR = ", c.fdr)}
    }
    if(verbose){message(" - scanning negative thresholds ...")}
    n.thresh <- c()
    for(i in n.vals){
        num.exp <- sum(n.exp$x < as.numeric(i))
        c.fdr <- num.exp/(neg.nexp)
        if(is.na(c.fdr)){
            c.fdr <- 0
        }
        n.thresh <- c(n.thresh, c.fdr)
        if(verbose){message(" - (-) correlation threshold = ", i, " | FDR = ", c.fdr)}
    }

    # select cut-offs
    p.threshold <- min(p.vals[which(p.thresh <= fdr)])
    n.threshold <- max(n.vals[which(n.thresh <= fdr)])

    # filter
    obs <- subset(obs, obs$x > p.threshold | obs$x < n.threshold)

    # verbose number of +/- linkages
    pos.links <- nrow(subset(obs, obs$x > 0))
    neg.links <- nrow(subset(obs, obs$x < 0))
    if(verbose){message(" - found ",pos.links, " + and ", neg.links," - ACR-ACR links ...")}

    # return
    return(obs)
}

