###################################################################################################
###################################################################################################
###################################################################################################
#' Filter singleton cells from the UMAP embedding
#'
#' This function removes cells from the UMAP embedding that are not well supported by other cells.
#' Can be run on a per cluster basis to filter cells that exceed a heuristic threshold of distances
#' to most other cells in the same cluster.
#'
#' @importFrom FNN get.knn
#'
#' @param pro data.frame containing UMAP coordinates m1, and m2. See m1, m2 arguments for details.
#' @param k numeric, number of nearest neighbors. Defaults to 50.
#' @param threshold numeric, Z-score threshold to remove cells greater than X standard standard
#' deviations from other cells.
#' @param type character, which embedding to use. Must be one of "umap", "pca", or "tsne".
#' @param m1 character, if type=="umap", then set to the column name of first component ("umap1").
#' Only used when running sub-clustering (umap embedding names for subset and regenerated UMAP
#' embeddings of a subset of cells).
#' @param m2 character, if type=="umap", then set to the column name of second component ("umap2").
#' #' Only used when running sub-clustering (umap embedding names for subset and regenerated UMAP
#' embeddings of a subset of cells).
filterSingle  <- function(pro,
                          k=50,
                          threshold=3,
                          type="umap",
                          m1="umap1",
                          m2="umap2"){

    # set column names
    if(type=="umap"){
        vars <- c(m1, m2)
    }else if(type=="pca"){
        vars <- colnames(pro)
    }

    # ensure that k is less than number of samples
    if(k > nrow(pro)){
        k <- nrow(pro)-1
    }

    # get nearest neighbors
    topk <- FNN::get.knn(pro[,vars], k=k)
    cell.dists <- as.matrix(topk$nn.dist)
    rownames(cell.dists) <- rownames(pro)
    colnames(cell.dists) <- paste0("k",seq(1:ncol(cell.dists)))
    aves <- apply(cell.dists, 1, mean)
    zscore <- as.numeric(scale(aves))
    names(zscore) <- rownames(pro)

    # thresholds
    p.zscore <- zscore[order(zscore, decreasing=T)]
    num.pass <- length(zscore[zscore < threshold])

    # filter
    prop.good <- zscore[zscore < threshold]
    ids <- names(prop.good)
    out <- pro[rownames(pro) %in% ids,]

    # return object
    return(out)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Filter singleton cells from the UMAP embedding by distance per cluster
#'
#' This function removes cells from the UMAP embedding that are not well supported by other cells.
#' Can be run on a per cluster basis to filter cells that exceed a heuristic threshold of distances
#' to most other cells in the same cluster.
#'
#' @importFrom FNN get.knn
#'
#' @param b meta data after clustering with Seurat
#' @param m1 character, if type=="umap", then set to the column name of first component ("umap1").
#' Only used when running sub-clustering (umap embedding names for subset and regenerated UMAP
#' embeddings of a subset of cells).
#' @param m2 character, if type=="umap", then set to the column name of second component ("umap2").
#' #' Only used when running sub-clustering (umap embedding names for subset and regenerated UMAP
#' embeddings of a subset of cells).
#' @param threshold2 numeric, Z-score threshold to remove cells greater than X standard standard
#' deviations from other cells.
filtDistClst <- function(b,
                         umap1="umap1",
                         umap2="umap2",
                         threshold2=2){

    # iterate over each cluster
    clusts <- unique(b$seurat_clusters)
    out <- lapply(clusts, function(x){
        b.sub <- subset(b, b$seurat_clusters == x)
        b.umap <- b.sub[,c(umap1, umap2)]
        out.b.umap <- filterSingle(b.umap, k=25, threshold=threshold2, m1=umap1, m2=umap2)
        return(rownames(out.b.umap))
    })
    out <- do.call(c, out)
    b.out <- b[rownames(b) %in% as.character(out),]
    message("   * total number of cells surviving subcluster filtering = ", nrow(b.out))
    return(b.out)

}


###################################################################################################
###################################################################################################
###################################################################################################
#' Cluster cells using Louvain, Leiden or other graph-based methods implemented by Seurat
#'
#' @import Seurat
#' @import SeuratObject
#' @import Matrix
#'
#' @param obj list, object containing 'PCA', 'UMAP', 'counts' and 'meta'. Required.
#' @param res numeric, resolution for Louvain/Leiden clustering. Typical values range from 0 to 2.
#' Defaults to 0.8.
#' @param k.near numeric, number of nearest neighbors for graph building. Defaults to 20.
#' @param clustOB character, which embedding to use for clustering. Can be one of c("umap", "svd").
#' It is highly recommended to use 'svd' for graph-based clustering. Defaults to "svd".
#' @param cname character, column name in meta cluster IDs. Defaults to "LouvainClusters".
#' @param min.reads numeric, minimum number of aggregated nSites (number of accessible peaks per
#' cell) to consider a cluster as valid. Note: that this does not consider the total number of reads
#' in each cell. Rather, we take the sum of column sums (colSum) for cells within a single cluster.
#' Defaults to 5e4.
#' @param m.clst numeric, minimum number of cells to consider a cluster as valid. Defaults to 25.
#' @param threshold numeric, Z-score threshold to filter out cells greater than X distance to other
#' cells on average. Uses a knn to find the top 50 nearest cells. See 'filterSingle' for more
#' details.
#' @param umap1 character, column name of UMAP first component. Only important for sub-clustering.
#' @param umap2 character, column name of UMAP second component. Only important for sub-clustering.
#' @param key character, name of UMAP key when creating Seurat object. Defaults to "UMAP_".
#' @param e.thresh numeric, Z-score threshold to filter out cells failing to co-localize with other
#' cells apart of the same cluster. Set to 1e6 to skip cluster refining. This step may also be
#' skipped by setting cleanCluster to FALSE. Default is 2.
#' @param cleanCluster logical, whether or not to implement cluster refinement by removing cells
#' that do not co-localize the the majority of cells apart of the same cluster. Uses the UMAP
#' embedding for distance estimates.
#' @param cl.method numeric, graph clustering algorithm. Numeric values range from 1-4, see
#' Seurat::FindClusters for more details.
#' @param svd_slotName character, character string for the SVD slot to use for clustering. Defaults
#' to "PCA".
#' @param umap_slotName character, character string for the UMAP slot to use for analysis. Defaults
#' to "UMAP.
#' @param cluster_slotName character, character string for naming the results of callClusters.
#' Defaults to "Clusters".
#' @param verbose logical. Defaults to FALSE.
#' @param ... additional arguments sent to Seurat::FindClusters
#'
#' @rdname callClusters
#' @export
callClusters  <- function(obj,
                          res=0.4,
                          k.near=30,
                          clustOB="svd",
                          cname="LouvainClusters",
                          min.reads=5e4,
                          m.clst=25,
                          threshold=5,
                          umap1="umap1",
                          umap2="umap2",
                          e.thresh=3,
                          cleanCluster=T,
                          cl.method=1,
                          svd_slotName="PCA",
                          umap_slotName="UMAP",
                          cluster_slotName="Clusters",
                          verbose=FALSE,
                          ...){

    # filter umap coordinates
    if(verbose){message(" - filtering outliers in UMAP manifold (z-score e.thresh = ", e.thresh, ") ...")}
    umap.original <- obj[[umap_slotName]]
    umap.filtered <- filterSingle(obj[[umap_slotName]], threshold=e.thresh, m1=umap1, m2=umap2)
    counts.filtered <- obj$counts[,rownames(umap.filtered)]
    meta.filtered <- obj$meta[rownames(umap.filtered),]
    pca.filtered <-  obj[[svd_slotName]][rownames(umap.filtered),]

    # run graph-based clustering
    if(verbose){message(" - creating seurat object for graph-based clustering ...")}

    # create Seurat object, find clusters
    sro <- SeuratObject::CreateSeuratObject(counts.filtered, min.cells=0, min.features=0)
    sro[["svd"]] <- CreateDimReducObject(embeddings = pca.filtered, key = "PC_", assay = DefaultAssay(sro))
    sro[["umap"]] <- CreateDimReducObject(embeddings=as.matrix(umap.filtered), key="UMAP_", assay=DefaultAssay(sro))
    sro <- AddMetaData(sro, meta.filtered)
    if(nrow(meta.filtered) > 250000){
        nn.eps.val <- 0.25
        n.starts <- 10
    }else{
        nn.eps.val <- 0
        n.starts <- 100
    }
    sro <- FindNeighbors(sro, dims = 1:ncol(sro[[clustOB]]), reduction=clustOB, 
                         nn.eps=nn.eps.val, k.param=k.near, annoy.metric="cosine")
    sro <- FindClusters(sro, resolution=res, n.start=n.starts, algorithm=cl.method, ...)
    sro.meta <- data.frame(sro@meta.data)
    sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
    if(verbose){message(" - finished graph-based clustering ...")}
    
    # remove temp
    rm(umap.filtered)
    rm(counts.filtered)
    rm(meta.filtered)
    rm(pca.filtered)
    suppressMessages(gc())

    # remove outliers?
    if(cleanCluster){

        # verbose
        if(verbose){message("   * removing low quality clusters ...")}

        # prep
        sro.umap <- obj[[umap_slotName]][rownames(sro.meta),]
        colnames(sro.umap) <- c(umap1, umap2)
        sro.meta <- cbind(sro.meta, sro.umap)

        # filter by cluster size and number of cells
        agg.reads <- aggregate(sro.meta$nSites~sro.meta$seurat_clusters, FUN=sum)
        colnames(agg.reads) <- c("clusters","readDepth")
        clust.cnts <- table(sro.meta$seurat_clusters)
        agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
        agg.pass <- subset(agg.reads, agg.reads$num_cells >= m.clst & agg.reads$readDepth >= min.reads)
        sro.filt <- sro.meta[as.character(sro.meta$seurat_clusters) %in% as.character(agg.pass$clusters),]
        sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
        sro.filt$seurat_clusters <- as.numeric(factor(sro.filt$seurat_clusters))

        # remove outliers in the embedding
        if(verbose){message("   * filtering per-cluster outliers (z-score filtDistClst2 = ", threshold, ") ...")}
        sro.meta <- filtDistClst(sro.filt, umap1=umap1, umap2=umap2, threshold=threshold)
        sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)

    }

    # filter by cluster size
    if(verbose){message(" - filtering clusters with low cell/read counts ...")}
    agg.reads <- aggregate(sro.meta$nSites~sro.meta$seurat_clusters, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(sro.meta$seurat_clusters)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells>=m.clst & agg.reads$readDepth>=min.reads)
    sro.filt <- sro.meta[sro.meta$seurat_clusters %in% agg.pass$clusters,]
    sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
    sro.filt$seurat_clusters <- as.numeric(as.factor(sro.filt$seurat_clusters))

    # rename output
    final.UMAP <- umap.original[rownames(sro.filt),]
    clusts <- sro.filt$seurat_clusters
    obj[[cluster_slotName]] <- obj$meta[rownames(sro.filt),]
    obj[[cluster_slotName]][,cname] <- factor(clusts)
    obj[[cluster_slotName]][,c(umap1)] <- final.UMAP[,c(umap1)]
    obj[[cluster_slotName]][,c(umap2)] <- final.UMAP[,c(umap2)]

    # return
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
# reclust       <- function(a, b, variates, r.variates, modtype="regModel", chunks=1000,
#                           type="response", subpeaks=1000, ncpus=1, doPlot=F, prefix="out",
#                           stdize=T, doCenter=F, link="logit", method="elasticNet", alpha=0.5,
#                           var.th=1, scaleVar=T, doL2=1, cap=F, useTop=NULL, n.pcs=50,
#                           batch=NULL, theta=1, lambda=1, useTop.th=0.75, min.reads=5e5,
#                           m.clst=25, threshold= -0.25, k.nearC=15, res=0.03, clustOB="umap",
#                           key="UMAP", ...){
#
#     # get clusters to iterate over
#     b$LouvainClusters <- factor(b$LouvainClusters)
#     LC <- mixedsort(unique(b$LouvainClusters))
#
#     # iterate
#     finalresults <- mclapply(LC, function(i){
#
#         # update
#         prefix1 <- paste0(prefix,"_",i)
#         clustname <- "LouvainCluster_sub"
#
#         # verbose
#         message(" ------------ ", i, " ------------")
#
#         # clust
#         clust <- rownames(subset(b, b$LouvainClusters==i))
#         aa <- a[,clust]
#         bb <- b[colnames(aa),]
#         aa <- aa[Matrix::rowSums(aa)>=(ncol(aa)*0.005),]
#         aa <- aa[,Matrix::colSums(aa)>50]
#         aa <- aa[Matrix::rowSums(aa)>0,]
#         aa <- aa[,Matrix::colSums(aa)>0]
#
#         if(nrow(bb) < 200 | sum(bb$unique) < 2e6){
#             message(" - skipping Louvain cluster ", i, ": too few cells or reads")
#             bb$umapsub_1 <- bb$umap1
#             bb$umapsub_2 <- bb$umap2
#             bb[,c(clustname)] <- factor(bb$LouvainClusters)
#             return(bb)
#         }
#
#         # rerun pipeline
#         dev <- runModel(aa, bb, variates, r.variates,
#                         subpeaks = subpeaks,
#                         chunks = 256,
#                         type = type,
#                         modtype = modtype,
#                         nthreads = ncpus,
#                         prefix = prefix1,
#                         doPlot = doPlot,
#                         stdize = stdize,
#                         link = link,
#                         doCenter = doCenter,
#                         method = method,
#                         alpha = alpha,
#                         var.th = var.th)
#
#         # std dev
#         dev <- stddev(dev)
#
#         # reduce dims
#         out.pcs  <- reduceDims(dev, bb,
#                                n.pcs = n.pcs,
#                                batch = batch,
#                                theta = theta,
#                                lambda = lambda,
#                                doL2 = doL2,
#                                cap = cap,
#                                prefix = prefix1,
#                                scaleVar = scaleVar,
#                                raw = aa,
#                                center = doCenter,
#                                useTop = useTop,
#                                useTop.th = useTop.th)
#
#         # project with UMAP
#         out.umap <- projectUMAP(out.pcs)
#         out.umap$umapsub_1 <- out.umap$umap1
#         out.umap$umapsub_2 <- out.umap$umap2
#         out.umap$umap1 <- NULL
#         out.umap$umap2 <- NULL
#
#         # call clusters
#         bb <- callClusters(aa, bb, out.pcs, out.umap,
#                            umap1="umapsub_1", umap2="umapsub_2",
#                            min.reads = min.reads,
#                            m.clst = m.clst,
#                            threshold = threshold,
#                            k.near = k.nearC,
#                            res = res,
#                            prefix = prefix1,
#                            clustOB = clustOB,
#                            cname = clustname,
#                            key = "umapsub_")
#
#         # make column factor
#         bb[,c(clustname)] <- factor(bb[,c(clustname)])
#         bb$sub_louvain <- paste(bb$LouvainClusters, bb[,c(clustname)], sep="_")
#
#         # write stats to shell
#         reportResults(bb, column=clustname)
#
#         # output
#         outData(out.pcs, bb, prefix=prefix1, dev=NULL)
#
#         # plot UMAP
#         plotUMAP(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1, column=clustname)
#         plotSTATS2(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1)
#
#         # update list
#         finalresults[[i]] <- bb
#
#         # check number of sub-clusters
#         if(length(unique(bb$sub_louvain))>1){
#
#             # load marker/gene activity
#             bb$umap1 <- bb$umapsub_1
#             bb$umap2 <- bb$umapsub_2
#             datt <- loadMData(bb, out.pcs, geneact, mark, clustID="sub_louvain")
#             bdat <- datt$b
#             activity <- datt$activity
#             h.pcs <- datt$h.pcs
#             marker.info <- datt$marker.info
#             print(head(bdat))
#             print(head(activity[,1:5]))
#             print(head(h.pcs))
#
#             # normalize per cell activity by cluster average and size factors
#             results <- normalize.activity(bdat, activity, output=prefix1, logTransform=F, scaleP=F)
#             activity <- results$norm.act
#             row.o <- results$row.o
#
#             # impute gene accessibility scores
#             impute.activity <- smooth.data(activity, k=15, step=3, npcs=s.pcs, df=NULL,
#                                            rds=h.pcs, cleanExp=F, output=prefix1)
#
#             # collect marker accessibility
#             plot.act.scores(bdat, acts=activity,
#                             info=marker.info,
#                             logT=T,
#                             outname=paste0(prefix1,".normalized.known.Markers.pdf"))
#
#             plot.act.scores(bdat, acts=impute.activity,
#                             info=marker.info,
#                             logT=F,
#                             outname=paste0(prefix1,".impute.known.Markers.pdf"))
#
#         }
#     }, mc.cores=nthreads)
#
#     # condense list
#     df <- do.call(rbind, finalresults)
#     rownames(df) <- df$cellID
#     df$sub_louvain <- paste(df$LouvainClusters, df$LouvainCluster_sub, sep="_")
#     df$tissue_cluster <- as.numeric(factor(df$sub_louvain, levels=mixedsort(unique(as.character(df$sub_louvain)))))
#     return(df)
#
# }


###################################################################################################
###################################################################################################
###################################################################################################
# reclust2      <- function(a, b, pcs, prefix="out", n.pcs=50, batch=NULL, theta=1, nthreads=1,
#                           lambda=1, min.reads=5e5, m.clst=25, threshold3=0.25, doL2=F,
#                           clustType="densityClust", k.nearC=15, res=0.03, clustOB="umap",
#                           cleanCluster=T, e.thresh=2, doSTDize=F){
#
#     # get clusters to iterate over
#     LC <- mixedsort(unique(as.character(b$LouvainClusters)))
#
#     # load data for marker analysis
#     mark <- "/scratch/apm25309/single_cell/ATACseq/v3/step1_clustering/model_based/markers.bed"
#     geneact <- "/scratch/apm25309/single_cell/ATACseq/v3/sparse/genes/gene_activity/all.ex.GBaccessibility.sparse"
#
#     # resolutions
#     #res <- rep(res, length(LC))
#     message(" - initializing sub-clustering ... ")
#     print(LC)
#
#     # adjust number of threads
#     if(length(LC) < nthreads){
#         nthreads <- length(LC)
#     }
#
#     # iterate
#     it <- 0
#     finalresults <- mclapply(LC, function(i){
#
#         # update
#         prefix1 <- paste0(prefix,"_",i)
#         clustname <- "LouvainCluster_sub"
#
#         # verbose
#         message("---------------------------------------------")
#         message("------------------ ", i, " ------------------")
#         message("---------------------------------------------")
#
#         # clust
#         clust <- rownames(subset(b, as.character(b$LouvainClusters)==i))
#         aa <- a[,colnames(a) %in% clust]
#         bb <- b[colnames(aa),]
#
#         # use minimum among the number of available PCs and requested components
#         if(ncol(pcs) < n.pcs){
#             n.pcs <- ncol(pcs)
#         }
#         out.pcs1 <- pcs[colnames(aa),1:n.pcs]
#
#         # standarize embeddings
#         if(doSTDize == T){
#             out.pcs1 <- t(apply(out.pcs1, 1, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T)))
#         }
#
#         # L2 norm embeddings
#         if(doL2 == T){
#             out.pcs1 <- t(apply(out.pcs1, 1, function(x) x/sqrt(sum(x^2))))
#         }
#
#         # verbose
#         message(" - loaded ", nrow(out.pcs1), " cells ...")
#
#         # do harmony batch correciton
#         if(!is.null(batch)){
#             out.pcs1 <- HarmonyMatrix(out.pcs1, bb, do_pca=F, vars_use=batch,
#                                       theta=2, lambda=1, nclust=100, tau=0,
#                                       sigma=0.1, block.size=0.01)
#         }
#
#         # project with UMAP
#         umap.met <- "euclidean"
#         out.umap <- as.data.frame(projectUMAP(out.pcs1, m.dist=0.1, k.near=k.nearC, metric=umap.met))
#         colnames(out.umap) <- c("umap1","umap2")
#         message(" - finished UMAP ...")
#         out.umap$umapsub_1 <- out.umap$umap1
#         out.umap$umapsub_2 <- out.umap$umap2
#         out.umap$umap1 <- NULL
#         out.umap$umap2 <- NULL
#
#         # call clusters
#         message(" - begin clustering ...")
#         bb <- callClusters(aa, bb, out.pcs1, out.umap,
#                            umap1="umapsub_1", umap2="umapsub_2",
#                            min.reads = min.reads,
#                            m.clst = m.clst,
#                            threshold3 = threshold3,
#                            k.near = k.nearC,
#                            res = res,
#                            prefix = prefix1,
#                            clustOB = clustOB,
#                            cname = clustname,
#                            clustType = clustType,
#                            dynamic = F,
#                            e.thresh = e.thresh,
#                            cleanCluster = cleanCluster,
#                            cl.method = 2)
#
#         # make column factor
#         bb$sub_louvain <- factor(paste(bb$LouvainClusters, bb[,c(clustname)], sep="_"))
#
#         # write stats to shell
#         reportResults(bb, column=clustname)
#
#         # plot UMAP
#         bb$LouvainClusters <- factor(bb$LouvainClusters)
#         bb[,clustname] <- factor(bb[,clustname])
#         plotUMAP(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1, column=clustname)
#         plotSTATS2(bb, m1="umapsub_1", m2="umapsub_2", prefix=prefix1)
#
#         # update list
#         dd <- bb
#         message(" - clusters contain ", nrow(bb), " cells ... (out of = ",nrow(out.pcs1),")")
#
#         # check number of sub-clusters
#         if(length(unique(bb$sub_louvain))>1){
#
#             # load marker/gene activity
#             bb$umap1 <- bb$umapsub_1
#             bb$umap2 <- bb$umapsub_2
#             datt <- loadMData(bb, out.pcs1, geneact, mark, clustID="sub_louvain")
#             bdat <- datt$b
#             activity1 <- datt$activity
#             h.pcs1 <- datt$h.pcs
#             marker.info1 <- datt$marker.info
#             print(head(bdat))
#             print(head(activity1[,1:5]))
#             print(head(h.pcs1))
#
#             # normalize per cell activity by cluster average and size factors
#             results <- normalize.activity(bdat, activity1, output=prefix1, logTransform=F, scaleP=F)
#             activity1 <- results$norm.act
#             row.o <- results$row.o
#
#             # impute gene accessibility scores
#             impute.activity1 <- smooth.data(activity1, k=15, step=3, npcs=ncol(h.pcs1), df=NULL,
#                                             rds=h.pcs1, cleanExp=F, output=prefix1)
#
#             # collect marker accessibility
#             plot.act.scores(bdat, acts=activity1,
#                             info=marker.info1,
#                             logT=T,
#                             outname=paste0(prefix1,".normalized.known.Markers.pdf"))
#
#             plot.act.scores(bdat, acts=impute.activity1,
#                             info=marker.info1,
#                             logT=F,
#                             outname=paste0(prefix1,".impute.known.Markers.pdf"))
#
#         }
#
#         # return results
#         return(dd)
#
#     }, mc.cores=nthreads)
#
#     # condense list
#     message(" - combining child processes ... ")
#     df <- as.data.frame(do.call(rbind, finalresults))
#     df <- df[!duplicated(df$cellID),]
#     df <- df[!is.na(df$cellID),]
#     rownames(df) <- df$cellID
#     df$tissue_cluster <- as.numeric(factor(df$sub_louvain, levels=mixedsort(unique(as.character(df$sub_louvain)))))
#     tissue.props <- prop.table(table(df$tissue_cluster,df$tissue),1)
#     tissue.calls <- apply(tissue.props, 1, function(z) names(z)[which.max(z)] )
#     df$top_tissue_cluster <- paste(tissue.calls[df$tissue_cluster], names(tissue.calls[df$tissue_cluster]), sep="_")
#     print(warnings())
#     return(df)
#
# }
