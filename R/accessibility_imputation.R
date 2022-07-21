###################################################################################################
###################################################################################################
###################################################################################################
#' accessibility imputation
#' this function will impute gene accessibility from a raw tn5 count matrix
#'
#'
#'

# libraries
# load libraries
library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
library(varistran)
library(edgeR)
library(parallel)
library(png)


# tfidf run
safe_tfidf         <- function(tf, 
                               idf,  
                               block_size=2000e6){
    result = tryCatch({
        result = tf * idf
        result
    }, error = function(e) {
        options(DelayedArray.block.size=block_size)
        DelayedArray:::set_verbose_block_processing(TRUE)
        
        tf = DelayedArray(tf)
        idf = as.matrix(idf)
        
        result = tf * idf
        result
    })
    return(result)
}


tfidf              <- function(bmat, 
                               frequencies=F, 
                               log_scale_tf=T, 
                               scale_factor=100000){
    
    # Use either raw counts or divide by total counts in each cell
    if (frequencies) {
        # "term frequency" method
        tf = t(t(bmat) / Matrix::colSums(bmat))
    } else {
        # "raw count" method
        tf = bmat
    }
    
    # Either TF method can optionally be log scaled
    if (log_scale_tf) {
        if (frequencies) {
            tf@x = log1p(tf@x * scale_factor)
        } else {
            tf@x = log1p(tf@x * 1)
        }
    }
    
    # IDF w/ "inverse document frequency smooth" method
    idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
    
    # TF-IDF
    tf_idf_counts = safe_tfidf(tf, idf)
    rownames(tf_idf_counts) = rownames(bmat)
    colnames(tf_idf_counts) = colnames(bmat)
    return(Matrix(tf_idf_counts, sparse=T))
}


estimate_sf_sparse <- function(counts, 
                               round_exprs=T, 
                               method="mean-geometric-mean-total"){
    if (round_exprs)
        counts <- round(counts)

    if(method == 'mean-geometric-mean-total') {
        cell_total <- Matrix::colSums(counts)
        sfs <- cell_total / exp(mean(log(cell_total)))
    }else if(method == 'mean-geometric-mean-log-total') {
        cell_total <- Matrix::colSums(counts)
        sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
    }

    sfs[is.na(sfs)] <- 1
    sfs
}


normalizeSF        <- function(adj.act, 
                               verbose=F){
    
    # verbose
    if(verbose==T){
        message(" - unnormalized activities:")
        print(head(adj.act[,1:5]))
    }
    
    # get per sample norm factors
    norm.factor <- estimate_sf_sparse(adj.act)
    if(verbose==T){
        message(" - normalization factors:")
        print(head(norm.factor))
    }
    
    # normalized counts
    norm.act <- adj.act %*% Diagonal(x=1/norm.factor)
    colnames(norm.act) <- colnames(adj.act)
    rownames(norm.act) <- rownames(adj.act)
    if(verbose==T){
        message(" - normalized activities:")
        print(head(norm.act[,1:5]))
    }
    
    # output
    return(list(norm.act=norm.act, norm.factor=norm.factor))
}

normalize.activity <- function(df,
                               acts,
                               tar_clust_colnm,
                               n.random=NULL,
                               output="output",
                               logTransform=F,
                               scaleP=F,
                               plotHeatmapRaw=F){

    # verbose
    message(" - normalizing libraries ...")

    ##set a Cluster in df
    df$Cluster <- df[[tar_clust_colnm]]

    # quick clean
    acts <- acts[Matrix::rowSums(acts>0)>0,]
    acts <- acts[,Matrix::colSums(acts>0)>0]

    # check data frames/matrices
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,colnames(acts) %in% rownames(df)]
    df <- df[colnames(acts),]
    acts.o <- acts

    # if select random
    if(!is.null(n.random)){
        acts <- acts[sample(nrow(acts), n.random),]
    }

    # find cluster means
    its <- 0
    message(" - extracting average gene activity across cells per cluster ...")
    clust.means <- lapply(unique(df$Cluster), function(x){
        clust.ids <- rownames(subset(df, df$Cluster==x))
        Matrix::rowMeans(acts[,clust.ids])
    })
    clust.means <- do.call(cbind, clust.means)
    colnames(clust.means) <- unique(df$Cluster)

    # plot raw cluster means
    clust.means <- as.matrix(clust.means)
    c.means <- clust.means[,mixedorder(colnames(clust.means))]
    row.o <- apply(c.means, 1, which.max)
    c.means <- c.means[order(row.o, decreasing=F),]
    row.o <- rownames(c.means)

    # if do plot heatmap
    if(plotHeatmapRaw==TRUE){
        pdf(paste0(output,".raw_heatmap.pdf"), width=5, height=6)
        heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
                  useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
        dev.off()
    }
   # write.table(c.means, file=paste0(output_dir,'/',output,".RawClusterMeans.txt"),
   #             quote=F, row.names=T, col.names=T, sep="\t")


    # fit 2-guassian mixture model, return higher threshold
    message(" - fitting mixture model ...")
    thresholds <- apply(clust.means, 2, function(x){
        mod <- Mclust(x, G=2, verbose=F)
        top <- names(mod$parameters$mean)[which.max(mod$parameters$mean)]
        upper.cells <- x[which(mod$classification==top)]
        val <- quantile(upper.cells[upper.cells>0], c(0.05))
        #val <- mean(upper.cells, na.rm=T)
        if(val == 0){
            val <- min(upper.cells[upper.cells>0])
        }
        return(val)
        #max(mod$parameters$mean)
    })

    # scale cell activity by cluster-specific threshold for expression
    message(" - scaling by cluster averages ...")
    adj.thresh <- thresholds[df$Cluster]
    print(head(adj.thresh, n=10))
    adj.act <- acts.o %*% Diagonal(x=1/adj.thresh)
    adj.act@x[is.infinite(adj.act@x)] <- 0
    adj.act@x[is.na(adj.act@x)] <- 0
    adj.act@x <- round(adj.act@x, digits=0)
    adj.act <- adj.act[Matrix::rowSums(adj.act>0)>0,]
    adj.act <- adj.act[,Matrix::colSums(adj.act>0)>0]

    # un-normalized counts
    ua.out <- as.data.frame(summary(adj.act))
    ua.out$i <- rownames(adj.act)[ua.out$i]
    ua.out$j <- colnames(adj.act)[ua.out$j]
    ua.out <- ua.out[ua.out$x>0,]
    #write.table(ua.out, file=paste0(output_dir,'/',output,".countsActivity.sparse"),
    #            quote=F, row.names=F, col.names=F, sep="\t")

    # normalize by size factors
    if(logTransform==F){
        message(" - lib size before size factors = ")
        print(head(Matrix::colSums(adj.act)))
    }

    # do normalization
    message(" - estimating normalization factors ...")
    results <- normalizeSF(adj.act, verbose=T)
    norm.act <- results$norm.act
    norm.factor <- results$norm.factor

    # log transform?
    if(logTransform==T){
        message(" - square-root transforming counts activity ...")
        norm.act <- sqrt(norm.act)
        norm.act@x[is.na(norm.act@x)] <- 0
        norm.act@x[is.infinite(norm.act@x)] <- 0
        message(" - lib size afer square-root transformation = ")
        print(head(Matrix::colSums(norm.act)))
    }

    # write size factors to meta data file
    df$size_factors <- norm.factor[rownames(df)]
    #write.table(df, file=paste0(output_dir,'/',output,".size_factors.txt"),
    #            quote=F, row.names=T, col.names=T, sep="\t")

    # verbose
    message(" - lib size after size factors = ")
    print(head(Matrix::colSums(norm.act)))

    # remove empty cells/genes
    norm.act <- norm.act[Matrix::rowSums(norm.act>0)>0,]
    norm.act <- norm.act[,Matrix::colSums(norm.act>0)>0]
    message(" - cells = ",ncol(norm.act), " | genes = ", nrow(norm.act))
    print(head(norm.act[,1:10]))

    # print output to disk
    ia.out <- as.data.frame(summary(norm.act))
    ia.out$i <- rownames(adj.act)[ia.out$i]
    ia.out$j <- colnames(adj.act)[ia.out$j]
    ia.out <- ia.out[ia.out$x>0,]

    #write.table(ia.out, file=paste0(output_dir,'/',output,".normalizedActivity.sparse"), quote=F, row.names=F,
    #            col.names=F, sep="\t")

    #saveRDS(norm.act,paste0(output_dir,'/',output,".normalizedActivity_mtx.rds"))

    # return
    message(" - returning normalized matrix ...")
    return(list(norm.act=norm.act, norm.factor=norm.factor, adj.act=adj.act, row.o=row.o))
}


acc.imputation      <- function(x,
                               meta_slot = "Clusters",
                               tar_clust_colnm='LouvainClusters',
                               imputation_slot = "impute.acc",
                               k=25,
                               step=3,
                               npcs=30,
                               cleanExp=F,
                               df=NULL,
                               rds=NULL,
                               prefix="output"){

    message(' - set the ipt files')

    ##read the objective
    ipt_obj <- x
    all.b <- ipt_obj[[meta_slot]]
    rownames(all.b) <- all.b$cellID
    all.hpcs <- ipt_obj$SVD

    ##extract the gene activity data
    gene_acc_sparse <- ipt_obj$sc_gene_ac


    gene_acc_sparse <- as.data.frame(unclass(gene_acc_sparse),stringsAsFactors=TRUE)
    ##transfer the sparse to the matrix
    all.activity <- sparseMatrix(i=as.numeric(gene_acc_sparse$gene_name),
                                 j=as.numeric(gene_acc_sparse$barcode),
                                 x=as.numeric(gene_acc_sparse$accessability),
                                 dimnames=list(levels(gene_acc_sparse$gene_name), levels(gene_acc_sparse$barcode)))
    all.activity <- as(all.activity, "dgCMatrix")

    ##align all ids
    all.activity <- all.activity[,rownames(all.b)]
    all.activity <- all.activity[Matrix::rowSums(all.activity)>0,]
    all.activity <- all.activity[,Matrix::colSums(all.activity)>0]
    ids <- intersect(rownames(all.b), colnames(all.activity))
    ids <- intersect(ids, rownames(all.hpcs))
    all.b <- all.b[ids,]
    all.activity <- all.activity[,ids]
    all.hpcs <- all.hpcs[ids,]

    ##step01 normalize the activity
    message(paste0('-the target clust column name is ',tar_clust_colnm))
    results <- normalize.activity(all.b,
                                  all.activity,
                                  tar_clust_colnm,
                                  output=prefix,
                                  logTransform=F,
                                  scaleP=F,
                                  plotHeatmapRaw=F)
    activity.norm <- results$norm.act
    print(head(activity.norm))

    ##line-up all ids again
    barcode.ids <- intersect(colnames(activity.norm), rownames(all.hpcs))
    barcode.ids <- intersect(barcode.ids, rownames(all.b))
    activity.norm <- activity.norm[,barcode.ids]
    all.hpcs <- all.hpcs[barcode.ids,]
    all.b <- all.b[barcode.ids,]

    message(" - imputing gene activity ...")

    # input
    data.use <- activity.norm
    x <- activity.norm

    # verbose
    if(!is.null(rds)){

        if(!is.null(df)){
            message("   * using UMAP manifold for smoothing ...")
            pcs <- df[,c("umap1","umap2")]
        }else{
            message("   * using prior PC space as manifold ...")
            pcs <- rds[colnames(x),c(1:npcs)]
        }
    }else{

        # LSI
        message("   * PC manifold set to NULL, running LSI (TFIDF)...")
        x[x>0] <- 1
        tf.idf <- tfidf(x)

        # get PCS
        message("   * PC manifold set to NULL, running LSI ...")
        pc <- irlba(t(tf.idf), npcs)
        pcs <- pc$u
        rownames(pcs) <- colnames(x)
        colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))

        # do l2-norm
        pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
    }

    # get KNN
    message("   * finding knn graph ...")
    knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
    j <- as.numeric(x = t(x = knn.graph))
    i <- ((1:length(x = j)) - 1) %/% k + 1
    edgeList = data.frame(i, j, 1)
    A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])

    # Smooth graph
    ##some questions
    ##1. what's meanning of A
    ##2. A/Matrix::rowSums(A)
    ##3. meanning of A%*%A
    message("   * smoothing graph ...")
    A = A + t(A) ##A now is symmetric
    A = A / Matrix::rowSums(A)
    step.size = step
    if(step.size > 1){
        for(i in 1:step.size){
            message("     ~ step ",i)
            A = A %*% A
        }
    }

    # smooth data
    message("   * smoothing activity ...")
    out <- t(A %*% t(data.use))
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    print(head(out[,1:10]))

    # find non-expressed genes
    if(cleanExp==T){
        message("   * finding activity thresholds ...")
        its <- 0
        new <- sparse_apply(out, 2, function(z){
            its <<- its + 1
            mod <- Mclust(z, G=2)
            top.parameter <- which.max(mod$parameters$mean)
            vals <- split(z, mod$classification)[[top.parameter]]
            threshold <- quantile(vals[vals>0], c(0.05))[1]
            if(its %% 100 == 0){message("   * iterated over ",its, " cells (min thresh = ",threshold,") ...")}
            z[z < threshold] <- 0
            return(z)
        }, convert_to_dense=F)
        impute.activity <- Matrix(new, sparse=T)
    }else{
        impute.activity <- out
    }

    # clean empty rows (if present) and round to two decimal places
    #print(head(impute.activity[,1:10]))
    message("   * clean after imputation ...")
    #impute.activity <- impute.activity[Matrix::rowSums(impute.activity)>0,]
    #impute.activity <- impute.activity[,Matrix::colSums(impute.activity)>0]
    #impute.activity <- impute.activity %*% Diagonal(x=1e6/Matrix::colSums(impute.activity))
    #impute.activity@x <- round(impute.activity@x, digits=2)

    # write to disk
    # ia.out <- as.data.frame(summary(impute.activity))
    # ia.out$i <- rownames(impute.activity)[ia.out$i]
    # ia.out$j <- colnames(impute.activity)[ia.out$j]
    # write.table(ia.out, file=paste0(output,".imputedActivity.sparse"), quote=F, row.names=F,
    #             col.names=F, sep="\t")

    # return sparse Matrix

    ##add the impute.activity to the obj
    ##return a sparse
    ia.out <- as.data.frame(summary(impute.activity))
    ia.out$i <- rownames(impute.activity)[ia.out$i]
    ia.out$j <- colnames(impute.activity)[ia.out$j]
    # write.table(ia.out, file=paste0(output,".imputedActivity.sparse"), quote=F, row.names=F,
    #             col.names=F, sep="\t")

    ipt_obj[[imputation_slot]] <- impute.activity

    return(ipt_obj)
}














