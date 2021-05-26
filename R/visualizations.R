###################################################################################################
###################################################################################################
###################################################################################################
#' Plot UMAP embeddings
#'
#' Function to plot UMAP embeddings using meta data to establish cell colors.
#'
#' @import RColorBrewer
#' @import viridis
#' @import scales
#'
#' @param obj list, object containing meta data.
#' @param column character, column header specifying how to color cells. Factors are plotted
#' catagorically, while continuous numeric values are plotted along a spectrum.
#' @param cex float, set the point size. Defaults to 0.3.
#' @param opaque float, set the transparency of point colors. Defaults to 1.
#' @param cluster_slotName character, string specifying the desired UMAP slot to use for plotting.
#' Deafults to "UMAP".
#' @param xlab character string for x-axis name.
#' @param ylab character string for y-axis name.
#' @param main character string for graph title.
#' @param ... other arguments accepted by 'plot'.
#'
#' @rdname plotUMAP
#' @export
plotUMAP <- function(obj,
                     column="LouvainClusters",
                     cex=0.3,
		             opaque=1,
                     cluster_slotName="Clusters",
                     xlab="umap1",
                     ylab="umap2",
                     main="",
                     ...){

    # set b as meta data
    if(is.null(obj[[cluster_slotName]])){
        stop(" - ERROR: final.meta slot from callClusters is missing from object ...")
    }
    b <- obj[[cluster_slotName]]

    # test if column is present
    if(!column %in% colnames(b)){
        stop(" - ERROR: column header, ", column, ", is missing from meta ...")
    }

    # cols
    if(is.factor(b[,c(column)])){
        b <- b[sample(nrow(b)),]
        cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
        colv <- cols[factor(b[,column])]
    }else if(is.character(b[,column])){
        b[,column] <- factor(b[,column])
        b <- b[sample(nrow(b)),]
        cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
        colv <- cols[factor(b[,column])]
    }else if(is.numeric(b[,column])){
        b <- b[order(b[,column], decreasing=F),]
        cols <- viridis(100)
        colv <- cols[cut(b[,column], breaks=101)]
    }

    # plot
    plot(b$umap1, b$umap2, pch=16, cex=cex, col=alpha(colv,opaque),
         xlab=xlab,
         ylab=ylab,
         main=main,
         xlim=c(min(b$umap1), max(b$umap1)+(abs(max(b$umap1))*0.5)),
         ...)

    if(is.factor(b[,column])){
        legend("right", legend=sort(unique(b[,column])),
               fill=cols[sort(unique(b[,column]))])
    }

}


###################################################################################################
###################################################################################################
###################################################################################################
plotUMAPstats <- function(x, column="LouvainClusters", palette="Paired", m1="umap1", m2="umap2"){

    # require viridis
    require(viridis)

    # update coordinates
    x$umap11 <- as.numeric(x[,c(m1)])
    x$umap22 <- as.numeric(x[,c(m2)])
    x <- x[,c("umap11","umap22",column)]
    x <- x[complete.cases(x),]

    # set up color palette
    if(is.factor(x[,column])){
        message("   * plotting data in column: ",column," | class: ",class(x[,column]))
        colc <- colorRampPalette(brewer.pal(12,palette)[1:8])(length(unique(x[,column])))
        colv <- colc[factor(x[,column])]
        ids.names <- levels(factor(x[,column]))
        color.names <- colc[factor(ids.names, levels=ids.names)]
        type <- "fac"
    }else if(is.character(x[,column])){
        message("   * plotting data in column: ",column," | class: ",class(x[,column]))
        colc <- colorRampPalette(brewer.pal(12,palette)[1:8])(length(unique(x[,column])))
        colv <- colc[factor(x[,column])]
        ids.names <- levels(factor(x[,column]))
        color.names <- colc[factor(ids.names, levels=ids.names)]
        type <- "fac"
    }else if(is.numeric(x[,column]) | is.integer(x[,column])){
        message("   * plotting data in column: ",column," | class: ",class(x[,column]))
        numeric.nums <- as.numeric(x[,column])
        message("   *  converted to numeric ...")
        colc <- colorRampPalette(c("grey75","darkorchid4"))(100)
        message("   * successfully specified viridis color palette ...")
        colv <- colc[cut(numeric.nums, breaks=101)]
        type <- "num"
    }else{
        message("   * class for selected column is unsupported: ", class(x[,column]))
    }

    # plot
    plot(x$umap11, x$umap22, pch=16, cex=0.2, col=colv,
         bty="n", xaxt='n', yaxt='n',
         xlab="UMAP1", ylab="UMAP2", main=column,
         xlim=c(min(x$umap11),max(x$umap11)+abs(max(x$umap11)*0.5)))

    # legend
    if(type=="fac"){
        legend("right", legend=ids.names, fill=color.names, border=NA, cex=0.5)
    }else if(type=="num"){
        cc <- x[,column]
        half <- signif((min(cc)+max(cc))/2, digits=2)
        minn <- signif(min(cc), digits=2)
        maxx <- signif(max(cc), digits=2)
        legend("right", legend=c(minn,half,maxx), fill=c(colc[1],colc[51],colc[100]), border=NA)
    }

    # axes
    axis(1, lwd.tick=0, labels=FALSE)
    axis(2, lwd.tick=0, labels=FALSE)

}
