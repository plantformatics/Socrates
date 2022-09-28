###################################################################################################
###################################################################################################
###################################################################################################
#' loadBEDandGenomeData
#'
#' This function loads obj 5 column BED file specifying the single-bp Tn5 integration sites.
#' The file should not have any headers, with the following format:
#' chromosome, start, end, barcode, strand
#'
#' The annotation file should be GFF and the chromosome sizes file should be the chromsome ID
#' followed by its total length (e.g. chr1    123523000) with one chromosome per line.
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF
#'
#' @param bed Path to Tn5 integration site bed file with barcode information.
#' Required.
#' @param ann Path to GFF file for estimating % reads mapping to TSSs. Required.
#' @param sizes Path to chromosome sizes file for downstream processing. Required.
#' @param verbose logical. Defaults to TRUE.
#' @param is.fragment logical. Set to TRUE if the provided BED file is the fragment.csv output from
#' cellranger-atac count. Defaults to FALSE. 
#'
#' @rdname loadBEDandGenomeData
#' @export
#'
loadBEDandGenomeData <- function(bed, ann, sizes, attribute="Parent", verbose=T, is.fragment=F){

    # load hidden pre-check function
    .preRunChecks <- function(bed, ann, sizes, verbose=T){

        # verbose
        if(verbose){message("Running pre-check on input files and executable paths ...")}

        # check that files exist
        if(!file.exists(bed)){
            message("BED file, ", bed," does not exist... exiting")
            quit(save="no")
        }else{
            if(verbose){message("BED file path = ", bed, " ... ok")}
        }
        if(!file.exists(ann)){
            message("GFF annotation file, ", ann," does not exist... exiting")
            quit(save="no")
        }else{
            if(verbose){message("GFF file path = ", ann, " ... ok")}
        }
        if(!file.exists(sizes)){
            message("Chromosome sizes file, ", sizes," does not exist... exiting")
            quit(save="no")
        }else{
            if(verbose){message("Chromosome sizes file path = ", sizes, " ... ok")}
        }

        # check for macs2
        prg <- "macs2"
        if(Sys.which(prg) == ""){
            stop(message("Please install ", prg, ". ACR calling will not work without MACS2 in your path."))
        }else{
            if(verbose){message("Macs2 is installed .... ok")}
        }


    }
    .convertFragment2BED <- function(csv){
        
        
        
    }
    
    # run pre-checks
    .preRunChecks(bed, ann, sizes, verbose=verbose)

    # save paths
    bedpath <- bed
    annpath <- ann
    chrpath <- sizes

    # load args
    if(verbose){message(" - loading data (this may take obj while for big BED files) ...")}
    if(grepl(".gz$", bed)){
        obj <- read.table(gzfile(as.character(bed)))
    }else{
        obj <- read.table(as.character(bed))
    }
    if(grepl(".gtf", ann)){
        anntype <- "gtf"
    }else{
        anntype <- "gff3"
    }

    # check if input bed is obj fragment file
    if(is.fragment){
        
        if(verbose){message(" - converting fragment file to single-bp resolution Tn5 insertions sites ...")}
        start.coordinates <- data.frame(V1=obj$V1, V2=(obj$V2), V3=(obj$V2+1), V4=obj$V4, V5="+")
        end.coordinates <- data.frame(V1=obj$V1, V2=(obj$V3-1), V3=obj$V3, V4=obj$V4, V5="-")
        all.coordinates <- rbind(start.coordinates, end.coordinates)
        all.coordinates <- all.coordinates[order(all.coordinates$V1, all.coordinates$V2, decreasing=F),]
        obj <- all.coordinates[!duplicated(all.coordinates),]
        
    }
    
    # load Gff
    gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format=anntype, dbxrefTag=attribute)))
    chrom <- read.table(as.character(sizes))

    if(verbose){message(" - finished loading data")}
    return(list(bed=obj, gff=gff, chr=chrom, bedpath=bedpath, annpath=annpath, chrpath=chrpath))

}



###################################################################################################
###################################################################################################
###################################################################################################
#' countRemoveOrganelle
#'
#' This functions takes in an the Tn5 bed file, as well as obj vector
#' which contains organelle scaffolds, and identifies Tn5 integrations events which occured
#' in organelles. These reads are then assigned to the object slot PtMT. Additional parameters
#' allow for the removal of these reads so they don't interfere with ACR calling downstream.
#'
#' @param obj Object output from loadBEDandGenomeData. Required.
#' @param org_scaffolds vector of organelle scaffold names
#' @param remove_reads Logical of whether to remove reads from the BED file or not. Defaults to TRUE
#'
#' @rdname countRemoveOrganelle
#'
#' @export
#'
countRemoveOrganelle <- function(obj,
                                 org_scaffolds=NULL,
                                 remove_reads=FALSE){
    if (missing(org_scaffolds)) {

        message("... No Organelles Given, Returning All reads")

    } else {

        `%notin%` <- Negate(`%in%`)

        #Declare organelle scaffolds
        organelle = org_scaffolds

        #ID regions within the organelle
        org_group <- (obj$bed$V1 %in% organelle)

        if (sum(org_group) == 0) {
            message("No organeller reads identified ...")
            message("Are you sure the given names were correct?")
        }else{
            message("Identified ", sum(org_group), " organeller reads ...")
        }

        #Subset bed to only those
        organelle_sites <-  obj$bed[org_group, ]

        #Count the number of organelle reads present per barcode
        count_ID_number <- table(orgaelle_sites$V4)

        #Gather All Names to ID read names missing
        take_all_names <- unique(obj$bed$V4)

        #Grab all names with zero values
        cells_no_organelle_reads <- take_all_names[!(take_all_names %in% rownames(count_ID_number))]

        #Generate array using Table
        no_organelle_reads <- table(cells_no_organelle_reads)

        #Re-assing all values to 1
        no_organelle_reads[no_organelle_reads == 1] <- 0

        #Combine all
        zz <- c(no_organelle_reads, count_ID_number)

        if (remove_reads == TRUE){

            #Subsample bed so organelle reads do not interfere with ACR calls
            keep_bed_group <- c(obj$bed$V1 %notin% organelle)
            final_bed <- subset(obj$bed, keep_bed_group)
            obj$bed <- final_bed
            obj$PtMt <- zz

        } else {

            message("... Keeping organelle reads. This may affect ACR calling...")
            obj$PtMt <- zz

        }

        #Assign
        obj$PtMt <- zz
    }

    return(obj)

    }



###################################################################################################
###################################################################################################
###################################################################################################
#' callACRs
#'
#' This function calls ACRs using the bulk Tn5 integration sites with MACS2
#'
#' @param obj Object output from loadBEDandGenomeData. Required.
#' @param genomesize Effective genome size parameter for MACS2 (defaults to 1.6e9 for Zea mays).
#' @param shift Shift the start coordinates upstream. Defaults to -50 bp.
#' @param extsize Extends the fragment length downstream. Defaults to 100 bp.
#' @param output Output file names for MACS2. Defaults to 'bulk_peaks'.
#' @param tempdir Directory name to store the output of MACS2. Defaults to "./macs2_temp".
#' @param verbose Logical. Whether to print information.
#' @param fdr Float. Significance cut-off for peak calling. Defaults to q-value 0.05.
#'
#' @rdname callACRs
#' @export
#'
callACRs <- function(obj, genomesize=1.6e9, shift= -50, extsize=100,
                     output="bulk_peaks", tempdir="./macs2_temp", verbose=T, fdr=0.05){

    # verbose
    if(verbose){message(" - running MACS2 on bulk BED file ...")}
    bed <- obj$bedpath

    # create temp dif
    mac2temp <- tempdir
    if(file.exists(mac2temp)){
        unlink(mac2temp)
        dir.create(mac2temp)
    }else{
        dir.create(mac2temp)
    }

    # build command
    cmdline <- paste0("macs2 callpeak -t ", bed, " -f BED -g ", genomesize, " --keep-dup all -n ", output,
                      " --nomodel --shift ",shift, " --extsize ", extsize, " --outdir ",mac2temp, " --qvalue ", fdr)

    # run macs2
    suppressMessages(system(cmdline))

    # load peaks
    peaks <- read.table(paste0(mac2temp,"/",output,"_peaks.narrowPeak"))

    # return object with ACRs
    obj$acr <- peaks

    # return acrs
    return(obj)

}



###################################################################################################
###################################################################################################
###################################################################################################
#' buildMetaData_OG
#'
#' This function builds the meta data information for each barcode input from the Tn5 BED file.
#' Execution of this function relies on several function from Bioconductor package GenomicRanges.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges promoters
#' @importFrom IRanges subsetByOverlaps
#'
#' @param obj Object output from loadBEDandGenomeData. Required.
#' @param tss.window Size of windows flanking TSS for estimating proportion of Tn5 integration sites
#' near gene TSSs.
#' @param organelle_scaffolds Organelle scaffolds vector. If given, remove annotations from Granges object
#' which appear on organelles. Could potentially alter FRiP scores downstream as well as TSS
#' @param verbose Logical. Whether to print information.
#'
#' @rdname buildMetaData_OG
#' @export
#'
buildMetaData <- function(obj,
                          tss.window=2000,
                          organelle_scaffolds=NULL,
                          verbose=T){

    # split object
    bed <- obj$bed

    # TO DO - Write catch in case the
    # user hasn't previously run countRemoveOrganelle
    if(is.null(obj$PtMt)){
        skip.PtMt <- T
        org_val <- obj$PtMt
        gff <- obj$gff
        acr <- obj$acr
    }else{
        skip.PtMt <- F
        org_val <- obj$PtMt
        gff <- obj$gff
        acr <- obj$acr
    }

    # count number of integration sites per barcode
    if(verbose){message(" - counting Tn5 sites per barcode ...")}
    ids <- unique(as.character(bed$V4))
    counts <- table(bed$V4)

    # convert bed to Granges
    bed.gr <- GRanges(seqnames=as.character(bed$V1),
                      ranges=IRanges(start=as.numeric(bed$V2),
                                     end=as.numeric(bed$V3)),
                      strand=as.character(bed$V5),
                      names=as.character(bed$V4))

    # convert ACRs to Granges
    acr.gr <- GRanges(seqnames=as.character(acr$V1),
                      ranges=IRanges(start=as.numeric(acr$V2),
                                     end=as.numeric(acr$V3)))

    # get 1-kb up/downstream of TSS
    tss <- promoters(gff, upstream=tss.window, downstream=tss.window)
    tss <- tss[!duplicated(tss),]

    if (missing(organelle_scaffolds) | skip.PtMt) {
        tss <- promoters(gff, upstream=tss.window, downstream=tss.window)
        tss <- tss[!duplicated(tss),]
    } else {
        if(verbose){message(" - removing organelle scaffolds from annotation ...")}
        tss <- promoters(gff, upstream=tss.window, downstream=tss.window)
        tss <- tss[!duplicated(tss),]
        tss <- dropSeqlevels(tss, organelle_scaffolds, pruning.mode="tidy")
    }

    # get reads that overlap tss
    if(verbose){message(" - counting Tn5 sites at TSSs per barcode ...")}
    tss.reads <- suppressWarnings(subsetByOverlaps(bed.gr, tss, ignore.strand=T))
    tss.reads <- as.data.frame(tss.reads)
    tss.counts <- table(tss.reads$names)

    # get reads overlapping ACRs
    if(verbose){message(" - counting Tn5 sites within ACRs per barcode ...")}
    acr.reads <- suppressWarnings(subsetByOverlaps(bed.gr, acr.gr, ignore.strand=T))
    acr.reads <- as.data.frame(acr.reads)
    acr.counts <- table(acr.reads$names)

    if(verbose){message(" - finalizing meta data creation ...")}

    # merge
    counts <- counts[ids]
    tss.counts <- tss.counts[ids]
    acr.counts <- acr.counts[ids]
    if(!skip.PtMt){
        org.counts <- org_val[ids]
    }
    counts[is.na(counts)] <- 0
    tss.counts[is.na(tss.counts)] <- 0
    acr.counts[is.na(acr.counts)] <- 0
    if(!skip.PtMt){
        org.counts[is.na(org.counts)] <- 0
    }else{
        org.counts <- rep(NA, length(ids))
        names(org.counts) <- ids
    }


    df <- data.frame(cellID=ids,
                     total=as.numeric(counts),
                     tss=as.numeric(tss.counts),
                     acrs=as.numeric(acr.counts),
                     ptmt=as.numeric(org.counts),
                     row.names=ids)

    # return
    if(verbose){message("   ~ returning metadata ...")}
    obj$meta <- df
    return(obj)
}



###################################################################################################
###################################################################################################
###################################################################################################
#' findCells
#'
#' This function filters cells based on read-depth using obj mixture of user thresholds and obj spline-
#' fitting algorithm.
#'
#' @importFrom MASS kde2d
#' @importFrom viridis magma
#'
#' @param obj Object output from buildMetaData. Required.
#' @param set.tn5.cutoff Override spline fitting to set minimum tn5 count per cell.
#' Defaults to NULL.
#' @param min.cells Lower limit on the number of identified cells. Defaults to 1000.
#' @param max.cells Upper limit on the number of identified cells. Defaults to 15000.
#' @param min.tn5 Lower threshold for the minimum number of Tn5 integration sites for retaining
#' obj barcode. Defaults to 1000.
#' @param filt.org Logical. Whether or not to filter barcodes on based on proportion Tn5 sites occuring
#' in an organelle. Defaults to FALSE. Expects organelles built off meta
#' @param org.filter.thresh Remove cells with an organelle ratio (Organalle/Total_reads) greater than N
#' @param filt.tss Logical. Whether or not to filter barcodes on based on proportion Tn5 sites over-
#' lapping gene TSS. Defaults to TRUE. Filtering barcodes with this parameter is highly recommended.
#' @param tss.min.freq Float. Minimum frequency of Tn5 sites near TSSs. Defaults to 0.2.
#' @param tss.z.thresh Numeric. Z-score threshold to remove barcodes below the mean. Defeaults to 3
#' standard deviations (3).
#' @param filt.frip Logical. Whether or not to filter barcodes based on the proportion of Tn5 sites over-
#' lapping bulk ACRs. Defaults to TRUE. Using this filter is highly recommended.
#' @param frip.min.freq Float. Minimum frequency of Tn5 sites within ACRs. Defaults to 0.2.
#' @param frip.z.thresh Numeric. Z-score threshold to remove barcodes with values below X standard
#' deviations from the mean. Defaults to 3.
#' @param doplot Logical. Whether or not to plot the density scatter plots for tss/frip values.
#' @param prefix Character. Prefix output name for plots. If changed from the default (NULL), this will
#' save the images to disk as obj PDF.
#'
#' @rdname findCells_OG
#' @export
#'
findCells <- function(obj,
                      set.tn5.cutoff=NULL,
                      min.cells=1000,
                      max.cells=15000,
                      min.tn5=1000,
                      filt.org=F,
                      org.filter.thresh=0.8,
                      filt.tss=T,
                      tss.min.freq=0.2,
                      tss.z.thresh=2,
                      filt.frip=T,
                      frip.min.freq=0.2,
                      frip.z.thresh=2,
                      doplot=F,
                      prefix=NULL){
    
    # filter functions
    .filterTSS <- function(obj, min.freq=0.2, z.thresh=3, doplot=F, main=""){
        
        # get meta
        x <- subset(obj$meta.v1, obj$meta.v1$total > 0)
        
        # get tss props
        x$prop <- x$tss/x$total
        x$prop[is.na(x$prop)] <- 0
        x$zscore <- as.numeric((x$prop-mean(x$prop, na.rm=T))/sd(x$prop, na.rm=T))
        x$zscore[is.na(x$zscore)] <- 0
        x$zscore[is.infinite(x$zscore) & x$zscore < 0] <- min(x$zscore[is.finite(x$zscore)])
        x$zscore[is.infinite(x$zscore) & x$zscore > 0] <- max(x$zscore[is.finite(x$zscore)])
        x <- x[order(x$zscore, decreasing=T),]
	score <- any(x$zscore < -1*z.thresh)
	if(score){
	    prop.z.min <- max(x$prop[x$zscore < -1*z.thresh])
	    if(prop.z.min < min.freq){
		thresh <- min.freq
	    }else{
		thresh <- prop.z.min
	    }
	}else{
            thresh <- min.freq
        }
        x <- x[order(x$total, decreasing=T),]
        n.cells <- nrow(subset(x, x$tss/x$total >= thresh))
        
        # if doplot is true
        if(doplot){
            
            # get density
            den <- kde2d(log10(x$total), x$prop,
                         n=300, h=c(0.2, 0.05),
                         lims=c(range(log10(x$total)), c(0,1)))
            image(den, useRaster=T, col=c("white", rev(magma(100))),
                  xlab="Tn5 integration sites per barcode (log10)", ylab="Fraction reads TSS",
                  main=main)
            grid(lty=1, lwd=0.5, col="grey90")
            abline(h=thresh, col="red", lty=2, lwd=2)
            legend("topright", legend=paste("# cells = ", n.cells, sep=""), fill=NA, col=NA, border=NA)
            box()
        }
        
        # filter
        x$prop <- NULL
        x$zscore <- NULL
        out <- subset(x, x$tss/x$total >= thresh)
        
        # rename
        obj$meta.v2 <- out
        
        # return
        return(obj)
        
    }
    .filterFRiP <- function(obj, min.freq=0.1, z.thresh=2, doplot=F, main=""){
        
        # get meta
        x <- subset(obj$meta.v2, obj$meta.v2$total > 0)
        
        # get FRiP props
        x$prop <- x$acr/x$total
        x$prop[is.na(x$prop)] <- 0
        x$zscore <- as.numeric((x$prop-mean(x$prop, na.rm=T))/sd(x$prop, na.rm=T))
        x$zscore[is.na(x$zscore)] <- 0
        x$zscore[is.infinite(x$zscore) & x$zscore < 0] <- min(x$zscore[is.finite(x$zscore)])
        x$zscore[is.infinite(x$zscore) & x$zscore > 0] <- max(x$zscore[is.finite(x$zscore)])
        x <- x[order(x$zscore, decreasing=T),]
	score <- any(x$zscore < -1*z.thresh)
        if(score){
            prop.z.min <- max(x$prop[x$zscore < -1*z.thresh])
            if(prop.z.min < min.freq){
                thresh <- min.freq
            }else{
                thresh <- prop.z.min
            }
        }else{
            thresh <- min.freq
        }
        x <- x[order(x$total, decreasing=T),]
        n.cells <- nrow(subset(x, x$acr/x$total >= thresh))
        
        # if doplot is true
        if(doplot){
            
            # get density
            den <- kde2d(log10(x$total), x$prop,
                         n=300, h=c(0.2, 0.05),
                         lims=c(range(log10(x$total)), c(0,1)))
            image(den, useRaster=T, col=c("white", rev(magma(100))),
                  xlab="Tn5 integration sites per barcode (log10)", ylab="Fraction Tn5 insertions in ACRs",
                  main=main)
            grid(lty=1, lwd=0.5, col="grey90")
            abline(h=thresh, col="red", lty=2, lwd=2)
            legend("topright", legend=paste("# cells = ", n.cells, sep=""), fill=NA, col=NA, border=NA)
            box()
        }
        
        # filter
        x$prop <- NULL
        x$zscore <- NULL
        out <- subset(x, x$acr/x$total >= thresh)
        
        # rename
        obj$meta.v3 <- out
        
        # return
        return(obj)
    }
    .filterOrganelle <- function(obj, cell_threshold=0.8, remove_cells=FALSE, doplot=F, main=""){
        
        x <- obj$meta.v1
        x$ptmt.ratio <- x$ptmt/x$total
        
        if (remove_cells == TRUE) {
            message("... Filtering Cells based of Oragnelle Reads")
            out_meta <- subset(x, x$ptmt.ratio <= cell_threshold)
            
        } else {
            message("... Not Filtering Cells on Organelle Ratio")
            out_meta <- x
        }
        
        n.cells = nrow(out_meta)
        if(doplot){
            hist(x$ptmt.ratio, xlim=c(0,1),
                 breaks=102,
                 xlab="Organelle Ratio",
                 ylab="Counts",
                 main=main)
            abline(v=cell_threshold, col="red", lty=2, lwd=2)
            legend("topright", legend=paste("# cells = ", n.cells, sep=""), fill=NA, col=NA, border=NA)
            
        }
        
        #write back
        obj$meta.v1 <- out_meta
        
        return(obj)
    }
    
    # get meta data
    x <- obj$meta
    # order DF by total Tn5 sites
    x <- x[order(x$total, decreasing=T),]
    rank <- log10(seq(1:nrow(x)))
    depth <- log10(x$total+1)
    df <- data.frame(rank=rank, depth=depth)
    fit <- smooth.spline(rank, depth, spar=0.1)
    
    # find local minima slope
    X <- data.frame(t=seq(min(rank),max(rank),length=nrow(df)))
    Y <- predict(fit, newdata=X, deriv=1)
    xvals <- Y$x
    yvals <- Y$y
    knee <- xvals[which.min(yvals[min.cells:max.cells])]
    cells <- which.min(yvals[min.cells:max.cells])
    reads <- (10^(min(df[1:cells,]$depth))) - 1
    
    # ensure reads > min.tn5
    if(reads < min.tn5){
        reads <- min.tn5
        cells <- nrow(subset(df, df$depth>log10(reads)))
        if(cells > max.cells){
            cells <- max.cells
        }
        knee <- log10(cells)
    }
    
    # if override spline fitting
    if(!is.null(set.tn5.cutoff)){
        reads <- set.tn5.cutoff
        cells <- nrow(subset(df, df$depth>log10(reads)))
        if(cells > max.cells){
            cells <- max.cells
        }
        knee <- log10(cells)
    }
    
    # if number of cells less than threshold
    if(cells < min.cells){
        cells <- min.cells
        knee <- log10(cells)
    }
    reads <- 10^(depth[cells])
    
    
    if(filt.org == T){plot_num = 4} else {plot_num = 3}
    
    # plot
    if(doplot){
        if(!is.null(prefix)){pdf(paste0(prefix,".QC_FIGURES.pdf"), width=12, height=4)}
        layout(matrix(c(1:plot_num), nrow=1))
        plot(rank[(cells+1):length(rank)], depth[(cells+1):length(rank)],
             type="l", lwd=2, col="grey75", main=prefix,
             xlim=range(rank), ylim=range(depth),
             xlab="Barcode rank (log10)",
             ylab="Tn5 integration sites per barcode (log10)")
        lines(rank[1:cells], depth[1:cells], lwd=2, col="darkorchid4")
        grid(lty=1, lwd=0.5, col="grey90")
        abline(v=knee, col="red", lty=2)
        abline(h=log10(reads), col="red", lty=2)
        legend("bottomleft", legend=paste("# cells=",cells,", # reads=",reads,sep=""), fill=NA, col=NA, border=NA)
    }
    
    # append meta
    obj$meta.v1 <- head(x, n=cells)
    
    message("Making Dotplot")
    # filter by Organelle
    if(filt.org == T){
        obj <- .filterOrganelle(obj, remove_cells=filt.org, cell_threshold = org.filter.thresh, doplot=doplot, main="Organelle Ratio")
    }
    
    # filter by TSS
    if(filt.tss){
        obj <- .filterTSS(obj, min.freq=tss.min.freq, z.thresh=tss.z.thresh, doplot=doplot, main="% TSS")
    }
    if(filt.frip){
        if(filt.tss){
            obj <- .filterFRiP(obj, min.freq=frip.min.freq, z.thresh=frip.z.thresh, doplot=doplot, main="FRiP")
        }else{
            obj$meta.v2 <- obj$meta.v1
            obj <- .filterFRiP(obj, min.freq=frip.min.freq, z.thresh=frip.z.thresh, doplot=doplot, main="FRiP")
        }
    }
    if(!filt.tss & !filt.frip){
        obj$meta.v2 <- obj$meta.v1
        obj$meta.v3 <- obj$meta.v2
    }
    
    
    if(doplot){
        if(!is.null(prefix)){dev.off()}
    }
    # return
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' isCell
#'
#' This function compares each barcode to obj background and cell bulk reference to identify 
#' barcodes representing ambient DNA or broken nuclei. Three columns are added to the metadata:
#' background, cellbulk, is_cell, which reflect the correlation of the cell's chromatin profile with
#' the background reference, the correlation with the predicted cellbulk reference, and whether a 
#' barcode is predicted to be a cell (1 = cell, 0 = background noise).
#'
#' @importFrom qlcMatrix corSparse
#'
#' @param obj Object output from generateMatrix. Requires counts and meta data (buildMetaData) slots populated.
#' @param num.test Number of barcodes to query (ranked by total # of unique Tn5 insertions). Set to NULL to select barcodes
#' to test using obj minimum total number Tn5 insertions. Defaults to 20,000
#' @param num.tn5 Set the minimum number of Tn5 insertions to select test cells. Overridden by num.test. Defaults to NULL.
#' @param num.ref Number of cells to use as the cell bulk reference (top X cells based on # unique Tn5 insertions)
#' @param background.cutoff Maximum unique Tn5 insertions to use for selecting barcodes for the background reference set
#' @param min.pTSS Minimum per cent TSS for good cells. Defaults to 0.2.
#' @param min.FRiP Minimum fraction reads in peaks for good cells. Defaults to 0.2.
#' @param min.pTSS.z Minimum z-score from per cent TSS for good cells. Defaults to -2.
#' @param min.FRiP.z Minimum z-score from fraction reads in peaks for good cells. Defaults to -2.
#' @param verbose Default to False. Set to TRUE to progress messages. 
#' @export
#'
isCell <- function(obj, 
                   num.test=20000, 
                   num.tn5=NULL, 
                   num.ref=1000, 
                   background.cutoff=100,
                   min.pTSS=0.2, 
                   min.FRiP=0.2, 
                   min.pTSS.z= -2, 
                   min.FRiP.z= -2, 
                   verbose=F){
    
    # hidden functions
    .RowVar <- function(x) {
        spm <- t(x)
        stopifnot(methods::is(spm, "dgCMatrix"))
        ans <- sapply(base::seq.int(spm@Dim[2]), function(j) {
            if (spm@p[j + 1] == spm@p[j]) {
                return(0)
            }
            mean <- base::sum(spm@x[(spm@p[j] + 1):spm@p[j +
                                                             1]])/spm@Dim[1]
         
   sum((spm@x[(spm@p[j] + 1):spm@p[j + 1]] - mean)^2) +
                mean^2 * (spm@Dim[1] - (spm@p[j + 1] - spm@p[j]))
        })/(spm@Dim[1] - 1)
        names(ans) <- spm@Dimnames[[2]]
        ans
    }
    
    sparse_count_matrix <- obj$counts
    # make sure bins/cells are factors
    if(verbose){message(" - converting triplet format to sparseMatrix")}
    sparse_count_matrix$V1 <- factor(sparse_count_matrix$V1)
    sparse_count_matrix$V2 <- factor(sparse_count_matrix$V2)


    # convert to sparseMatrix format
    sparse_count_matrix <- Matrix::sparseMatrix(i=as.numeric(sparse_count_matrix$V1),
                              j=as.numeric(sparse_count_matrix$V2),
                              x=as.numeric(sparse_count_matrix$V3),
                             dimnames=list(levels(sparse_count_matrix$V1),levels(sparse_count_matrix$V2)))


    
    
    # select same cells 
    shared <- intersect(rownames(obj$meta), colnames(sparse_count_matrix))
    sparse_count_matrix <- sparse_count_matrix[,shared]
    obj$meta <- obj$meta[shared,]
    
    # generate stats
    obj$meta <- obj$meta[order(obj$meta$nSites, decreasing=T),]
    obj$meta$pTSS <- obj$meta$tss/obj$meta$total
    obj$meta$FRiP <- obj$meta$acrs/obj$meta$total
    obj$meta$pOrg <- obj$meta$ptmt/obj$meta$total
    
    # set initial thresholds
    if(verbose){message(" - setting filters")}
    obj$meta$tss_z <- (obj$meta$pTSS - mean(obj$meta$pTSS))/sd(obj$meta$pTSS)
    obj$meta$acr_z <- (obj$meta$FRiP - mean(obj$meta$FRiP))/sd(obj$meta$FRiP)
    obj$meta$sites_z <- (log10(obj$meta$nSites) - mean(log10(obj$meta$nSites)))/sd(log10(obj$meta$nSites))
    obj$meta$tss_z[is.na(obj$meta$tss_z)] <- -10
    obj$meta$acr_z[is.na(obj$meta$acr_z)] <- -10
    obj$meta$sites_z[is.na(obj$meta$sites_z)] <- -10 
    obj$meta$qc_check <- ifelse(obj$meta$tss_z < min.pTSS.z | obj$meta$pTSS < min.pTSS, 0, 
                              ifelse(obj$meta$acr_z < min.FRiP.z | obj$meta$FRiP < min.FRiP, 0, 1))
    
    # cells to test
    if(is.null(num.test)){
        test.set <- subset(obj$meta, obj$meta$total > num.tn5)
    }else{
        test.set <- head(obj$meta[order(obj$meta$total, decreasing=T),], n=num.test)
    }
    
    # select good cell reference
    if(verbose){message(" - parsing initial boundaries")}
    good.cells <- rownames(subset(obj$meta, obj$meta$qc_check==1))
    if(length(good.cells) > num.ref){
        good.cells <- head(good.cells, n=num.ref)
    }
    
    # select bad cell reference
    bad.cells <- rownames(subset(obj$meta, obj$meta$qc_check==0 & obj$meta$total < background.cutoff))
    
    # construct references
    gg <- sparse_count_matrix[,colnames(sparse_count_matrix) %in% good.cells]
    bb <- sparse_count_matrix[,colnames(sparse_count_matrix) %in% bad.cells]
    
    # select sites to use 
    sites <- Matrix::rowMeans(gg > 0)
    sites <- sites[order(sites, decreasing=T)]
    num.sites <- min(c(max(obj$meta$nSites), 250000))
    if(length(sites) < num.sites){
        num.sites <- length(sites)
    }
    
    # filter reference matrices
    gg <- gg[names(sites)[1:num.sites],]
    gg <- gg[,Matrix::colSums(gg) > 0]
    bb <- bb[rownames(gg),]
    bb <- bb[,Matrix::colSums(bb) > 0]
    shared <- intersect(rownames(gg), rownames(bb))
    gg <- gg[shared,]
    bb <- bb[shared,]
    gg <- gg[,Matrix::colSums(gg) > 0]
    bb <- bb[,Matrix::colSums(bb) > 0]
    
    # normalize references
    if(verbose){message(" - normalizing distributions and creating references")}
    #precent duplicate cells from being in both bad sample and test.set
    grab_cells <- unique(c(colnames(gg), colnames(bb), rownames(test.set)))
    sub.counts <- sparse_count_matrix[,grab_cells]
    sub.counts <- sub.counts[rownames(gg),]

    #Remove cells if they have no integrations around selected sites (Functionally cells with no data)
    # Addeed 9/27/2022 PM 
    cells <- Matrix::colSums(sub.counts > 0 )
    usable_cells <- cells[cells > 0]
    sub.counts <- sub.counts[,names(usable_cells)]
    
    all.res <- tfidf(list(counts=sub.counts), doL2=T)$residuals
    bb.norm <- all.res[,colnames(bb)]
    gg.norm <- all.res[,colnames(gg)]
    test.tfidf <- all.res[,names(usable_cells)] 
    
    # pick sites
    if(verbose){message(" - performing feature selection (this step is a bottle-neck and may take a while to complete)")}
    res.ave <- Matrix::rowMeans(gg.norm)
    res.res <- .RowVar(gg.norm)
    resis <- res.res/res.ave
    #resis <- loess(res.res~res.ave)$residuals
    names(resis) <- rownames(gg.norm)
    resis <- resis[order(resis, decreasing=T)]

    #If low number of sites - set to top.sites max
    if (length(resis) < 100000) {
        top.sites <- names(resis)[1:length(resis)]
    } else if (length(resis) >= 100000) {
        top.sites <- names(resis)[1:100000]
    } else {
        top.sites <- names(resis)
    }

    bb.norm <- bb.norm[top.sites,]
    gg.norm <- gg.norm[top.sites,]
    test.tfidf <- test.tfidf[top.sites,]
    
    # make bulk references, l2 norm
    bad.ref <- Matrix::rowMeans(bb.norm)
    good.ref <- Matrix::rowMeans(gg.norm)
    bad.ref <- Matrix(matrix(c(bad.ref / (sqrt(sum(bad.ref^2)))), ncol=1), sparse=T)
    good.ref <- Matrix(matrix(c(good.ref / (sqrt(sum(good.ref^2)))),ncol=1), sparse=T)
    ref <- cbind(bad.ref, good.ref)
    colnames(ref) <- c("background", "cellbulk")
    
    # prep test cells for comparisons 
    if(verbose){message(" - estimating correlations")}
    ref.cor <- corSparse(test.tfidf, ref)
    ref.cor <- as.data.frame(ref.cor)
    colnames(ref.cor) <- c("background", "cellbulk")
    rownames(ref.cor) <- colnames(test.tfidf)
    ref.cor$is_cell <- ifelse(ref.cor$cellbulk > ref.cor$background, 1, 0)
    
    # update meta data
    shared <- intersect(rownames(obj$meta), rownames(ref.cor))
    meta <- obj$meta[shared,]
    ref.cor <- ref.cor[shared,]
    updated.meta <- cbind(meta, ref.cor)
    nonrefs <- obj$meta[!rownames(obj$meta) %in% rownames(ref.cor),]
    nonrefs$background <- NA
    nonrefs$cellbulk <- NA
    nonrefs$is_cell <- NA
    updated <- rbind(updated.meta, nonrefs)
    
    # return
    obj$meta <- updated
    return(obj)
} 
###################################################################################################
###################################################################################################
###################################################################################################
#' generateMatrix
#'
#' This function generates the sparse matrix from equally sized genomic bins or ACRs.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges tileGenome
#' @importFrom IRanges subsetByOverlaps
#'
#' @param obj Object output from findCells or buildMetaData. Required.
#' @param filtered Logical. Whether or not to use the filtered set of cells. Defaults to TRUE.
#' @param windows Integer. Window size to build bins. If the 'peaks' parameter is set to TRUE,
#' this argument is over-ridden.
#' @param peaks Logical. If TRUE, use ACRs to build sparse matrix instead of genomic tiles.
#' Default is set to FALSE.
#' @param blacklist in bed format. If given removes black list regions from either peaks or 
#' generated bins. Default is set to null.
#' Default is set to FALSE.
#' @param verbose Logical. Whether or not to print progress.
#'
#' @rdname generateMatrix
#' @export
#'
generateMatrix <- function(obj,
                           filtered=T,
                           windows=1000,
                           peaks=FALSE,
                           blacklist=NULL,
                           organelle_scaffolds = NULL,
                           verbose=T){

    # convert tn5 bed to Granges
    tn5.gr <- GRanges(seqnames=as.character(obj$bed$V1),
                      ranges=IRanges(start=as.numeric(obj$bed$V2),
                                     end=as.numeric(obj$bed$V3)),
                      strand=as.character(obj$bed$V5),
                      names=as.character(obj$bed$V4))
    
    
    # Remove Organell Scaffolds if given
    if(is.null(organelle_scaffolds) == FALSE) {
        
        tn5.gr <- dropSeqlevels(tn5.gr, organelle_scaffolds)

    } else {
        tn5.gr <- tn5.gr
    }

    
    # Read in baclist if given
    if(is.null(blacklist) == FALSE) {
        blacklist_r <- read.table(as.character(blacklist))

        blacklist.gr <- GRanges(seqnames=as.character(blacklist_r$V1),
          
                ranges=IRanges(start=as.numeric(blacklist_r$V2),
                                         end=as.numeric(blacklist_r$V3)),
                          names=as.character(blacklist_r$V4))
    } else {
        blacklist.gr <- NULL
    }

    # use filtered barcodes?
    if(filtered){
        use <- obj$meta.v3
    }else{
        use <- obj$meta
    }


    # generate intervals
    if(!peaks){
        # build bins from specified tile length
        chr.seq.lengths <- as.numeric(obj$chr$V2)
        names(chr.seq.lengths) <- obj$chr$V1
        intervals <- tileGenome(chr.seq.lengths, tilewidth=windows, cut.last.tile.in.chrom=TRUE)

        #Remove if black list included
        #Remove procedure learned from: https://www.biostars.org/p/263214/
       if (is.null(blacklist.gr) == FALSE){

            intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),] 
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")


        }else{
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
        }


    }else{

        # generate intervals from ACRs
        intervals <- GRanges(seqnames=as.character(obj$acr$V1),
                             ranges=IRanges(start=as.numeric(obj$acr$V2),
                                            end=as.numeric(obj$acr$V3)))

        if (is.null(blacklist.gr) == FALSE){
            intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),] 
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
        }else{
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
        }
    }


    # get intervals overlapping Tn5 sites by barcode
    hits <- as.data.frame(findOverlaps(tn5.gr, intervals))
    df <- data.frame(regions=regions[hits$subjectHits], barcodes=as.character(obj$bed$V4)[hits$queryHits])
    df <- df[!duplicated(df),]
    df$binary <- 1
    colnames(df) <- c("V1","V2","V3")
    
    
    #9/26/2022 include for sake of calculation of isCell 
    # make sure nSites is calculated
    #Integration Sites
    a <- df
    
    #Meta data to interset
    b <- use
    a$V1 <- factor(a$V1)
    a$V2 <- factor(a$V2)

    #Generate sparse matrix
    a <- Matrix::sparseMatrix(i=as.numeric(a$V1),
                              j=as.numeric(a$V2),
                              x=as.numeric(a$V3),
                              dimnames=list(levels(a$V1),levels(a$V2)))

    # align barcodes
    both <- intersect(rownames(b), colnames(a))
    a <- a[,both]
    b <- b[both,]

    # make sure nSites is calculated
    b$nSites   <- Matrix::colSums(a)
    b$log10nSites <- log10(b$nSites)

    # return
    obj$counts <- df
    obj$meta <- b 

    return(obj) 
    
    }





###################################################################################################
###################################################################################################
###################################################################################################
#' writeReadME
#'
#' This function generates obj README of the output files from quality control
#'
#' @param fileID Character string specifying the output file name.
#'
#' @rdname writeReadME
#' @export
#'
writeReadME <- function(fileID="QC_README"){

    if(file.exists(fileID)){
        unlink(fileID)
    }
    fileConn <- file(fileID)
    writeLines(c("#########################################################################################################",
                 "# All intermediate steps are saved in the .rds file, which can be loaded in R and further analysed/partitioned:",
                 "",
                 "suppressWarnings(suppressMessages(library(GenomicRanges)))",
                 "suppressWarnings(suppressMessages(library(GenomicFeatures)))",
                 "",
                 'obj <- readRDS("youroutput.rds")',
                 "",
                 "# obj$bed contains the Tn5 insertion and barcode data",
                 "# obj$gff contains the annotation information",
                 "# obj$acr contains the narrowPeaks bulk ACRs output from MACS2 (also found in the directory called macs2_temp)",
                 "# obj$meta contains the raw barcode meta information (no barcode filtering)",
                 "# obj$meta.v1 contains the read depth filtered barcodes (step1)",
                 "# obj$meta.v2 contains the TSS proportion filtered barcodes (step2)",
                 "# obj$meta.v3 contains the FRiP filtered barcodes (step3)",
                 "",
                 '# to view these data, you can use commands on the slots, such as head(obj$bed)',
                 "",
                 "",
                 "",
                 "#########################################################################################################",
                 "# obj sparse matrix is also generated for further analysis (file_bins.sparse), compatiable with Socrates.",
                 "# The format of this matrix is as follows:",
                 "",
                 "# Column 1 = window/ACR coordinates",
                 "# Column 2 = barcode",
                 "# Column 3 = 1 # binary value indicating accessibility (0's are ommitted to save space)",
                 "",
                 "# these matrices can be loaded in R with the package 'Matrix':",
                 "",
                 'library(Matrix)',
                 'obj <- read.table("file_bins.sparse")',
                 'obj <- sparseMatrix(i=as.numeric(obj$V1), j=as.numeric(obj$V2), x=as.numeric(obj$V3), dimnames=list(levels(obj$V1),levels(obj$V2)))',
                 'head(obj[,1:5])'), fileConn)
    close(fileConn)



}


###################################################################################################
###################################################################################################
###################################################################################################
#' mergeQCdataRDS
#'
#' This function merges QC data RDS objects. Useful for joining replicates and multiple samples
#' for downstream analysis. obj simple way to aggregate multiple samples into obj single socrates
#' object is to move all the samples of interest into obj directory and load the paths with
#' list.files. IMPORTANT - only merge RDS objects where the same bins/peaks were used for each of
#' the objects. The easiest way to ensure mergeable objects is to use the same window size (instead
#' of ACRs) for the function `generateMatrix`. Use this function for rds objects saved from
#' QC analysis. The function `mergeSocratesRDS` should be specifically used for objects produced
#' using `convertSparseData`.
#'
#' # specify input files
#' > input.files <- list.files(pattern="*.rds", path="/path/to/rds/directory")
#' > input.files
#' [1] "rep1.rds" "rep2.rds"
#'
#' # name the input files
#' > names(input.files) <- c("rep1", "rep2")
#' > input.files
#'       rep1       rep2
#' "rep1.rds" "rep2.rds"
#'
#' # merge reps 1 and 2
#' > merged.obj <- mergeQCdatRDS(filenames=input.files)
#'
#' # check data structure
#' > str(merged.obj)
#'
#' @param filenames Named vector of filepaths.
#'
#' @rdname mergeQCdataRDS
#' @export
#'
mergeQCdataRDS <- function(filenames){

    # load RDS objects
    all.rds <- lapply(filenames, function(x){
        readRDS(x)
    })
    names(all.rds) <- names(filenames)

    # load filtered meta files
    meta.files <- lapply(names(filenames), function(x){
        all.rds[[x]]$meta.v3
    })
    meta.files <- as.data.frame(do.call(rbind, meta.files))

    # merge ACRs
    acrs <- lapply(names(filenames), function(x){
        all.rds[[x]]$acr
    })
    names(acrs) <- names(filenames)

    # merge BED
    beds <- lapply(names(filenames), function(x){
        all.rds[[x]]$bed
    })
    names(beds) <- names(filenames)

    # get gff
    gff <- all.rds[[1]]$gff

    # get chr
    chr <- all.rds[[1]]$chr

    # merge counts data
    cnts <- lapply(names(filenames), function(x){
        ct <- all.rds[[x]]$counts
        colnames(ct) <- c("V1", "V2", "V3")
        ct$V1 <- factor(ct$V1)
        ct$V2 <- factor(ct$V2)
        ct$V3 <- as.numeric(ct$V3)
        ct
    })
    cnts <- as.data.frame(do.call(rbind, cnts))
    cnts$V1 <- factor(cnts$V1)
    cnts$V2 <- factor(cnts$V2)
    cnts$V3 <- as.numeric(cnts$V3)
    cnts <- sparseMatrix(i=as.numeric(cnts$V1),
                         j=as.numeric(cnts$V2),
                         x=as.numeric(cnts$V3),
                         dimnames=list(levels(cnts$V1), levels(cnts$V2)))

    # shared ids only
    shared.ids <- intersect(rownames(meta.files),
                            colnames(cnts))
    cnts <- cnts[,shared.ids]
    cnts <- cnts[Matrix::rowSums(cnts)>0,]
    meta.files <- meta.files[shared.ids,]

    # build new object
    m.obj <- list(meta=meta.files,
                  counts=cnts,
                  bed=beds,
                  gff=gff,
                  chr=chr,
                  samples=names(filenames))

    # return
    return(m.obj)

}


###################################################################################################
###################################################################################################
###################################################################################################
#' mergeSocratesRDS
#'
#' This function merges Socrates data RDS objects. Useful for joining replicates and multiple samples
#' for downstream analysis. obj simple way to aggregate multiple samples into obj single socrates
#' object is to move all the samples of interest into obj directory and load the paths with
#' list.files. IMPORTANT - only merge RDS objects where the same bins/peaks were used for generating
#' the sparse matrices. The easiest way to ensure mergeable objects is to use the same window size
#' (instead of ACRs) for the function `generateMatrix` from each sample.
#'
#' # specify input files
#' > input.files <- list.files(pattern="*.rds", path="/path/to/rds/directory")
#' > input.files
#' [1] "rep1.rds" "rep2.rds"
#'
#' # name the input files
#' > names(input.files) <- c("rep1", "rep2")
#' > input.files
#'       rep1       rep2
#' "rep1.rds" "rep2.rds"
#'
#' # merge reps 1 and 2
#' > merged.obj <- mergeSocratesRDS(filenames=input.files)
#'
#' # check data structure
#' > str(merged.obj)
#'
#' # one can also load using obj list of socrates objects.
#' > all.soc.obj <- list(rep1=soc.obj1, rep22=soc.obj2)
#' > merged.obj <- mergeSocratesRDS(obj.list=all.soc.obj)
#'
#' @param filenames Named vector of filepaths.
#' @param obj.list Named list of Socrates objects output from `convertSparseData`.
#'
#' @rdname mergeSocratesRDS
#' @export
#'
mergeSocratesRDS <- function(filenames=NULL, obj.list=NULL){

    # checks
    if(is.null(filenames)){
        if(is.null(obj.list)){
            stop("arguments 'filenames' or 'obj.list' must not be NULL ...")
        }
    }

    # load file RDS
    if(!is.null(filenames)){

        # load RDS objects
        all.rds <- lapply(filenames, function(x){
            readRDS(x)
        })
        names(all.rds) <- names(filenames)
    }else{
        all.rds <- obj.list
        filenames <- names(all.rds)
        names(filenames) <- names(all.rds)
        rm(obj.list)
    }

    # load filtered meta files
    row.ids <- c()
    its <- 0
    meta.files <- lapply(filenames, function(x){
        its <<- its + 1
        dat <- all.rds[[x]]$meta
        dat$sampleID <- names(filenames)[its]
        row.ids <<- c(row.ids, rownames(dat))
        dat
    })
    meta.files <- as.data.frame(do.call(rbind, meta.files))
    rownames(meta.files) <- row.ids

    # merge counts data
    cnts <- lapply(filenames, function(x){
        ct <- as.data.frame(summary(all.rds[[x]]$counts))
        colnames(ct) <- c("V1", "V2", "V3")
        ct$V1 <- rownames(all.rds[[x]]$counts)[ct$V1]
        ct$V2 <- colnames(all.rds[[x]]$counts)[ct$V2]
        ct$V1 <- factor(ct$V1)
        ct$V2 <- factor(ct$V2)
        ct
    })
    cnts <- as.data.frame(do.call(rbind, cnts))
    cnts$V1 <- as.factor(cnts$V1)
    cnts$V2 <- as.factor(cnts$V2)
    cnts$V3 <- as.numeric(cnts$V3)
    cnts <- sparseMatrix(i=as.numeric(cnts$V1),
                         j=as.numeric(cnts$V2),
                         x=as.numeric(cnts$V3),
                         dimnames=list(levels(cnts$V1), levels(cnts$V2)))

    # shared ids only
    shared.ids <- intersect(rownames(meta.files),
                            colnames(cnts))
    cnts <- cnts[,shared.ids]
    cnts <- cnts[Matrix::rowSums(cnts)>0,]
    meta.files <- meta.files[shared.ids,]

    # build new object
    m.obj <- list(meta=meta.files,
                  counts=cnts)

    # return
    return(m.obj)

}



