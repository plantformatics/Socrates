###################################################################################################
###################################################################################################
###################################################################################################
#' loadBEDandGenomeData
#'
#' This function loads a 5 column BED file specifying the single-bp Tn5 integration sites.
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
#'
#' @rdname loadBEDandGenomeData
#' @export
#'
loadBEDandGenomeData <- function(bed, ann, sizes, verbose=T){

    # load hidden pre-check function
    .preRunChecks <- function(bed, ann, chr, verbose=T){

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
        if(!file.exists(chr)){
            message("Chromosome sizes file, ", chr," does not exist... exiting")
            quit(save="no")
        }else{
            if(verbose){message("Chromosome sizes file path = ", chr, " ... ok")}
        }

        # check for macs2
        prg <- "macs2"
        if(Sys.which(prg) == ""){
            stop(message("Please install ", prg, ". ACR calling will not work without MACS2 in your path."))
        }else{
            if(verbose){message("Macs2 is installed .... ok")}
        }


    }

    # run pre-checks
    .preRunChecks(bed, ann, chr, verbose=verbose)

    # save paths
    bedpath <- bed
    annpath <- ann
    chrpath <- chr

    # load args
    if(verbose){message(" - loading data (this may take a while for big BED files) ...")}
    if(grepl(".gz$", bed)){
        a <- read.table(gzfile(as.character(bed)))
    }else{
        a <- read.table(as.character(bed))
    }

    # load Gff
    gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format="gff3", dbxrefTag="Parent")))
    chrom <- read.table(as.character(sizes))

    if(verbose){message(" - finished loading data")}
    return(list(bed=a, gff=gff, chr=chrom, bedpath=bedpath, annpath=annpath, chrpath=chrpath))

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
#' buildMetaData
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
#' @param verbose Logical. Whether to print information.
#'
#' @rdname buildMetaData
#' @export
#'
buildMetaData <- function(obj, tss.window=2000, verbose=T){

    # split object
    bed <- obj$bed
    gff <- obj$gff
    acr <- obj$acr

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

    # merge
    counts <- counts[ids]
    tss.counts <- tss.counts[ids]
    acr.counts <- acr.counts[ids]
    counts[is.na(counts)] <- 0
    tss.counts[is.na(tss.counts)] <- 0
    acr.counts[is.na(acr.counts)] <- 0
    df <- data.frame(cellID=ids,
                     total=as.numeric(counts),
                     tss=as.numeric(tss.counts),
                     acrs=as.numeric(acr.counts),
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
#' This function filters cells based on read-depth using a mixture of user thresholds and a spline-
#' fitting algorithm.
#'
#' @importFrom MASS kde2d
#' @importFrom viridis magma
#'
#' @param obj Object output from buildMetaData. Required.
#' @param set.tn5.cutoff Override spline fitting to set minimum tn5 count per cell.
#' Defaults to NULL.
#' @param min.cells Lower limit on the number of identified cells. Defaults to 100.
#' @param max.cells Upper limit on the number of identified cells. Defaults to 15000.
#' @param min.tn5 Lower threshold for the minimum number of Tn5 integration sites for retaining
#' a barcode. Defaults to 1000.
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
#' save the images to disk as a PDF.
#'
#' @rdname findCells
#' @export
#'
findCells <- function(obj,
                      set.tn5.cutoff=NULL,
                      min.cells=100,
                      max.cells=15000,
                      min.tn5=1000,
                      filt.tss=T,
                      tss.min.freq=0.2,
                      tss.z.thresh=3,
                      filt.frip=T,
                      frip.min.freq=0.1,
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
        if(min(x$zscore) < (-1*z.thresh) & min(x$prop) > min.freq){
            below <- subset(x, x$zscore < (-1*z.thresh) & x$prop > min.freq)
            thresh <- max(below$prop, na.rm=T)
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
        if(min(x$zscore) < (-1*z.thresh) & min(x$prop) > min.freq){
            below <- subset(x, x$zscore < (-1*z.thresh) & x$prop > min.freq)
            thresh <- max(below$prop, na.rm=T)
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
    knee <- xvals[which.min(yvals[1:max.cells])]
    cells <- which.min(yvals[1:max.cells])
    reads <- min(df[1:cells,]$total)

    # ensure reads > min.tn5
    if(reads < min.tn5){
        reads <- min.tn5
        cells <- nrow(subset(df, df$depth>log10(reads)))
        knee <- log10(cells)
    }

    # if override spline fitting
    if(!is.null(set.tn5.cutoff)){
        reads <- set.tn5.cutoff
        cells <- nrow(subset(df, df$depth>log10(reads)))
        knee <- log10(cells)
    }

    # if number of cells less than threshold
    if(cells < min.cells){
        cells <- min.cells
        knee <- log10(cells)
    }
    reads <- 10^(depth[cells])

    # plot
    if(doplot){
        if(!is.null(prefix)){pdf(paste0(prefix,".QC_FIGURES.pdf"), width=12, height=4)}
        layout(matrix(c(1:3), nrow=1))
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
#' @param verbose Logical. Whether or not to print progress.
#'
#' @rdname generateMatrix
#' @export
#'
generateMatrix <- function(obj,
                           filtered=T,
                           windows=1000,
                           peaks=FALSE,
                           verbose=T){

    # convert tn5 bed to Granges
    tn5.gr <- GRanges(seqnames=as.character(obj$bed$V1),
                      ranges=IRanges(start=as.numeric(obj$bed$V2),
                                     end=as.numeric(obj$bed$V3)),
                      strand=as.character(obj$bed$V5),
                      names=as.character(obj$bed$V4))

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
        regions <- as.data.frame(intervals)
        regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")

    }else{

        # generate intervals from ACRs
        intervals <- GRanges(seqnames=as.character(obj$acr$V1),
                             ranges=IRanges(start=as.numeric(obj$acr$V2),
                                            end=as.numeric(obj$acr$V3)))
        regions <- as.data.frame(intervals)
        regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
    }

    # get intervals overlapping Tn5 sites by barcode
    hits <- as.data.frame(findOverlaps(tn5.gr, intervals))
    df <- data.frame(regions=regions[hits$subjectHits], barcodes=as.character(obj$bed$V4)[hits$queryHits])
    df <- df[!duplicated(df),]
    df$binary <- 1
    colnames(df) <- c("V1","V2","V3")

    # return
    obj$counts <- df
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' writeReadME
#'
#' This function generates a README of the output files from quality control
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
                 'a <- readRDS("youroutput.rds")',
                 "",
                 "# a$bed contains the Tn5 insertion and barcode data",
                 "# a$gff contains the annotation information",
                 "# a$acr contains the narrowPeaks bulk ACRs output from MACS2 (also found in the directory called macs2_temp)",
                 "# a$meta contains the raw barcode meta information (no barcode filtering)",
                 "# a$meta.v1 contains the read depth filtered barcodes (step1)",
                 "# a$meta.v2 contains the TSS proportion filtered barcodes (step2)",
                 "# a$meta.v3 contains the FRiP filtered barcodes (step3)",
                 "",
                 '# to view these data, you can use commands on the slots, such as head(a$bed)',
                 "",
                 "",
                 "",
                 "#########################################################################################################",
                 "# A sparse matrix is also generated for further analysis (file_bins.sparse), compatiable with Socrates.",
                 "# The format of this matrix is as follows:",
                 "",
                 "# Column 1 = window/ACR coordinates",
                 "# Column 2 = barcode",
                 "# Column 3 = 1 # binary value indicating accessibility (0's are ommitted to save space)",
                 "",
                 "# these matrices can be loaded in R with the package 'Matrix':",
                 "",
                 'library(Matrix)',
                 'a <- read.table("file_bins.sparse")',
                 'a <- sparseMatrix(i=as.numeric(a$V1), j=as.numeric(a$V2), x=as.numeric(a$V3), dimnames=list(levels(a$V1),levels(a$V2)))',
                 'head(a[,1:5])'), fileConn)
    close(fileConn)



}


###################################################################################################
###################################################################################################
###################################################################################################
#' mergeQCdataRDS
#'
#' This function merges QC data RDS objects. Useful for joining replicates and multiple samples
#' for downstream analysis. A simple way to aggregate multiple samples into a single socrates
#' object is to move all the samples of interest into a directory and load the paths with
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
#' for downstream analysis. A simple way to aggregate multiple samples into a single socrates
#' object is to move all the samples of interest into a directory and load the paths with
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
#' # one can also load using a list of socrates objects.
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



