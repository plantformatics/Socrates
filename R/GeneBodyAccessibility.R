###################################################################################################
###################################################################################################
###################################################################################################
#' Reload a bed file of Tn5 integration events into a socrates obj
#'
#' This function re-inserts the bed file back into the socrates object if it
#' has been missing in prior analysis
#'
#' @param obj Socrates object 
#' @param bed file with tn5 integrations mapped to 1bp
reload_bed_file <- function(obj, bed) {
    
    if(grepl(".gz$", bed)){
        obj$bed <- read.table(gzfile(as.character(bed)))
    }else{
        obj$bed <- read.table(as.character(bed))}
    
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Reload a gtf or gff3 file into a socrates object
#'
#' This function re-inserts the gff3 file back into the socrates object if it
#' has been missing in prior analysis
#'
#' @param obj Socrates object 
#' @param ann Annotation file in gtf or gff3 format.
reload_gff_file <- function(obj, ann) {
    
    if(grepl(".gtf", ann)){
        anntype <- "gtf"
    }else{
        anntype <- "gff3"
    }
        
    gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format=anntype, dbxrefTag="Parent")))
    
    obj$gff <- gff
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Find the genes near Tn5 insertion
#'
#' This function calculates the accessibilities near annotated genes per cells.  
#' Can be run after loadBEDandGenomeData function.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @param obj Object output from loadBEDandGenomeData. Required.
#' @param FeatureName character, Must be one of "gene" or "transcript". Default is "gene".
#' @param nRange numeric, Extension length for annotation. It will increase the finding regions backward and forward 'nRange'bp of genes.


GeneBodyAccessibility <- function(obj, FeatureName="gene", nRange=500){

    ## Load Tn5 file with GRanges
    Tn5_Grange <-  GRanges(seqnames = obj$bed$V1,
                ranges = IRanges(start = obj$bed$V2,
                                 end = obj$bed$V3,
                                 names = obj$bed$V4))
    
    ## Select Features
    if (FeatureName == "gene"){
      Ann <- genes(obj$gff)
    } else if(FeatureName == "transcript"){
      Ann <- transcripts(obj$gff)
    } else{message("ERROR: Feature name should be 'gene' or 'transcript")}

    
    ## Broaden annotation
    Start_new <- c((ranges(Ann)+nRange)@start)-1 ## To convert to the bed format coordinate
    Width_new <- c((ranges(Ann)+nRange)@width)+1 ## To convert to the bed format coordinate
    BroadRange_Ann <- GRanges(seqnames=Ann@seqnames,
                              ranges= 
                                IRanges(start=Start_new,
                                        width=Width_new,
                                        names=names(Ann)))
    #obj$ann_broad <- BroadRange_Ann

    ## Find overlap
    hits_Within <- findOverlaps(Tn5_Grange,  BroadRange_Ann,minoverlap=1,
                                type=c("within"),select="all",ignore.strand = TRUE)
    
    Intersect <- paste(names(BroadRange_Ann)[hits_Within@to],
                       names(Tn5_Grange)[hits_Within@from],sep="/_Com_/")
    
    Intersect <- table(Intersect)
    Intersect<- data.frame(gene_name = as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                 split="/_Com_/"), "[", 1)),
                        barcode= as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                              split="/_Com_/"), "[", 2)),
                        accessability = as.character(Intersect))
    
    obj$sc_gene_ac <- Intersect

    return(obj)}

