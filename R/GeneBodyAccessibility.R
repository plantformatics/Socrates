## The obj should have "bed" ,"ann" (gtf) information from loadBEDandGenomeData function.

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

