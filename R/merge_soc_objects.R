###################################################################################################
###################################################################################################
###################################################################################################
#' mergeSocObjects
#'
#' This function will merge newly formed Socrates objects across multiple samples. 
#' The function assumes that the meta data for each soc object has a column with a unique identifier
#' specifying the sample ID, such as `library`. Will remove all slots except `meta` and `counts`. 
#' This should be run prior to TFIDF or regModel. 
#'
#' @import Matrix
#' @param obj.list a list of Socrates objects. 
#'
#' @rdname mergeSocObjects
#' @export
#'
mergeSocObjects <- function(obj.list){
    
    # functions
    .merge.sparse <- function(cnt.list) {
        
        cnnew <- character()
        rnnew <- character()
        x <- vector()
        i <- numeric()
        j <- numeric()
        
        for (M in cnt.list) {
            
            cnold <- colnames(M)
            rnold <- rownames(M)
            
            cnnew <- union(cnnew,cnold)
            rnnew <- union(rnnew,rnold)
            
            cindnew <- match(cnold,cnnew)
            rindnew <- match(rnold,rnnew)
            ind <- summary(M)
            i <- c(i,rindnew[ind[,1]])
            j <- c(j,cindnew[ind[,2]])
            x <- c(x,ind[,3])
        }
        
        sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
    }
    
    # separate counts and meta
    counts <- lapply(obj.list, function(x){
        x$counts
    })
    counts <- .merge.sparse(counts)
    
    # meta
    metas <- lapply(obj.list, function(x){
        x$meta
    })
    metas <- do.call(rbind, metas)
    rownames(metas) <- metas$cellID
    #metas$library <- data.frame(do.call(rbind, strsplit(rownames(metas), "-")))[,2]
    metas <- metas[colnames(counts),]
    
    # new object
    new.obj <- list(counts=counts, meta=metas)
    return(new.obj)
    
}

