###################################################################################################
###################################################################################################
###################################################################################################
#' coAccess
#'
#' This function loads the sparse matrix and meta data and returns a list - counts & meta
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#' @importFrom cicero
#'
#' @param input Path to sparse cell x peak text file in triplet format. Input data can be gzipped.
#' Required.
#' @param meta Path to meta data for cells in sparse cell x peak data set in tsv format. Headers
#' and barcode/cellID rownames are required. Meta data can be gzipped. Providing meta data is
#' highly recommended, but not required.
#' @param verbose logical. Defaults to FALSE.
#'
#' @rdname coAccess
#' @export
#'
