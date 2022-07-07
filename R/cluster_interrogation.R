###################################################################################################
###################################################################################################
###################################################################################################
#' Run a modified wilcoxon-rank-sum test on the chosen clusters. Reports back the most descriptive genes
#' associated with each cluster. Find detials here: https://github.com/immunogenomics/presto
#'
#'
#' @param meta_obj Socrates object. Should have meta data and cluster labels, as well as accessability as gene accessability 
#' @param meta_slot which slot of meta data associated with the Socrates object you would like to run Socrates on.
#' @param cluster_name The cluster slot name which you wish to run presto on. Needs to be in the meta data to be used).
#' @param counts_raw Accessability counts slow as calculated by GeneBodyAccessibility function.
run_presto <- function(meta_obj,
                       meta_slot = "Clusters", 
                       cluster_name = "LouvainClusters", 
                       counts_raw="acr_counts_raw"){
  ### Pull the correct meta datafrom the Socrates Object
  sparse_matrix <- meta_obj[[counts_raw]]
  meta_data <- meta_obj[[meta_slot]]
  
  
  ### Read the Giant Sparse Matrix
  ### Old from when calculations were done with python script
  #loaded_sparse_matric <- read_delim(sparse_matrix, delim='\t', 
  #                                   col_names = c("gene_name", "barcode", "accessability"), 
  #                                   col_types = "ccn")
  
    
  message("Loading Sparse Matrix for marker testing")
  loaded_sparse_matric <- sparse_matrix
  
  combined_large_w_sparse <- loaded_sparse_matric  %>% 
    filter(gene_name != "Annotation")  %>% 
    dplyr::select(gene_name, barcode, accessability)  %>% 
    left_join(., meta_data, by = c("barcode" = "cellID"))  %>% 
    filter(is.na(!!sym(cluster_name)) != TRUE)  %>% 
    arrange(!!sym(cluster_name))
  
  #order <- combined_large_w_sparse[!!sym(cluster_name)]
  
  combined_large_w_sparse <- combined_large_w_sparse  %>% 
    dplyr::select(gene_name, barcode, accessability)  %>% 
    mutate(across(accessability, as.numeric))
  
  
  gene_names <- unique(combined_large_w_sparse$gene_name)
  barcodes <- unique(combined_large_w_sparse$barcode)
  
  combined_large_w_sparse$row <- match(combined_large_w_sparse$gene_name, gene_names)
  combined_large_w_sparse$col <- match(combined_large_w_sparse$barcode, barcodes)
  
  
  UIMatrix <- sparseMatrix(i = combined_large_w_sparse$row,
                           j = combined_large_w_sparse$col,
                           x = combined_large_w_sparse$accessability,
                           dimnames=list(gene_names, barcodes))
  
  barcode_tibble <- tibble(barcode = unlist(UIMatrix@Dimnames[2])) %>% 
    left_join(., as_tibble(meta_data), by = c("barcode" = "cellID"))  %>% 
    filter(is.na(total) != TRUE)
  
  clust_name <- c(cluster_name)
  garbbed_louv_order <- pull(barcode_tibble,!!sym(clust_name))
  
  sparse_matrix_test <- wilcoxauc(UIMatrix, garbbed_louv_order)
  
  meta_obj[["presto_marker"]] <- sparse_matrix_test
  
  return(meta_obj)
}

