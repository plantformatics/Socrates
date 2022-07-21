###################################################################################################
###################################################################################################
###################################################################################################
#'  Visualized imputed markers 
#'
#' This function calculates the accessibilities near annotated genes per cells.  
#' Can be run after loadBEDandGenomeData function.
#'
#' @param obj Socrates object
#' @param meta_slot slot holding meta data to use for plotting
#' @param imputed_slot slot holding hte imputed data outbput by function acc.imputation
#' @param marker_slot marker info in table
#' @param marker_slot markers either the string "all" or vector of marker genes
#' by name which to plot. Note that for ALL a pdf file output is reccomended as
#' the output is LARTGE

viz_markers <- function(soc.obj,
                        meta_slot = "Clusters",
                        imputed_slot= "impute.acc",
                       marker_slot = "markers", 
                       markers="all",
                       output_file = NULL,
                       output_dir = "."){
    
    #Check if the imputation slot exits
    if ( rlang::has_name(soc.obj, imputed_slot) == TRUE) {
        message("Imputed data present.. Working....") 
        imputed_sparse <- (soc.obj[[imputed_slot]])
    } else if ( rlang::has_name(soc.obj, imputed_slot) == FALSE) {
        message(paste0("Slot name: ", imputed_slot, " Does not exits in the passed object... Cancelling..."))
        stop()
    } else {
        message("Error...")
    }
    
    #Check if marker slot exits
    if ( rlang::has_name(soc.obj, marker_slot) == TRUE) {
        message("Markers present... Working....")
        gene_markers <- soc.obj[[marker_slot]]
    } else if ( rlang::has_name(soc.obj, marker_slot) == FALSE) {
        message(paste0("Slot name: ", marker_slot, " Does not exits in the passed object... Need markers to visualize"))
        stop()
    } else {
        message("Error...")
    }
                
    #Check meta_slot
    if ( rlang::has_name(soc.obj, meta_slot) == TRUE) {
        message("Meta data present... Working....")
        meta_data <- soc.obj[[meta_slot]]
    } else if ( rlang::has_name(soc.obj, meta_slot) == FALSE) {
        message(paste0("Slot name: ", meta_slot, " Does not exits in the passed object... Is this the correct meta data slot?"))
        stop()
    } else {
        message("Error...")
    }
    
        
    gene_markers <- gene_markers  %>%
        dplyr::arrange(type)

    # Select markers from the imputed data. If all - use all makers stored in the slot, if not, isolate down.
    if (markers == "all") {
        selected_markers <- gene_markers$geneID
        selected_names <- gene_markers$name
        message("Plotting all markers found in marker slot")
        selected_markers <- selected_markers[selected_markers %in% rownames(imputed_sparse)]
        selected_markers <- unique(selected_markers[selected_markers %in% rownames(imputed_sparse)])
     } else {
        
        filtered_markers <- gene_markers  %>% 
            dplyr::filter(name %in% markers)
        
        selected_markers <- filtered_markers$geneID
        selected_names <- filtered_markers$name
        marker_string <- paste(selected_names)
        message(paste0("Plotting a subset of markers given: ", marker_string))
        
        selected_markers <- selected_markers[selected_markers %in% rownames(imputed_sparse)]
        #selected_markers <- unique(selected_markers[selected_markers %in% rownames(imputed_sparse)])
        print(selected_markers)
    }   
    
    
    gathered_list <- list()
    take_length_markers <- seq(1,length(selected_markers))
    
    
    # If running all markers - adding to PDF image and save output using prefix                         
    print(take_length_markers)
    for ( z in take_length_markers ) {
        i <- selected_markers[[z]]
        message(paste0("Plotting Marker", i))

        storage_character <- as.character(z)
        grab_gene_info <- gene_markers[gene_markers$geneID == i, ]
        pull_sparse_matrix_gene <- imputed_sparse[i,]
        imputer_sparse_2_filtered <- as_tibble(pull_sparse_matrix_gene, rownames = "cellID")


        imputated_sparse_2_join <- left_join(meta_data, imputer_sparse_2_filtered, by = c("cellID"), copy = TRUE) %>%
            replace(is.na(.), 0)

        #cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
        cols <- colorRampPalette(c("grey80","#DD3497", "#49006A"), interpolate="linear", bias = .5)(100)
        upper.lim <- quantile(imputated_sparse_2_join$value, .90, na.rm = TRUE) + 1e-6

        graphable_sparse_2 <- imputated_sparse_2_join  %>%
        dplyr::mutate(final_ac = case_when(value >= upper.lim ~ upper.lim,
                                           TRUE ~ value))


        grad <- scale_colour_gradientn(colours = cols,
            limits = c(0, max(graphable_sparse_2$final_ac)))

        generate_title <- str_c(grab_gene_info$name, grab_gene_info$geneID, sep ="\n")
        arranged_imputation_plot <- graphable_sparse_2  %>%
        arrange(final_ac)  %>%
        ggplot(., aes(umap1, umap2, colour = final_ac)) +
            geom_jitter(position = "jitter", size = .5, stroke = 0, shape = 16) +
            scale_fill_continuous(type="viridis") +
            grad +
            theme_classic() +
            #theme(legend.position="none") +
            ggtitle(generate_title)

        gathered_list[[storage_character]] <- arranged_imputation_plot

    }


    print("Plotting All Markers in Grid...")
    captured_final_plot <- plot_grid(plotlist = gathered_list, ncol = 6)
    width_cal <- 6 * 5
    length_cal <- (length(selected_markers)/6 * 5)

    if ( is.null(output_file) == TRUE ) {
        captured_final_plot
        return(captured_final_plot)
    } else {
        
        output_file_name = paste0(output_dir, "/", output_file)
        #output_name <- str_c(prefix, "pdf", sep = ".")

        ggsave(output_file_name, plot = captured_final_plot,
            width = width_cal, height = length_cal,
            units = c('in'), limitsize = FALSE, dpi = 300)
        
        
                
    }
}
