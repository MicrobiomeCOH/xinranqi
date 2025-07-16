
run_differential_abundance <- function(dt_input_feature,
                                       dt_input_metadata,
                                       method = method, #c("deseq2", "ancombc", "lefse")
                                       group_field = group_field,
                                       save_path = save_path, 
                                       is.relative.abundance = FALSE) {
  # Load required packages
  library(data.table)
  library(tidyr)
  library(ANCOMBC)
  library(DESeq2)
  library(lefser)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  
  # Step 1: Reshape dt_input_feature to wide format
  abundance_wide <- as.data.frame(dcast(dt_input_feature, feature_type ~ sampleid, value.var = "feature_value", fill = 0))
  abundance_wide$feature_type[is.na(abundance_wide$feature_type)] <- "NA"
  rownames(abundance_wide) <- abundance_wide$feature_type
  #abundance_wide <- abundance_wide[, -1]
  
  # Step 2: Format metadata for ANCOM-BC
  metadata <- dt_input_metadata[, .(sampleid, group = get(group_field))]
  metadata <- unique(metadata)
  metadata <- metadata[sampleid %in% colnames(abundance_wide)]
  
  # Ensure ordering
  abundance_wide <- abundance_wide[, metadata$sampleid, drop = FALSE]
  #abundance_wide$Feature <- rownames(abundance_wide)
  
  # Step 3: Convert to phyloseq object
  library(phyloseq)
  
  OTU <- if (is.relative.abundance) {otu_table(as.matrix(abundance_wide/100), taxa_are_rows = TRUE)} else {otu_table(as.matrix(abundance_wide), taxa_are_rows = TRUE)}
  SAM <- sample_data(data.frame(row.names = metadata$sampleid, group = metadata$group))
  physeq <- phyloseq(OTU, SAM)
  
  # Extract genus names from OTU table row names
  taxa_names <- taxa_names(physeq)
  tax_table_mat <- matrix(taxa_names, ncol = 1)
  rownames(tax_table_mat) <- taxa_names
  colnames(tax_table_mat) <- "Genus"  # this can be "Feature" or "Taxon" too
  
  # Add dummy taxonomy table to physeq
  tax_table(physeq) <- tax_table(tax_table_mat)
  
  if (method == "ancombc") {
    # Step 4: Run ANCOM-BC
    if (is.relative.abundance) {
    out <- ancombc2(
      data = physeq,
      assay_name = "relab",         # because you're using relative abundance
      tax_level = "Genus",        # or "Genus" if your taxa are annotated
      fix_formula = "group",        # matches your sample metadata column
      rand_formula = NULL,
      group = "group",   
      p_adj_method = "fdr",
      prv_cut = 0.10,               # you can adjust this
      lib_cut = 0,
      struc_zero = TRUE,
      neg_lb = TRUE,
      pseudo = 1e-06,
      pseudo_sens = TRUE
    )} else {
      out <- ancombc2(
        data = physeq,
        tax_level = "Genus",        # or "Genus" if your taxa are annotated
        fix_formula = "group",        # matches your sample metadata column
        rand_formula = NULL,
        group = "group",   
        p_adj_method = "fdr",
        prv_cut = 0.10,               # you can adjust this
        lib_cut = 0,
        struc_zero = TRUE,
        neg_lb = TRUE,
        pseudo = 1e-06,
        pseudo_sens = TRUE
      )
    }

    # Step 5: Extract and save results
    res_df <- out$res
    
    # Extract relevant columns
    pq_df <- res_df[, sapply(colnames(res_df), function(x){x=="taxon" | {any(unlist(strsplit(x, "_")) %in% c("p","q")) & all(unlist(strsplit(x, "_")) != "(Intercept)")}})]
    
    # Optionally rename the columns for clarity
    colnames(pq_df) <- c("taxon", "p_value", "q_value")
    
    # Save result table
    fwrite(pq_df, file.path(save_path, "ancombc_differential_abundance_results.tsv"), sep = "\t")
    
    # Step 6: Volcano plot
    plot_volcano <- function(res_df, p_threshold = 0.05) {
      # Build the column names dynamically
      lfc_col <- colnames(res_df)[sapply(colnames(res_df), function(x){any(unlist(strsplit(x, "_"))=="lfc") & all(unlist(strsplit(x, "_")) != "(Intercept)") })]
      p_col   <- colnames(res_df)[sapply(colnames(res_df), function(x){any(unlist(strsplit(x, "_"))=="p") & all(unlist(strsplit(x, "_")) != "(Intercept)") })]
      ss_col  <- colnames(res_df)[sapply(colnames(res_df), function(x){any(unlist(strsplit(x, "_"))=="passed") & all(unlist(strsplit(x, "_")) != "(Intercept)") })]
      
      # Create significance flag
      res_df$significant <- with(res_df, res_df[[p_col]] < p_threshold & res_df[[ss_col]])
      
      # Generate plot
      p <- ggplot(res_df, aes(x = .data[[lfc_col]], y = -log10(.data[[p_col]]), color = significant)) +
        geom_point(alpha = 0.8, size = 2) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dotted", color = "red") +
        scale_color_manual(values = c("gray", "red")) +
        labs(
          title = paste("Volcano Plot:", unlist(strsplit(p_col, "group"))[2]),
          x = paste0("Log Fold Change (", unlist(strsplit(p_col, "group"))[2], " vs. reference)"),
          y = paste0("-log10(p-value)"),
          color = "Significant"
        ) +
        theme_minimal()
      
      return(p)
    }
    
    volcano <- plot_volcano(res_df)
    
    ggsave(file.path(save_path, "ancombc_differential_abundance_volcano_plot.png"), plot = volcano, width = 8, height = 6, dpi = 300)
  }
  
  if (method == "deseq2" & is.relative.abundance==FALSE) {
    # Convert to DESeq2 object (replace "Group" with your variable)
    dds <- phyloseq_to_deseq2(physeq, ~ group)
    
    # Estimate size factors and run DESeq
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds, fitType = "parametric")
    
    # Get differential abundance results (replace levels as needed)
    res <- results(dds)
    
    # Convert to data frame and add taxonomic info
    res_df <- as.data.frame(res)
    res_df$taxa <- rownames(res_df)

    # Add significance column
    res_df <- res_df %>%
      mutate(Significance = case_when(
        pvalue < 0.05 & log2FoldChange > 1 ~ "Up",
        pvalue < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "Not Significant"
      ))
    
    volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
      geom_point(alpha = 0.8, size = 2) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "gray")) +
      theme_minimal() +
      labs(title = "Volcano Plot of Differential Abundance",
           x = "Log2 Fold Change",
           y = "-log10(p-value)") +
      theme(legend.position = "right")
    
    ggsave(file.path(save_path, "deseq2_differential_abundance_volcano_plot.png"), plot = volcano, width = 8, height = 6, dpi = 300)
    
    fwrite(res_df, file.path(save_path, "deseq2_differential_abundance_results.tsv"), sep = "\t")
  }
  
  if (method == "lefse" & is.relative.abundance==FALSE) {
    # Step 1: Extract the OTU matrix (features = rows, samples = columns)
    otu <- as(otu_table(physeq), "matrix")
    if (!taxa_are_rows(physeq)) {
      otu <- t(otu)
    }
    
    # Step 2: Extract sample metadata
    meta <- as.data.frame(sample_data(physeq))
    
    # Step 3: Build SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(counts = otu),
      colData = meta
    )
    
    # Normalize
    se <- relativeAb(se)
    
    # Run LEfSe (assuming your group variable is called 'group')
    res <- lefser(se, class = "group", assay = "rel_abs", lda.threshold = 2.0)
    
    # Plot
    lda <- lefserPlot(res)
    ggsave(file.path(save_path, "lefse_differential_abundance_LDAscore_plot.png"), plot = lda, width = 8, height = 6, dpi = 300)
    
    # Save results
    fwrite(as.data.frame(res), file.path(save_path, "lefse_differential_abundance_results.tsv"), sep = "\t")
    }
  
  message("All results saved to: ", save_path)
}
