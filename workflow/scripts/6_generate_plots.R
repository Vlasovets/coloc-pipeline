#!/usr/bin/env Rscript
# Stage 6: Generate Locus Zoom Plots
# This script creates visualization of colocalization results

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Get parameters from Snakemake
trait <- snakemake@wildcards[["trait"]]
input_results <- snakemake@input[["results"]]
output_plots <- snakemake@output[["plots"]]
log_file <- snakemake@log[[1]]

# Redirect output to log
sink(log_file, append = FALSE, split = TRUE)
cat("Stage 6: Generate Plots\n")
cat("Trait:", trait, "\n")
cat("Started at:", as.character(Sys.time()), "\n\n")

tryCatch({
  
  # Load aggregated results
  cat("Loading colocalization results...\n")
  results <- fread(input_results)
  cat("Loaded", nrow(results), "results\n\n")
  
  # Filter for high-confidence colocalizations
  sig_results <- results[PP.H4.abf > 0.8]
  cat("Found", nrow(sig_results), "results with PP.H4 > 0.8\n\n")
  
  if (nrow(sig_results) > 0) {
    # Create output directory
    dir.create(output_plots, recursive = TRUE, showWarnings = FALSE)
    
    # TODO: Implement locus zoom plots
    # For now, create summary plot
    cat("Creating summary plots...\n")
    
    # PP.H4 distribution by tissue
    p1 <- ggplot(results, aes(x = tissue, y = PP.H4.abf)) +
      geom_boxplot() +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(title = paste("Colocalization Results -", trait),
           x = "Tissue", y = "PP.H4") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_plots, paste0(trait, "_pp4_distribution.pdf")),
           p1, width = 8, height = 6)
    
    cat("Saved summary plot\n")
    
  } else {
    cat("No significant colocalizations to plot\n")
  }
  
  cat("\nStage 6 completed at:", as.character(Sys.time()), "\n")
  
}, error = function(e) {
  cat("\nERROR in Stage 6:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})

sink()
