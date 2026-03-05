suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# Source modular functions with absolute paths
script_dir <- "/home/itg/oleg.vlasovets/projects/coloc-pipeline/workflow/scripts"
cat("Sourcing coloc_helpers.R...\n")
source(file.path(script_dir, "coloc_helpers.R"))
cat("Sourcing data_loader.R...\n")
source(file.path(script_dir, "data_loader.R"))
cat("Sourcing qtl_processor.R...\n")
source(file.path(script_dir, "qtl_processor.R"))
cat("All sourced successfully!\n")
