#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript step16_trait_split.R <mlm_stats.txt> <output_dir>", call. = FALSE)
}

mlm_path <- args[1]
output_dir <- args[2]

if (!file.exists(mlm_path)) {
  stop(paste("MLM statistics file not found:", mlm_path))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

mlm_stats <- read.table(mlm_path, header = TRUE, sep = "\t", check.names = FALSE)
traits <- unique(mlm_stats$Trait)
cat("Traits found:", paste(traits, collapse = ", "), "\n")

sanitize <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", x)
}

for (trait in traits) {
  trait_data <- subset(mlm_stats, Trait == trait)
  safe_name <- sanitize(trait)
  output_file <- file.path(output_dir, paste0("mlm_stats_", safe_name, ".csv"))
  write.csv(trait_data, output_file, row.names = FALSE, quote = TRUE)
  cat("? Saved:", output_file, "(", nrow(trait_data), "rows)\n")
}
