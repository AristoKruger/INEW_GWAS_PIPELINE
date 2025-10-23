#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript step17_trait_adjust.R <input_dir> <pattern> <output_dir> <summary_file> <n_tests> <top_n>", call. = FALSE)
}

input_dir <- args[1]
pattern <- args[2]
output_dir <- args[3]
summary_file <- args[4]
n_tests <- as.numeric(args[5])
top_n <- as.integer(args[6])

if (!dir.exists(input_dir)) {
  stop(paste("Input directory does not exist:", input_dir))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
summary_dir <- dirname(summary_file)
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
}

bonf_threshold <- 0.05 / n_tests
bonf_threshold_log10 <- -log10(bonf_threshold)

trait_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(trait_files) == 0) {
  stop(paste("No trait files found in", input_dir, "matching pattern", pattern))
}

signif_summary <- data.frame()
cat("Running correction using", n_tests, "tests per trait\n")
cat("Bonferroni threshold =", format(bonf_threshold, scientific = TRUE), "\n")

for (file in trait_files) {
  trait_name <- sub("mlm_stats_(.*)\\.csv", "\\1", basename(file))
  data <- read.csv(file, check.names = FALSE)
  data$p <- as.numeric(data$p)

  clean <- data %>% filter(!is.na(p) & is.finite(p) & p > 0)

  adj <- clean %>%
    mutate(
      p_Bonferroni = pmin(p * n_tests, 1),
      p_FDR = p.adjust(p, method = "fdr"),
      log10p = -log10(p),
      log10_threshold = bonf_threshold_log10,
      significant = ifelse(p < bonf_threshold, "yes", "no")
    ) %>%
    select(Trait, Marker, Chr, Pos, p, p_Bonferroni, p_FDR, log10p, log10_threshold, significant)

  adj_file <- file.path(output_dir, paste0("adj_p_", trait_name, ".csv"))
  write.csv(adj, adj_file, row.names = FALSE, quote = TRUE)
  cat("? Saved:", adj_file, "(", nrow(adj), "rows)\n")

  signif_summary <- bind_rows(signif_summary, adj %>% filter(significant == "yes"))

  top_file <- file.path(output_dir, paste0("top_", top_n, "_", trait_name, ".csv"))
  top_snps <- adj %>% arrange(p) %>% head(top_n)
  write.csv(top_snps, top_file, row.names = FALSE, quote = TRUE)
}

write.csv(signif_summary, summary_file, row.names = FALSE, quote = TRUE)
cat("Summary saved:", summary_file, "with", nrow(signif_summary), "significant SNPs total\n")
