#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript step13_pcoa.R <dissimilarity_csv> <labels_txt> <coordinates_csv> <plot_png>", call. = FALSE)
}

dissimilarity_path <- args[1]
labels_path <- args[2]
coordinates_path <- args[3]
plot_path <- args[4]

if (!file.exists(dissimilarity_path)) {
  stop(paste("Dissimilarity matrix not found:", dissimilarity_path))
}
if (!file.exists(labels_path)) {
  stop(paste("Label file not found:", labels_path))
}

make_dir <- function(path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}

make_dir(coordinates_path)
make_dir(plot_path)

raw <- read.csv(dissimilarity_path, check.names = FALSE)
if (ncol(raw) < 2) {
  stop("Dissimilarity matrix must have at least two columns")
}

if (!is.numeric(raw[[1]])) {
  row_labels <- raw[[1]]
  raw <- raw[-1]
} else {
  row_labels <- NULL
}

matrix_data <- as.matrix(raw)
labels <- readLines(labels_path)
if (length(labels) != nrow(matrix_data)) {
  warning("Label count does not match matrix rows; using available labels")
  labels <- labels[seq_len(min(length(labels), nrow(matrix_data)))]
}

fit <- cmdscale(matrix_data, eig = TRUE, k = min(5, nrow(matrix_data) - 1))
coords <- as.data.frame(fit$points)
rownames(coords) <- labels

write.csv(coords, file = coordinates_path, row.names = TRUE)

png(plot_path, width = 800, height = 800)
plot(coords[, 1], coords[, 2], xlab = "PCo 1", ylab = "PCo 2", main = "Principal Coordinate Analysis", pch = 19)
text(coords[, 1], coords[, 2], labels = labels, cex = 0.7, pos = 3)
dev.off()
