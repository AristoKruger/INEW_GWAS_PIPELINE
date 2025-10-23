#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    paste(
      "Usage:",
      "Rscript step10_ld_decay_enhanced.R",
      "<ld_summary_path>",
      "<ld_bins_tsv>",
      "<ld_smooth_tsv>",
      "<metrics_tsv>",
      "<plot_png>",
      "<bin_size_bp>",
      "<unused_loess_span>"
    ),
    call. = FALSE
  )
}

ld_summary_path <- args[[1]]
bins_path <- args[[2]]
smooth_path <- args[[3]]
metrics_path <- args[[4]]
plot_path <- args[[5]]
bin_size_bp <- suppressWarnings(as.numeric(args[[6]]))

ld_threshold <- 0.2

if (!file.exists(ld_summary_path)) {
  stop(sprintf("LD summary file not found: %s", ld_summary_path), call. = FALSE)
}
if (!is.finite(bin_size_bp) || bin_size_bp <= 0) {
  stop("bin_size_bp must be a positive numeric value.", call. = FALSE)
}

dir.create(dirname(bins_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(smooth_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)

ld_df <- tryCatch(
  read.delim(
    ld_summary_path,
    sep = "",
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  error = function(e) {
    stop(sprintf("Failed to read LD summary '%s': %s", ld_summary_path, e$message), call. = FALSE)
  }
)

if (ncol(ld_df) < 2) {
  stop("Expected LD summary with at least two columns (distance, r2).", call. = FALSE)
}

ld_df <- ld_df %>%
  transmute(
    dist = as.numeric(V1),
    rsq = as.numeric(V2)
  ) %>%
  filter(is.finite(dist), is.finite(rsq))

if (nrow(ld_df) == 0) {
  stop("No finite distance/r^2 pairs available for LD decay.", call. = FALSE)
}

breaks <- seq(from = max(0, min(ld_df$dist, na.rm = TRUE) - 1),
              to = max(ld_df$dist, na.rm = TRUE) + 1,
              by = bin_size_bp)
if (tail(breaks, n = 1L) < max(ld_df$dist, na.rm = TRUE)) {
  breaks <- c(breaks, max(ld_df$dist, na.rm = TRUE) + bin_size_bp)
}

ld_df <- ld_df %>%
  mutate(dist_bin = cut(dist, breaks = breaks, include.lowest = TRUE))

bin_stats <- ld_df %>%
  group_by(dist_bin) %>%
  summarise(
    mean_rsq = mean(rsq),
    median_rsq = median(rsq),
    .groups = "drop"
  ) %>%
  mutate(
    start = as.numeric(str_extract(dist_bin, "(?<=\\[|\\()(.*?)(?=,)")),
    end = as.numeric(str_extract(dist_bin, "(?<=,)(.*?)(?=\\]|\\))")),
    mid = (start + end) / 2
  ) %>%
  arrange(start)

write.table(
  bin_stats %>%
    select(start_bp = start, end_bp = end, mid_bp = mid, mean_r2 = mean_rsq, median_r2 = median_rsq),
  bins_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  bin_stats %>%
    select(distance_bp = mid, mean_r2 = mean_rsq),
  smooth_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

mean_curve <- bin_stats %>%
  filter(!is.na(mean_rsq)) %>%
  arrange(mid)

threshold_crossing_bp <- NA_real_
if (nrow(mean_curve) >= 2) {
  for (idx in 2:nrow(mean_curve)) {
    y_prev <- mean_curve$mean_rsq[idx - 1L]
    y_curr <- mean_curve$mean_rsq[idx]
    if (!is.finite(y_prev) || !is.finite(y_curr)) {
      next
    }
    if (y_prev > ld_threshold && y_curr <= ld_threshold) {
      x_prev <- mean_curve$mid[idx - 1L]
      x_curr <- mean_curve$mid[idx]
      if (!is.finite(x_prev) || !is.finite(x_curr)) {
        next
      }
      if (y_curr == y_prev) {
        threshold_crossing_bp <- x_curr
      } else {
        threshold_crossing_bp <- x_prev + (ld_threshold - y_prev) * (x_curr - x_prev) / (y_curr - y_prev)
      }
      break
    }
  }
}

metrics <- data.frame(
  metric = c(
    "r2_threshold",
    "bin_size_bp",
    "d0.2_bin_bp",
    "d0.2_loess_bp",
    "total_pairs",
    "total_bins",
    "warnings"
  ),
  value = c(
    ld_threshold,
    bin_size_bp,
    ifelse(is.finite(threshold_crossing_bp), threshold_crossing_bp, NA_real_),
    NA_real_,
    nrow(ld_df),
    nrow(bin_stats),
    "none"
  ),
  stringsAsFactors = FALSE
)

write.table(metrics, metrics_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_df <- bin_stats %>%
  filter(!is.na(mean_rsq)) %>%
  mutate(distance_mb = mid / 1e6)

distance_range <- range(plot_df$distance_mb, na.rm = TRUE)
if (!is.finite(distance_range[1]) || !is.finite(distance_range[2])) {
  distance_range <- c(0, 1)
}
distance_breaks <- pretty(distance_range, n = 6)
if (!0 %in% distance_breaks) {
  distance_breaks <- sort(unique(c(0, distance_breaks)))
}

p <- ggplot(plot_df, aes(x = distance_mb, y = mean_rsq)) +
  geom_point(size = 0.4, colour = "grey20") +
  geom_line(aes(y = mean_rsq), linewidth = 0.8, alpha = 0.5, colour = "grey40") +
  geom_hline(yintercept = ld_threshold, linetype = "dashed", colour = "#cc4c02", linewidth = 0.7) +
  labs(
    x = "Distance (Mb)",
    y = expression(LD~(r^2)),
    title = "LD decay"
  ) +
  scale_x_continuous(
    breaks = distance_breaks,
    limits = range(distance_breaks, na.rm = TRUE)
  ) +
  theme_bw()

if (is.finite(threshold_crossing_bp)) {
  p <- p +
    geom_vline(
      xintercept = threshold_crossing_bp / 1e6,
      colour = "#2ca02c",
      linetype = "dotted",
      linewidth = 0.8
    )
}

ggsave(plot_path, plot = p, width = 8, height = 6, dpi = 300)
