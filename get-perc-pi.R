#!/usr/bin/env Rscript

# handling args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
   stop("Usage: <input_file.pi> <target_row_index> <target_region_windows>")
}

file_path <- args[1]
target_row <- as.numeric(args[2])
region_windows <- as.numeric(args[3])

# reading data
df <- read.table(file_path, header = TRUE)
pi_vals <- df[[5]]

row_pi <- pi_vals[target_row]
pi_vals <- pi_vals[!is.na(pi_vals) & is.finite(pi_vals)]

# getting the fraction of diveristy values which are larger than our row.
p_raw <- mean(pi_vals >= row_pi, na.rm = TRUE)
# correcting for number of windows in our region.
p_corrected <- min(p_raw * region_windows, 1)

# output
cat(sprintf("PI value: %.6f\n", row_pi))
cat(sprintf("Percentile: %.2f%%\n", mean(pi_vals <= row_pi, na.rm = TRUE) * 100))
cat(sprintf("p (raw): %.4f\n", p_raw))
cat(sprintf("p (Bonferroni): %.4f\n", p_corrected))

if (abs(p_corrected) < 0.05) {
  cat("Result: Statistically significant (p < 0.05)\n")
} else {
  cat("Result: Not significant at alpha 0.05\n")
}
