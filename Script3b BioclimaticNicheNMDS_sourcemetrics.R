setwd('')

library(tidyverse)

cat("Starting bootstrap source metrics...\n")

load('fullclimate.RData')

args <- commandArgs(trailingOnly = TRUE)
i <- args[1]

n_boot <- 1000
min_sources <- 12

print(paste0("Analyzing species: ", i))

climate.list <- split(full_climate, f = full_climate$Species, drop = TRUE)
climate_data <- as.data.frame(climate.list[[i]])
colnames(climate_data) <- colnames(full_climate)

frm_df <- climate_data[climate_data$Group == "C", ]

n_sources <- integer(n_boot)

for (b in seq_len(n_boot)) {

  f_idx <- sample(
    seq_len(nrow(frm_df)),
    size    = nrow(frm_df),
    replace = TRUE,
    prob    = frm_df$Weight
  )

  n_sources[b] <- length(unique(f_idx))
}

n_low_sources <- sum(n_sources < min_sources)
prop_valid    <- mean(n_sources >= min_sources)

out <- data.frame(
  Species = i,
  N_boot = n_boot,
  Min_sources_threshold = min_sources,
  Mean_unique_sources = mean(n_sources),
  Median_unique_sources = median(n_sources),
  Q2.5_unique_sources = quantile(n_sources, 0.025),
  Q97.5_unique_sources = quantile(n_sources, 0.975),
  N_low_source_bootstraps = n_low_sources,
  Prop_valid_bootstraps = prop_valid
)

out_file <- paste0("SourceMetrics_", gsub(" ", "_", i), ".txt")

write.table(out,
            out_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
