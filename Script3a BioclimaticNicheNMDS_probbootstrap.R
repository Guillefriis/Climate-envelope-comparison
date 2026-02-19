setwd('')

library(sf)
library(tidyverse)
library(vegan)

## Bootstrapped NMDS + PERMANOVA
#-----------------------------------------------------------------------------------------------------------
cat("Starting NMDS + PERMANOVA analysis...\n")

load('fullclimate.RData')

args <- commandArgs(trailingOnly = TRUE)
i <- args[1]  # Species name passed from the bash script

n_boot <- 1000        # Number of bootstraps
r_perm <- 1000        # Number of permutations in PERMANOVA
d_try  <- 10

## thresholds
min_sources <- 12
min_prop_ok <- 0.80

print(paste0("Analyzing species: ", i))

climate.list <- split(full_climate, f = full_climate$Species, drop = TRUE)
climate_data <- as.data.frame(climate.list[[i]])
colnames(climate_data) <- colnames(full_climate)

wild_df <- climate_data[climate_data$Group == "W", ]
frm_df  <- climate_data[climate_data$Group == "C", ]

## --------------------------------------------------------------------------------
## PHASE 1: bootstrap only seed source composition
## --------------------------------------------------------------------------------

cat("Phase 1: bootstrapping seed source composition...\n")

n_sources <- integer(n_boot)

for (b in seq_len(n_boot)) {

  f_idx <- sample(
    seq_len(nrow(frm_df)),
    size     = nrow(frm_df),
    replace  = TRUE,
    prob     = frm_df$Weight
  )

  n_sources[b] <- length(unique(f_idx))
}

n_low_sources <- sum(n_sources < min_sources)
prop_valid    <- mean(n_sources >= min_sources)

cat("Bootstraps with < ", min_sources, " sources: ", n_low_sources, "\n", sep = "")
cat("Proportion valid bootstraps: ", round(prop_valid, 3), "\n", sep = "")

## If species fails threshold, stop here and write summary
if (prop_valid < min_prop_ok) {

  cat("Species failed 80% rule. Skipping NMDS + PERMANOVA.\n")

  boot_overlap_summary <- data.frame(
    Species = i,
    Min_sources_threshold = min_sources,
    Prop_valid_bootstraps = prop_valid,
    N_low_source_bootstraps = n_low_sources,
    Passed_threshold = FALSE,

    Mean_wild_area = NA,
    Mean_stand_area = NA,
    Mean_overlap = NA,

    Mean_percent_represented = NA,
    SD_percent_represented   = NA,
    SE_percent_represented   = NA,
    CI_repr_low  = NA,
    CI_repr_high = NA,

    Mean_percent_outside = NA,
    Median_PERMANOVA_p = NA,
    Rejection_rate_PERMANOVA = NA,
    Median_PERMANOVA_R2 = NA,
    CI_PERMANOVA_R2_low = NA,
    CI_PERMANOVA_R2_high = NA,
    Median_Dispersion_p = NA
  )

  out_file <- paste0("NMDSperm_", gsub(" ", "_", i), ".txt")

  write.table(
    boot_overlap_summary,
    out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Done.\n")
  quit(save = "no")
}

## --------------------------------------------------------------------------------
## PHASE 2: NMDS + PERMANOVA (valid species only)
## --------------------------------------------------------------------------------

cat("Phase 2: running NMDS + PERMANOVA...\n")

wild_area  <- seed_area  <- overlap_area <- numeric(n_boot)
permanova_R2 <- permanova_p <- disp_p <- numeric(n_boot)
valid_boot <- logical(n_boot)

for (b in seq_len(n_boot)) {

  f_idx <- sample(
    seq_len(nrow(frm_df)),
    size     = nrow(frm_df),
    replace  = TRUE,
    prob     = frm_df$Weight
  )

  if (length(unique(f_idx)) < min_sources) {
    valid_boot[b] <- FALSE
    wild_area[b]  <- NA
    seed_area[b]  <- NA
    overlap_area[b] <- NA
    permanova_p[b]  <- NA
    permanova_R2[b] <- NA
    disp_p[b]       <- NA
    next
  }

  valid_boot[b] <- TRUE

  frm_samp <- frm_df[f_idx, ]
  boot_df  <- rbind(frm_samp, wild_df)

  climate_dist <- boot_df %>%
    select(-Species, -Weight, -Group) %>%
    dist(method = "euclidean")

  set.seed(1)
  boot_mat <- climate_dist %>%
    metaMDS(
      autotransform = FALSE,
      distance = "euclidean",
      trymax = d_try,
      trace = FALSE,
      k = 2
    ) %>%
    scores(display = "sites") %>%
    as.data.frame()

  boot_mat$Group <- boot_df$Group

  ## PERMANOVA
  perm <- vegan::adonis2(climate_dist ~ Group, data = boot_mat, permutations = r_perm)
  permanova_R2[b] <- perm$R2[1]
  permanova_p[b]  <- perm[["Pr(>F)"]][1]

  ## Dispersion
  bd <- vegan::betadisper(climate_dist, boot_mat$Group)
  disp_p[b] <- anova(bd)$`Pr(>F)`[1]

  ## Convex hull
  hull_obj <- ordihull(boot_mat[, 1:2], boot_mat$Group, draw = "none")

  wpoly_out2 <- rbind(hull_obj[["W"]], hull_obj[["W"]][1, , drop = FALSE])
  cpoly_out2 <- rbind(hull_obj[["C"]], hull_obj[["C"]][1, , drop = FALSE])

  if (nrow(wpoly_out2) >= 3 && nrow(cpoly_out2) >= 3) {

    w_poly <- st_sfc(st_polygon(list(as.matrix(wpoly_out2))))
    c_poly <- st_sfc(st_polygon(list(as.matrix(cpoly_out2))))

    wild_area[b] <- st_area(w_poly)
    seed_area[b] <- st_area(c_poly)

    inter <- tryCatch(
      st_intersection(w_poly, c_poly),
      error = function(e) NA_real_
    )

    overlap_area[b] <- if (is.na(inter) || length(inter) == 0)
      NA_real_ else st_area(inter)

  } else {
    wild_area[b]    <- NA
    seed_area[b]    <- NA
    overlap_area[b] <- NA
  }
}

## --------------------------------------------------------------------------------
## SUMMARY STATISTICS (valid bootstraps only)
## --------------------------------------------------------------------------------

percent_represented <- 100 * overlap_area / wild_area
percent_outside     <- 100 * (seed_area - overlap_area) / seed_area

mean_pr  <- mean(percent_represented[valid_boot], na.rm = TRUE)
sd_pr    <- sd(percent_represented[valid_boot], na.rm = TRUE)
se_pr    <- sd_pr / sqrt(sum(valid_boot))
ci_pr_lo <- quantile(percent_represented[valid_boot], 0.025, na.rm = TRUE)
ci_pr_hi <- quantile(percent_represented[valid_boot], 0.975, na.rm = TRUE)

median_perm_p   <- median(permanova_p[valid_boot], na.rm = TRUE)
reject_rate_perm <- mean(permanova_p[valid_boot] < 0.05, na.rm = TRUE)
median_perm_R2  <- median(permanova_R2[valid_boot], na.rm = TRUE)
CI_R2_low       <- quantile(permanova_R2[valid_boot], 0.025, na.rm = TRUE)
CI_R2_high      <- quantile(permanova_R2[valid_boot], 0.975, na.rm = TRUE)
median_disp_p   <- median(disp_p[valid_boot], na.rm = TRUE)

boot_overlap_summary <- data.frame(
  Species = i,
  Min_sources_threshold = min_sources,
  Prop_valid_bootstraps = prop_valid,
  N_low_source_bootstraps = n_low_sources,
  Passed_threshold = TRUE,

  Mean_wild_area  = mean(wild_area[valid_boot], na.rm = TRUE),
  Mean_stand_area = mean(seed_area[valid_boot], na.rm = TRUE),
  Mean_overlap    = mean(overlap_area[valid_boot], na.rm = TRUE),

  Mean_percent_represented = mean_pr,
  SD_percent_represented   = sd_pr,
  SE_percent_represented   = se_pr,
  CI_repr_low  = ci_pr_lo,
  CI_repr_high = ci_pr_hi,

  Mean_percent_outside = mean(percent_outside[valid_boot], na.rm = TRUE),
  Median_PERMANOVA_p = median_perm_p,
  Rejection_rate_PERMANOVA = reject_rate_perm,
  Median_PERMANOVA_R2 = median_perm_R2,
  CI_PERMANOVA_R2_low = CI_R2_low,
  CI_PERMANOVA_R2_high = CI_R2_high,
  Median_Dispersion_p = median_disp_p
)

out_file <- paste0("NMDSperm_", gsub(" ", "_", i), ".txt")

write.table(
  boot_overlap_summary,
  out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")




























