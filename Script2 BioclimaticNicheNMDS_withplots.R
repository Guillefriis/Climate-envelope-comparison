setwd('')

load('fullclimate.RData')

library(sf)
library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)
i <- args[1]

n_boot <- 5
r_perm <- 100
d_try  <- 5

## thresholds
min_sources <- 12
min_prop_ok <- 0.80

cols <- c('#009E73', '#D55E00')

print(paste0("Analyzing species: ", i))

climate.list <- split(full_climate, f = full_climate$Species, drop = TRUE)
climate_data <- as.data.frame(climate.list[[i]])
colnames(climate_data) <- colnames(full_climate)

wild_df <- climate_data[climate_data$Group == "W", ]
frm_df  <- climate_data[climate_data$Group == "C", ]

## --------------------------------------------------------------------------------
## PHASE 1: bootstrap only seed source composition
## --------------------------------------------------------------------------------

n_sources <- integer(n_boot)

for (b in seq_len(n_boot)) {
  f_idx <- sample(
    seq_len(nrow(frm_df)),
    size = nrow(frm_df),
    replace = TRUE,
    prob = frm_df$Weight
  )
  n_sources[b] <- length(unique(f_idx))
}

n_low_sources <- sum(n_sources < min_sources)
prop_valid    <- mean(n_sources >= min_sources)

if (prop_valid < min_prop_ok) {

  boot_overlap_summary <- data.frame(
    Species = i,
    Min_sources_threshold = min_sources,
    Prop_valid_bootstraps = prop_valid,
    N_low_source_bootstraps = n_low_sources,
    Passed_threshold = FALSE
  )

  out_file <- paste0("NMDSperm_", gsub(" ", "_", i), ".txt")
  write.table(boot_overlap_summary, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no")
}

## --------------------------------------------------------------------------------
## PHASE 2: NMDS + PERMANOVA
## --------------------------------------------------------------------------------

pdf_file <- paste0("NMDS_bootstrap_", gsub(" ", "_", i), ".pdf")
pdf(pdf_file, width = 10, height = 4)

plot_list <- vector("list", n_boot)

wild_area  <- seed_area  <- overlap_area <- numeric(n_boot)
permanova_R2 <- permanova_p <- disp_p <- numeric(n_boot)
valid_boot <- logical(n_boot)

for (b in seq_len(n_boot)) {

  f_idx <- sample(
    seq_len(nrow(frm_df)),
    size = nrow(frm_df),
    replace = TRUE,
    prob = frm_df$Weight
  )

  if (length(unique(f_idx)) < min_sources) {
    valid_boot[b] <- FALSE
    plot_list[[b]] <- ggplot() + ggtitle(paste("Bootstrap", b, "< 12 sources")) + theme_void()
    next
  }

  valid_boot[b] <- TRUE

  frm_samp <- frm_df[f_idx, ]
  boot_df  <- rbind(frm_samp, wild_df)

  climate_dist <- boot_df %>%
    select(-Species, -Weight, -Group) %>%
    dist(method = "euclidean")

  set.seed(1)
  mds <- metaMDS(
    climate_dist,
    autotransform = FALSE,
    trymax = d_try,
    trace = FALSE,
    k = 2
  )

  group.MDS <- scores(mds, display = "sites") %>% as.data.frame()
  group.MDS$Group <- boot_df$Group

  centroids <- group.MDS %>%
    group_by(Group) %>%
    summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

  h <- ordihull(group.MDS[, 1:2], group.MDS$Group, draw = "none")
  hdf <- do.call(rbind, lapply(names(h), function(g)
    transform(as.data.frame(h[[g]]), Group = g)))
  colnames(hdf)[1:2] <- c("NMDS1", "NMDS2")

  plot_list[[b]] <- ggplot(group.MDS, aes(NMDS1, NMDS2, colour = Group)) +
    geom_path(data = hdf, aes(NMDS1, NMDS2, colour = Group), inherit.aes = FALSE) +
    geom_point(size = 1.5) +
    geom_point(
      data = centroids,
      aes(NMDS1, NMDS2, fill = Group),
      size = 4,
      shape = 21,
      colour = "black",
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = cols,
      breaks = c("C", "W"),
      labels = c("Seed sources", "Wild distribution"),
      name = ""
    ) +
    scale_fill_manual(
      values = cols,
      breaks = c("C", "W"),
      labels = c("Seed sources", "Wild distribution"),
      name = ""
    ) +
    ggtitle(paste(i, "- bootstrap", b)) +
    theme(
      plot.title = element_text(face = "italic"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "none"
    )

  ## PERMANOVA
  perm <- vegan::adonis2(climate_dist ~ Group, data = group.MDS, permutations = r_perm)
  permanova_R2[b] <- perm$R2[1]
  permanova_p[b]  <- perm[["Pr(>F)"]][1]

  ## Dispersion
  bd <- vegan::betadisper(climate_dist, group.MDS$Group)
  disp_p[b] <- anova(bd)$`Pr(>F)`[1]

  ## Convex hull areas
  wpts <- h[["W"]]
  cpts <- h[["C"]]

  if (!is.null(wpts) && !is.null(cpts) && nrow(wpts) >= 3 && nrow(cpts) >= 3) {
    w_poly <- st_sfc(st_polygon(list(rbind(wpts, wpts[1, ]))))
    c_poly <- st_sfc(st_polygon(list(rbind(cpts, cpts[1, ]))))
    wild_area[b] <- st_area(w_poly)
    seed_area[b] <- st_area(c_poly)
    inter <- st_intersection(w_poly, c_poly)
    overlap_area[b] <- if (length(inter) == 0) NA_real_ else st_area(inter)
  }
}

gridExtra::grid.arrange(grobs = plot_list, nrow = 1)
dev.off()

## --------------------------------------------------------------------------------
## SUMMARY STATISTICS
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
  Mean_percent_represented = mean_pr,
  Median_PERMANOVA_p = median_perm_p,
  Median_Dispersion_p = median_disp_p
)

out_file <- paste0("NMDSperm_", gsub(" ", "_", i), ".txt")
write.table(boot_overlap_summary, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
