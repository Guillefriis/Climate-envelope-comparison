setwd('')

library(tidyverse)
library(readxl)
library(dplyr)
library(vegan)
library(janitor)
library(broom)
library(scales)
library(ggplot2)

# Load PERMANOVA results
##--------------------------------------------------------------------------------------------------------
xls_path <- "BioclimaticNiche_Tables.xlsx"
bcn_table <- read_excel(xls_path, sheet = 'WeightedPERM') %>% janitor::clean_names()

names(bcn_table)

cover_df <- bcn_table %>%
  select(Species = species,
         Sources = seed_sources_count,
         Mean_percent_represented = mean_percent_represented,
         )
         
# FRM SpeciesClimate table
##--------------------------------------------------------------------------------------------------------
frm_table <- read_excel("FRM_SpeciesClimate.xlsx")

frm_quantity <- frm_table %>% 
  select(Species = Species, Quantity = Quantity)

frm_quantity <- frm_quantity %>%
  group_by(Species) %>%
  filter(n() > 5) %>%
  ungroup()

seed_metrics <- frm_quantity %>%
  group_by(Species) %>%
  mutate(p = Quantity / sum(Quantity, na.rm = TRUE)) %>%
  summarise(
    n_sources   = n(),
    total_qty   = sum(Quantity, na.rm = TRUE),
    H_shannon   = -sum(p * log(p), na.rm = TRUE),
    eff_sources = exp(H_shannon),
    evenness    = ifelse(n_sources > 1, H_shannon / log(n_sources), NA_real_)
  ) %>%
  ungroup()


# Correlations wit rejection rate
##--------------------------------------------------------------------------------------------------------
## Join + IDs and highlight flag
analysis_df <- cover_df %>%
  left_join(seed_metrics, by = "Species") %>%
  arrange(Species) %>%
  mutate(
    ID = row_number(),
  )

writexl::write_xlsx(analysis_df, path = "WeightedShannon_Table.xlsx")

species_legend <- analysis_df %>%
  arrange(ID) %>%
  transmute(line = paste0(ID, " = ", Species)) %>%
  pull(line) %>%
  paste(collapse = "\n")

## Regression test
# Plot: effective sources vs coverage
df2 <- analysis_df %>%
  select(x = eff_sources, y = Mean_percent_represented, ID) %>%
  tidyr::drop_na()

m2  <- lm(y ~ x, data = df2)
r2_2 <- summary(m2)$r.squared
p2   <- summary(m2)$coefficients["x", "Pr(>|t|)"]
lab2 <- paste0("R² = ", round(r2_2, 3), ",  p = ", formatC(p2, format = "g", digits = 3))

fig2 <- ggplot(df2, aes(x, y)) +
  geom_point() +
  geom_text(aes(label = ID), vjust = -0.6, size = 3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, colour = "blue") +
  scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "red3"),
                     name = "PERMANOVA\nRejection > 0.95") +
  labs(x = "Effective number of seed sources (exp(H'))",
       y = "Mean % wild climate represented") +
  # stats (top-right)
  annotate("text",
           x = max(df2$x, na.rm = TRUE),
           y = max(df2$y, na.rm = TRUE),
           label = lab2, hjust = 1, vjust = 1) +
  # species legend box (top-left)
  annotate("label",
           x = min(df2$x, na.rm = TRUE),
           y = max(df2$y, na.rm = TRUE),
           label = species_legend, hjust = 0, vjust = 1,
           size = 3, label.size = 0.25)

print(fig2)

ggsave("EffSources_Coverage.pdf", plot = fig2, width = 8, height = 6, units = "in")
