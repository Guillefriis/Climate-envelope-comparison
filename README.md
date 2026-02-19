# Bioclimatic niche comparison: wild distribution vs seed sources

R scripts to compare the climatic niche represented by registered seed source stands against the climatic niche occupied by wild populations for each species. The aim is to test whether seed sources underrepresent wild climatic conditions and to visualise and summarise any bias.

Project context: Friis et al. 2026 (in preparation)

## Scripts and run order

### 1) Data processing
**`Script1 DataProcessing.R`**  
Prepares the core datasets used by the downstream analyses, including:
- wild occurrence or distribution inputs
- seed source stand inputs
- extraction or assembly of bioclimatic variables
- creation of matched tables used by NMDS and summaries

### 2) NMDS niche ordination and plots
**`Script2 BioclimaticNicheNMDS_withplots.R`**  
Runs NMDS ordination on bioclimatic variables to represent niche space and produces ordination plots comparing wild vs seed sources.

### 3) Probabilistic bootstrap on niche representation
**`Script3a BioclimaticNicheNMDS_probbootstrap.R`**  
Implements bootstrap based inference to quantify uncertainty in niche representation and underrepresentation metrics (wild vs seed sources) in ordination space.

### 4) Seed source metrics
**`Script3b BioclimaticNicheNMDS_sourcemetrics.R`**  
Computes summary metrics describing how seed sources cover wild climatic space (coverage, distances, dispersion, and related summaries used in reporting) based on bootstrap replicates.

### 5) Regression and diagnostic plots
**`Script4 RegressionPlots.R`**  
Generates regression style plots and related diagnostics used to interpret relationships between niche coverage metrics and explanatory variables.

## Inputs and outputs

Inputs are expected to include:
- bioclimatic predictor layers or extracted values
- wild occurrence or distribution data
- seed source stand data

Outputs include:
- NMDS objects and ordination figures
- bootstrap summaries
- seed source niche representation metrics tables
- regression figures used in reporting

## Requirements

R ≥ 4.2 recommended. Common packages for this workflow typically include vegan, tidyverse, and plotting and spatial packages depending on how climate variables are handled locally.

## License

MIT License (see `LICENSE`).

