setwd('')

library(readxl)
library(stringr)
library(dplyr)
library(tidyverse)
library(writexl)
library(scales)

#-----------------------------------------------------------------------------------------------------------
# Process raw data
##------------------------------------------------------------------------------------------------------
# Read and tidy
frm.df <- read_xlsx(
  path = "MyData.xlsx",
  sheet = "Analysis", guess_max = 3000
)

# Rename columns 5 and 9, then remove spaces from all column names
names(frm.df)[c(5, 9)] <- c("Location", "Unit")
names(frm.df) <- gsub(" ", "", names(frm.df))

# Coerce Quantity to numeric early
frm.df <- frm.df %>%
  mutate(Quantity = as.numeric(Quantity))

# Create Site identifier
frm.df <- frm.df %>%
  mutate(
    Genus   = sub(" .*", "", Species),
    Epithet = sub(".* ", "", Species),
    Prefix  = paste0(substr(Genus, 1, 1), substr(Epithet, 1, 3)),
    Site    = paste0(Prefix, ave(seq_along(Species), Species, FUN = seq_along))
  ) %>%
  dplyr::select(-Genus, -Epithet, -Prefix)

# Extract UK OS grid refs into a list column
os_pattern <- "\\b([A-Z]{2})\\s*(?:\\d{3}\\s?\\d{3}|\\d{4}\\s?\\d{4}|\\d{5}\\s?\\d{5})\\b"

frm.df <- frm.df %>%
  mutate(
    Location_clean = toupper(gsub("[\r\n]+", " ", Location)),
    OSGrid_List = str_extract_all(Location_clean, os_pattern),
    OSGrid_List = lapply(OSGrid_List, function(v) unique(gsub("\\s+", "", v)))
  ) %>%
  dplyr::select(-Location_clean)

# One row per grid ref, divide Quantity by refs-per-row
frm.long <- frm.df %>%
  mutate(n_refs_row = lengths(OSGrid_List)) %>%
  filter(n_refs_row > 0) %>%
  mutate(Quantity = Quantity / n_refs_row) %>%
  unnest_longer(OSGrid_List, values_to = "OSGrid_Ref") %>%
  dplyr::select(-n_refs_row) %>%
  # Pop: unique per row based on the alpha prefix of Site
  mutate(Site_base = sub("\\d+$", "", Site)) %>%
  group_by(Site_base) %>%
  mutate(Pop = paste0(Site_base, row_number())) %>%
  ungroup() %>%
  dplyr::select(-Site_base)

sum(is.na(frm.long$OSGrid_Ref))

# Write both versions
write_xlsx(
  list(Full = frm.df, ByGridRef = frm.long),
  path = "FRMdata_FULL.xlsx"
)

rm(frm.df, frm.long, os_pattern)

#-----------------------------------------------------------------------------------------------------------
# Climatic datasets
##------------------------------------------------------------------------------------------------------

library(sf)
library(terra)
library(raster)
library(rgdal)

## Loading in climate data
#-----------------------------------------------------------------------------------------------------------
climpath <- 'Climate'
nc_files <- list.files(path = climpath, pattern = "\\.nc$", full.names = T)
nc_files <- c(nc_files, 'Climate/wc2.1_30s_elev.tif')

nc_names <- c("days_frost",
              "mean_humidity",
              "sealev_pressure",
              "vapour_press",
              "total_rain",
              "wind_speed",
              "days_snow",
              "sun_dura",
              "mean_temp",
              "max_temp",
              "min_temp",
              "elevation")

for (i in seq_along(nc_files)) {
  file_path <- nc_files[i]
  var_name <- nc_names[i]
  assign(var_name, raster(file_path))
}

# Stack climate data
climate <- stack(total_rain, mean_temp, max_temp, min_temp, sun_dura, wind_speed,
                 sealev_pressure, mean_humidity, vapour_press, days_frost, days_snow)

rasterblank <- projectExtent(climate, elevation)
climate_new <- projectRaster(climate, rasterblank)
inbt1 <- crop(elevation, extent(climate_new))
inbt2 <- resample(inbt1, climate_new, method= "ngb")
climate_stack <- stack(climate_new, inbt2)
plot(climate_stack)

# Import seed zone map
sz_outline <- readOGR(dsn = "SeedZone", layer = "Forest_Reproductive_Materials_Regions_Of_Provenance_GB")
sz_shape <- spTransform(sz_outline, CRS(proj4string(climate_stack$wc2.1_30s_elev)))
plot(sz_shape)

# Rasterize seed zones
ext <- extent(sz_shape)
r <- raster(ext, res = 0.1)
sz_raster <- rasterize(sz_shape, r, field = sz_shape$SEED_ZONES)
crs(sz_raster) <- crs(climate_stack$wc2.1_30s_elev)

# Reproject climate data to match seed zone raster
sz_new <- projectRaster(climate_stack, sz_raster)

# Crop and resample seed zone raster to match climate data
sz_raster_cropped <- crop(sz_raster, extent(sz_new))
sz_raster_resampled <- resample(sz_raster_cropped, sz_new, method = "ngb")

# Stack reprojected climate data and resampled seed zone raster
sz_stack <- stack(sz_new, sz_raster_resampled)
plot(sz_stack)

#Rename layers
names(sz_stack)

names(sz_stack) <- c("Total_precip",
                     "Mean_air_temp",
                     "Max_air_temp",
                     "Min_air_temp",
                     "Sunshine_hours",
                     "Wind_speed_10m",
                     "Pressure_sealev",
                     "Rel_humidity",
                     "Watervap_pressure",
                     "Days_ground_frost",
                     "Days_snow_9am",
                     "Altitude",
                     "Seed_zone")

rm(total_rain, mean_temp, max_temp, min_temp, sun_dura, wind_speed,sealev_pressure, climate,
   mean_humidity, vapour_press, days_frost, days_snow, inbt1, inbt2, climpath, climate_new,
   rasterblank, file_path, i, nc_files, nc_names, var_name, ext, r, sz_new, sz_raster, sz_raster_cropped,
   sz_raster_resampled)


## FRM data table (Seed Sources)
#-----------------------------------------------------------------------------------------------------------
frm_table <- read_excel("FRMdata_FULL.xlsx", sheet = "ByGridRef")

frm_location <- frm_table %>% 
  dplyr::select(Species_bin = Species, lat = Lat, lon = Lon,
                Seed_zone = SeedZone, Grid_reference = OSGrid_Ref,
                Category = Category, Quantity = Quantity)

Batch1_species <- unique(frm_location$Species_bin)
writeLines(sort(Batch1_species), "taxa_list.txt")
species_location <- frm_location

## Heat-map map
sz_stread <- st_read(dsn = "SeedZone", layer = "Forest_Reproductive_Materials_Regions_Of_Provenance_GB")

seed_lookup <- tribble(
  ~Species_bin,            ~seeds_per_kg,
  "Acer campestre",         18000,
  "Acer pseudoplatanus",     8500,
  "Alnus glutinosa",       192000,
  "Aria edulis",            53200,
  "Betula pendula",       1500000,
  "Betula pubescens",     4800000,
  "Carpinus betulus",       24200,
  "Castanea sativa",          239,
  "Cornus sanguinea",       20200,
  "Corylus avellana",         797,
  "Crataegus laevigata",    11200,
  "Crataegus monogyna",     11200,
  "Euonymus europaeus",     29000,
  "Fagus sylvatica",         4600,
  "Frangula alnus",         57600,
  "Ilex aquifolium",        35800,
  "Juniperus communis",     80300,
  "Larix decidua",          95000,
  "Malus sylvestris",       61000,
  "Picea abies",           136000,
  "Picea sitchensis",      230000,
  "Pinus sylvestris",       45000,
  "Populus spp.",         8500000,
  "Prunus avium",            5100,
  "Prunus padus",           20100,
  "Prunus spinosa",          5400,
  "Quercus petraea",          316,
  "Quercus robur",            273,
  "Salix aurita",            7000,
  "Salix caprea",            7000,
  "Salix cinerea",           7000,
  "Salix fragilis",          7040,
  "Salix lanata",            7000,
  "Salix lapponum",          7000,
  "Salix myrsinifolia",      7000,
  "Salix myrsinites",        7000,
  "Salix pentandra",         5720,
  "Salix phylicifolia",      7000,
  "Salix purpurea",          7000,
  "Salix repens",            7000,
  "Salix viminalis",         7000,
  "Sambucus nigra",        350000,
  "Sorbus aucuparia",      290000,
  "Sorbus torminalis",      45300,
  "Taxus baccata",          17000,
  "Tilia cordata",          30100,
  "Ulmus glabra",           89300
)

frm_with_seeds <- frm_location %>%
  left_join(seed_lookup, by = "Species_bin") %>%
  mutate(Seeds_est = Quantity * seeds_per_kg)

zone_seeds <- frm_with_seeds %>%
  group_by(Seed_zone) %>%
  summarise(Seeds_est = sum(Seeds_est, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct = Seeds_est / sum(Seeds_est, na.rm = TRUE))

# join polygons
sz_stread <- sz_stread %>%
  left_join(zone_seeds, by = c("SEED_ZONES" = "Seed_zone")) %>%
  mutate(.zone_id = SEED_ZONES)

# labels: "<zone>\n<percent>"
labels_sf <- st_point_on_surface(sz_stread) %>%
  mutate(lbl = paste0(.zone_id, "\n",
                      label_percent(accuracy = 0.1)(pct)))

# plot
ggplot() +
  geom_sf(data = sz_stread, aes(fill = Seeds_est)) +
  scale_fill_gradient(
    low = "yellow", high = "red", na.value = "grey90",
    trans = "log10",
    labels = label_number(scale_cut = cut_si(""))
  ) +
  geom_sf_text(data = labels_sf, aes(label = lbl),
               size = 2.2, lineheight = 0.9, fontface = "bold") +
  theme_minimal() +
  labs(
    title = "Seed Zone Collection Heatmap",
    fill = "Estimated seeds"
  )

write_xlsx(as.data.frame(sz_stread[, 2:6, drop = TRUE]),
                        path = "Heatmap_Table.xlsx")

rm(sz_stread, seed_lookup, frm_with_seeds, zone_seeds)

## Plots
#-----------------------------------------------------------------------------------------------------------
cols <- rainbow(length(Batch1_species))
layout_matrix <- matrix(c(1, 4,
                          2, 5,
                          3, 6), nrow = 3, ncol = 2, byrow = TRUE)

pdf("FRM_seed_sources_maps.pdf", width = 10, height = 14)

# Loop over each species
for (i in 1:length(Batch1_species)) {
  
  # start a new page every 6 species
  if ((i - 1) %% 6 == 0) {
    layout(mat = layout_matrix, heights = c(1, 1, 1), widths = c(1, 1))
    par(mar = c(1, 1, 1, 1))
  }
  
  species <- Batch1_species[i]
  species_split <- species_location[grep(species, species_location$Species_bin), ]
  species_grid <- subset(species_split, select = c(Species_bin, lat, lon))
  species_grid <- species_grid[complete.cases(species_grid$lat, species_grid$lon), ]  # avoid NA coords
  
  plot(sz_shape)
  if (nrow(species_grid) > 0) {
    points(species_grid$lon, species_grid$lat, cex = 1.5, pch = 21, bg = cols[i])
  }
  title(main = bquote(italic(.(species))), line = -0.4, adj = 0.45)
  
}

dev.off()

# Extract climate data for seed sources
sp_grid <- subset(species_location, select = c(Species_bin, lat, lon, Quantity))
sp_climate <- terra::extract(sz_stack, sp_grid[ , c(3, 2)])
sp_climate <- as.data.frame(cbind(sp_climate, sp_grid[, c('Quantity', 'Species_bin')]))

colnames(sp_climate)[ncol(sp_climate)] <- "Species"
sp_climate <- sp_climate[complete.cases(sp_climate), ]

## Tidy up
species_climate <- subset(sp_climate, select = -c(Altitude, Seed_zone)) # We have used the Elevation data for compiling layers
                                                                        # but we get rid of it here because is not climate-informative

species_climate <- species_climate %>% mutate_at(vars(1:12), as.numeric)
table(species_climate$Species)

## Clump identical sites and filter species with less than 3 sources, IMPORTANT FOR NMDS spaces computing! 
clim_vars <- c("Total_precip","Mean_air_temp","Max_air_temp","Min_air_temp",
               "Sunshine_hours","Wind_speed_10m","Pressure_sealev","Rel_humidity",
               "Watervap_pressure","Days_ground_frost","Days_snow_9am")

species_climate_pre <- species_climate

summary_C_pre <- species_climate_pre %>%
  group_by(Species) %>%
  summarise(
    n_total            = n(),
    n_distinct_C_sites = nrow(dplyr::distinct(pick(all_of(clim_vars)))),
    .groups = "drop"
  )

pre_lt3_total <- species_climate_pre %>%
  semi_join(summary_C_pre %>% filter(n_total < 3), by = "Species")

pre_lt3_distinct <- species_climate_pre %>%
  semi_join(summary_C_pre %>% filter(n_distinct_C_sites < 3), by = "Species")

writexl::write_xlsx(
  list(
    pre_lt3_total    = pre_lt3_total,
    pre_lt3_distinct = pre_lt3_distinct,
    pre_summary      = summary_C_pre
  ),
  "species_C_preclump_review.xlsx"
)

# Clumping
species_climate <- species_climate_pre %>%
  group_by(Species, across(all_of(clim_vars))) %>%
  summarise(Quantity = sum(Quantity, na.rm = TRUE), .groups = "drop")

# drop species with < 3 sites after clumping
summary_C_post <- species_climate %>%
  group_by(Species) %>%
  summarise(n_C_sites = n(), .groups = "drop")

to_drop <- summary_C_post %>%
  filter(n_C_sites < 3) %>%
  pull(Species)

# keep only species that pass the threshold
species_climate <- species_climate %>%
  filter(!Species %in% to_drop)

writexl::write_xlsx(
  species_climate, '2025_11_04 FRM_SpeciesClimate.xlsx'
 )

rm(layout_matrix, species_grid, species_split, Batch1_species, cols, i, species, sp_grid, sp_climate, clim_vars,
   summary_C_pre, pre_lt3_total, pre_lt3_distinct, summary_C_post, to_drop)

## Grab and format Wild data - Plant ATLAS2020
#-----------------------------------------------------------------------------------------------------------
## Load and format data
files <- list.files(path = "Distribution", pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(files, function(file) {
  
  df <- read.csv(file)
  df$filename <- tools::file_path_sans_ext(basename(file))
  df
  
})

atlas.df <- do.call(rbind, data_list)
atlas.df <- subset(atlas.df, gr.projection == 'gb', select = c(filename, gr, lon, lat))
colnames(atlas.df) <- c('Species', 'GR', 'Lon', 'Lat')
table(atlas.df$Species)

## Extract climate
atlas_climate <- as.data.frame(terra::extract(sz_stack, atlas.df[ , 3:4]))
atlas_climate <- atlas_climate %>% mutate_all(as.numeric)
atlas_climate <- cbind(atlas.df$Species, atlas_climate)
colnames(atlas_climate)[1] <- "Species"
atlas_climate <- atlas_climate[complete.cases(atlas_climate), ]

table(atlas_climate$Species)

rm(files, data_list, atlas.df)

## Format final dataset
#-----------------------------------------------------------------------------------------------------------
atlas_climate$Group <- 'W'
atlas_climate <- atlas_climate %>%
  dplyr::mutate(Weight = 1) %>%
  dplyr::select(-Altitude, -Seed_zone) %>%
  dplyr::select(Species, Group, Weight, everything())

species_climate2 <- species_climate %>%
  filter(Species %in% atlas_climate$Species) %>%
  mutate(Group = "C")

# Exploring and normalizing quantities
ggplot(species_climate2, aes(x = Species, y = Quantity)) +
  geom_boxplot(outlier.color = "red", fill = "lightblue") +
  facet_wrap(~ Species, scales = "free") +
  theme_minimal() +
  ylab("Quantity (kg)") +
  theme(axis.text.x = element_blank())

species_climate2 <- species_climate2 %>%
  group_by(Species) %>%
  mutate(Weight = scales::rescale(Quantity, to = c(1e-6, 0.999999))) %>%
  dplyr::select(-Quantity) %>%
  dplyr::select(Species, Group, Weight, everything()) %>%
  ungroup()

summary(species_climate2)

# Lost and found
setdiff(
  species_climate  %>% distinct(Species) %>% pull(Species),
  species_climate2 %>% distinct(Species) %>% pull(Species)
)

atlas_climate <- atlas_climate %>%
  filter(Species %in% species_climate2$Species)

full_climate <- bind_rows(atlas_climate, species_climate2) %>%
  mutate(Group = factor(Group)) %>%
  as.data.frame()

levels(factor(full_climate$Species))
nlevels(factor(full_climate$Species))
table(full_climate$Species)
table(full_climate$Species[full_climate$Group == 'W'])
table(full_climate$Species[full_climate$Group == 'C'])
tibble::glimpse(full_climate)

writeLines(sort(unique(full_climate$Species)), "species_list.txt")
save.image("fullclimate.RData")
