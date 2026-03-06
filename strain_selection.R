# C. briggsae strain selection and labelling by country
# Using 20250626_c_briggsae_strain_data.csv, available on CaeNDR
# - Clean coordinates
# - Assign countries
# - Select one strain per isotype
# - Generate files for PLINK and further analyses
# --------------------

# Load necessary libraries
library(sf)             # for spatial data
library(rnaturalearth)  # for country polygons
library(dplyr)          # data manipulation
library(tidyverse)      # includes readr, ggplot2, etc.

# 1. Load data
Cbriggsae <- read_csv(file.choose())

# 2. Clean data
# Keep only rows with longitude and latitude
Cbriggsae_clean <- Cbriggsae %>%
  filter(!is.na(longitude) & !is.na(latitude))

# 3. Convert to spatial object
df_sf <- st_as_sf(
  Cbriggsae_clean,
  coords = c("longitude", "latitude"),
  crs = 4326  # WGS84
)

# 4. Load countries and assign to strains
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  rename(country = name)  # rename for clarity

Cbriggsae_result <- st_join(df_sf, world["country"])

# 5. Save the full result (all strains with country info)
write_csv(
  st_drop_geometry(Cbriggsae_result),
  "Cbriggsae_result.csv"
)

# 6. Select one strain per isotype
unique_isotypes_sf <- Cbriggsae_result %>%
  filter(!is.na(isotype)) %>%         # remove missing isotypes
  group_by(isotype) %>%               # group by isotype
  slice(1) %>%                        # select first strain per isotype
  ungroup()

# Note: slice(1) selects the first row in each group based on the original row order
# This ensures reproducibility if the input CSV order is consistent; otherwise, check strain list provided

# Generate PLINK keep file 
keep_file <- unique_isotypes_sf %>%
  st_drop_geometry() %>%             # drop geometry for table output
  select(isotype) %>%
  mutate(FID = isotype,
         IID = isotype) %>%
  select(FID, IID)

write.table(
  keep_file,
  "unique_isotypes_keep.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)

# 8. Save one-row-per-isotype table with country info
write_csv(
  st_drop_geometry(unique_isotypes_sf),
  "unique_isotypes.csv"
)

