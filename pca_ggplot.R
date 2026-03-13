# PCA plot - coloured by country (no labels)

library(tidyverse)
library(RColorBrewer)

# Define directories
analysis_dir <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
results_dir  <- "/mnt/loki/hartfield/caenorhabditis/results/Darcey/Cbriggsae"

# Read PCA output (NEW FILE NAMES)
eigenvec <- read_table2(
  file.path(analysis_dir, "cbriggsae_pca2.eigenvec"),
  col_names = FALSE
)

eigenval <- scan(
  file.path(analysis_dir, "cbriggsae_pca2.eigenval")
)

# Prepare PCA dataframe
pca_df <- eigenvec %>%
  select(FID = X1, IID = X2, PC1 = X3, PC2 = X4) %>%
  mutate(
    PC1_var = round(eigenval[1] / sum(eigenval) * 100, 1),
    PC2_var = round(eigenval[2] / sum(eigenval) * 100, 1)
  )

# Read metadata
meta_file <- file.path(analysis_dir, "unique_isotypes.csv")
meta <- read_csv(meta_file)

# Merge metadata with PCA results
pca_df <- pca_df %>%
  left_join(meta, by = c("IID" = "strain"))

# Create colour column using country name
pca_df$color <- pca_df$name

# Export strains from two clusters
# (This is to identify & generate strain lists for my India and Australia populations)
# Cluster 1: PC1 > 0.10 & PC2 > 0
cluster1 <- pca_df %>%
  filter(PC1 > 0.10, PC2 > 0) %>%
  select(strain = IID, country = name)

# Cluster 2: PC1 > 0.10 & PC2 < 0
cluster2 <- pca_df %>%
  filter(PC1 > 0.10, PC2 < 0) %>%
  select(strain = IID, country = name)

# Write CSVs
write_csv(cluster1, file.path(results_dir, "cluster1_strains.csv"))
write_csv(cluster2, file.path(results_dir, "cluster2_strains.csv"))

# Generate colour palette
num_colors <- length(unique(pca_df$color))
palette_colors <- colorRampPalette(
  brewer.pal(12, "Paired")
)(num_colors)

# Create PCA plot (NO LABELS, ONLY LEGEND)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = color)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = palette_colors) +
  labs(
    x = paste0("PC1 (", unique(pca_df$PC1_var), "%)"),
    y = paste0("PC2 (", unique(pca_df$PC2_var), "%)"),
    color = "Country"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "lines")
  ) 

# Save plot
png(
  file.path(results_dir, "cbriggsae_pca2_plot.png"),
  width = 2500, height = 2000, res = 300
)
print(p)
dev.off()
