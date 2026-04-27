#!/usr/bin/env Rscript 
library(data.table)
library(ggplot2)
library(RColorBrewer)

# Paths
input_dir <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
output_dir <- "/mnt/loki/hartfield/caenorhabditis/results/Darcey/Cbriggsae"

# Chromosomes
chroms <- c("I","II","III","IV","V","X")
chrom_colors <- brewer.pal(6, "Set2")
names(chrom_colors) <- chroms

# Combine all chromosomes
all_bins <- rbindlist(lapply(chroms, function(chr) {

  bins <- fread(file.path(input_dir, paste0("ibd_bins_australia_", chr, "_ibd.txt")))

# Ensure numeric
# Freq is the fraction of strains that share AT LEAST one segment IBD 
# with ANY other strain in that genomic bin 
bins[, freq := as.numeric(freq)]
bins[, start := as.numeric(start)]
bins[, end   := as.numeric(end)]

bins[, chr := chr]
bins
}))

# Order chromosomes
all_bins[, chr := factor(chr, levels = chroms)]

# Plot

p <- ggplot(all_bins) +
  geom_segment(aes(x = start, xend = end, y = freq, yend = freq, color = chr),
               size = 2, alpha = 0.8) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  facet_wrap(~chr, scales = "free_x", ncol = 2) +
  labs(x = "Genomic position (Mb)",
       y = "Fraction of strains sharing IBD",
       color = "Chromosome") +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )
# Save
ggsave(file.path(output_dir, "ibd_bp_dashes.pdf"), p, width = 10, height = 8)
