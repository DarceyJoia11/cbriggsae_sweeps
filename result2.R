#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)

input_dir <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
output_dir <- "/mnt/loki/hartfield/caenorhabditis/results/Darcey/Cbriggsae"

# Chromsomes
chroms <- c("I","II","III","IV","V","X")
chrom_colors <- brewer.pal(6, "Set2")
names(chrom_colors) <- chroms

# Combine GERMLINE files
all_ibd <- rbindlist(lapply(chroms, function(chr) {

file <- file.path(input_dir, paste0("australia_", chr, "_ibd"))
if (!file.exists(file) || file.info(file)$size == 0) return(NULL)

dt <- fread(file, header = FALSE)

dt <- dt[V1 != V2]

dt[, length_bp := V4 - V3]

dt <- dt[length_bp > 0]

dt[, length_Mb := length_bp / 1e6]
dt[, chr := chr]

dt
}))

all_ibd[, chr := factor(chr, levels = chroms)]
#Plot
p2 <- ggplot(all_ibd, aes(x = length_Mb, fill = chr)) +
  geom_histogram(bins = 50, color = NA) +
  scale_fill_manual(values = chrom_colors) +
  facet_wrap(~chr, ncol = 2) +   # ← FIXED HERE
  scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "IBD tract length (Mb, log scale)",
       y = "Number of IBD tracts",
       fill = "Chromosome") +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.spacing = unit(1.2, "lines")
  )
#Save
ggsave(file.path(output_dir, "ibd_length_histograms_logMb.pdf"),
p2, width = 10, height = 8)
