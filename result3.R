library(data.table)
library(ggplot2)
library(scales)

#Paths
input_dir <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
output_dir <- "/mnt/loki/hartfield/caenorhabditis/results/Darcey/Cbriggsae"

#Genome size
genome_size <- 106181468

#Chromosomes
chroms <- c("I","II","III","IV","V","X")

#Combine GERMLINE files
all_ibd <- rbindlist(lapply(chroms, function(chr) {

  file <- file.path(input_dir, paste0("australia_", chr, "_ibd"))
  if (!file.exists(file) || file.info(file)$size == 0) return(NULL)

  dt <- fread(file, header = FALSE)
  dt <- dt[V1 != V2]

  dt[, pair := ifelse(V1 < V2,
                      paste(V1, V2, sep = "_"),
                      paste(V2, V1, sep = "_"))]

  dt[, length_bp := V4 - V3]
  dt <- dt[length_bp > 0]

  dt
}))

#Summarise per pair
pair_stats <- all_ibd[, .(
  total_ibd_bp = sum(length_bp)
), by = pair]

pair_stats[, prop_percent := (total_ibd_bp / genome_size) * 100]

#Plot
p <- ggplot(pair_stats, aes(x = prop_percent)) +

  geom_histogram(
    bins = 50,              # keep your thin bars
    boundary = 0,           # ✔ align bins to start at 0
    fill = "steelblue",     # ✔ restore original colour
    color = NA              # ✔ no borders
  ) +

  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +

  labs(
    x = "Proportion of genome shared IBD",
    y = "Number of strain pairs"
  ) +

  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),                     # remove grid
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # box
  )

#Save
ggsave(
  file.path(output_dir, "proportion_genome_ibd.pdf"),
  p,
  width = 7.5,
  height = 6
)
