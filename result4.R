library(data.table)
library(ggplot2)
library(scales)

# PATHS
input_dir <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
output_dir <- "/mnt/loki/hartfield/caenorhabditis/results/Darcey/Cbriggsae"

chroms <- c("I","II","III","IV","V","X")

genome_size <- 106181468

# LOAD GERMLINE DATA
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

# PAIRWISE SUMMARY
pair_stats <- all_ibd[, .(
  total_ibd_bp = sum(length_bp),
  n_segments = .N
), by = pair]

pair_stats[, prop_percent := (total_ibd_bp / genome_size) * 100]

# BINNING

x_breaks <- seq(0, 100, 10)
pair_stats[, x_bin := cut(prop_percent,
                          breaks = x_breaks,
                          include.lowest = TRUE,
                          right = FALSE)]

max_y <- ceiling(max(pair_stats$n_segments) / 10) * 10
y_breaks <- seq(0, max_y, 10)

pair_stats[, y_bin := cut(n_segments,
                          breaks = y_breaks,
                          include.lowest = TRUE,
                          right = FALSE)]

# GRID COUNTS
grid_counts <- pair_stats[, .N, by = .(x_bin, y_bin)]

# convert bins → numeric midpoints 
grid_counts[, x := x_breaks[as.numeric(x_bin)] + 5]
grid_counts[, y := y_breaks[as.numeric(y_bin)] + 5]

grid_counts <- na.omit(grid_counts)

# DISCRETE COLOUR BINS 
grid_counts[, N_class := cut(
  N,
  breaks = c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, Inf),
  right = FALSE,
  include.lowest = TRUE
)]

legend_levels <- levels(grid_counts$N_class)

legend_labels <- c("<5", "5", "10", "20", "30", "40", "50", "60", "70", "80")

grid_counts[, N_class := factor(N_class,
                                levels = legend_levels,
                                labels = legend_labels)]

# COLOUR SCALE
fill_cols <- c(
  "<5"  = "#dbe9f6",
  "5"   = "#cfe1f2",
  "10"  = "#b8d4ea",
  "20"  = "#a1c6e2",
  "30"  = "#8fbce0",
  "40"  = "#74aedb",
  "50"  = "#5fa6d6",
  "60"  = "#4a98cf",
  "70"  = "#2f89c5",
  "80" = "#0b4f7a"
)

# PLOT

p <- ggplot(grid_counts, aes(x = x, y = y, fill = N_class)) +

  geom_tile(width = 10, height = 10, color = "white", linewidth = 0.4) +

  scale_fill_manual(values = fill_cols, name = "Strain pairs") +

  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    limits = c(0, max_y),
    breaks = seq(0, max_y, 10),
    expand = c(0, 0)
  ) +

  labs(
    x = "Proportion of genome shared IBD",
    y = "Number of IBD tracts per pair"
  ) +

  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(face = "bold")
  )

# SAVE

ggsave(
  file.path(output_dir, "heatmap.pdf"),
  p,
  width = 7.5,
  height = 6
)
