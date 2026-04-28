library(data.table)

input_dir <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
chroms <- c("I","II","III","IV","V","X")
genome_size <- 106181468

# Load and combine IBD data
all_ibd <- rbindlist(lapply(chroms, function(chr) {

  file <- file.path(input_dir, paste0("australia_", chr, "_ibd"))
  if (!file.exists(file) || file.info(file)$size == 0) return(NULL)

  dt <- fread(file, header = FALSE)

  # remove self comparisons
  dt <- dt[V1 != V2]

  # add chromosome
  dt[, chr := chr]

  # tract length
  dt[, length_bp := V4 - V3]
  dt <- dt[length_bp > 0]

  # consistent pair ID
  dt[, pair := ifelse(V1 < V2,
                      paste(V1, V2, sep = "_"),
                      paste(V2, V1, sep = "_"))]

  dt
}))

# Per-chromosome summary
ibd_chr_summary <- all_ibd[, .(
  mean_length   = mean(length_bp),
  median_length = median(length_bp),
  max_length    = max(length_bp)
), by = chr]

cat("\n=== IBD TRACT LENGTHS PER CHROMOSOME ===\n")
print(ibd_chr_summary)

# Pairwise IBD summary
pair_stats <- all_ibd[, .(
  total_ibd_bp = sum(length_bp),
  n_segments   = .N
), by = pair]

pair_stats[, prop_shared := total_ibd_bp / genome_size]

# Summary across pairs
pair_summary <- pair_stats[, .(
  mean_prop   = mean(prop_shared),
  median_prop = median(prop_shared),
  min_prop    = min(prop_shared),
  max_prop    = max(prop_shared)
)]

cat("\n=== PAIRWISE GENOME SHARING SUMMARY ===\n")
print(pair_summary)

# Number of pairs
cat("\nNumber of unique strain pairs:\n")
print(uniqueN(pair_stats$pair))

# Segment counts
segment_summary <- pair_stats[, .(
  mean_segments   = mean(n_segments),
  median_segments = median(n_segments)
)]

cat("\n=== IBD SEGMENTS PER PAIR ===\n")
print(segment_summary)

# Individual level IBD coverage
# For each strain, sum all IBD it participates in
indiv_ibd <- all_ibd[, .(
  ibd_bp = sum(length_bp)
), by = .(individual = V1)]

# Each strain appears twice (as V1 and V2), so include both
indiv_ibd2 <- rbind(
  all_ibd[, .(individual = V1, length_bp)],
  all_ibd[, .(individual = V2, length_bp)]
)

indiv_summary <- indiv_ibd2[, .(
  total_ibd_bp = sum(length_bp)
), by = individual]

indiv_summary[, prop_genome_covered := total_ibd_bp / genome_size]

indiv_final_summary <- indiv_summary[, .(
  mean_individual_coverage   = mean(prop_genome_covered),
  median_individual_coverage = median(prop_genome_covered),
  min_individual_coverage    = min(prop_genome_covered),
  max_individual_coverage    = max(prop_genome_covered)
)]

cat("\n=== INDIVIDUAL GENOME COVERAGE BY IBD ===\n")
print(indiv_final_summary)

cat("\nDONE\n")
