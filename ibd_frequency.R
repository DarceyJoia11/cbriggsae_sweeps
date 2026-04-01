#this package is fast and efficient for reading in larger datasets:
#lets me use "fread" later instead of regular read.table, which is more efficient
library(data.table)

#Define paths:
#Folder containing your per-chromosome IBD files
data_folder <- "/mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae"
#Output prefix
output_prefix <- file.path(data_folder, "ibd_bins")

#Bin size (in bp), 50 kb
bin_size <- 50000

#List all chromosome files from germline output
chrom_files <- list.files(
  path = data_folder,
  pattern = "^australia_.*_ibd$",
  full.names = TRUE
)

#some large numbers were converted to scientific notation, remove these
options(scipen = 999)

#Loop through each chromosome file
for (file in chrom_files) {
  
  # Skip empty files
#(My script kept crashing earlier)
  if (file.info(file)$size == 0) {
    cat("Skipping empty file:", file, "\n")
    next
  }
  
  #Use filename as chromosome name (this is for naming output files later)
  chrom_name <- gsub(".*\\/|\\..*$", "", file)  # strips path and extension
  cat("Processing chromosome file:", chrom_name, "\n")
  
  #Read IBD data
  #Columns: V1=ind1, V2=ind2, V3=start, V4=end, V5=length_cM, V6=SNPs, V7=score
#fread is faster than read
  dt <- fread(file, header = FALSE)
  
  #Remove self-comparisons
  dt <- dt[V1 != V2]
  
  #Total unique strains in this chromosome
	#(Needed for later when calculating proportion of individuals with IBD)
  all_inds <- unique(c(dt$V1, dt$V2))
  N <- length(all_inds)
  
  #Assign start and end bins
#divides genomic position by 5000bp (bin size) to work out which bin we are in
  dt[, start_bin := floor(V3 / bin_size)]
  dt[, end_bin   := floor(V4 / bin_size)]
  

  # Expand segments into rows per bin for each individual
	#For each IBD segment, find the start and end (which bins it spans)
	#Then for the 2 strains sharing that segment, create a row with a column for the bin number and the individual
	#(The idea is to create a list of individuals and bins that they are found in)
  expanded_list <- lapply(1:nrow(dt), function(i){
    bins <- seq(dt$start_bin[i], dt$end_bin[i])
    inds <- c(dt$V1[i], dt$V2[i])
    data.table(bin = rep(bins, each = 2), individual = rep(inds, times = length(bins)))
  })
  #merge all of this into one big table
  expanded <- rbindlist(expanded_list)
  
  # Keep unique individuals per bin (remove duplicates so that each individual is counted only once per bin)
  expanded_unique <- unique(expanded[, .(bin, individual)])
  
  #Then, count unique individuals per bins
	#The idea is to find out how many strains share this segment IBD with another strain
	#And calculate proportion also
  counts <- expanded_unique[, .(n_individuals = .N), by = bin]
  counts[, freq := n_individuals / N]
  
  #The output files were not sorted in any particular order:
	#to fix this (so that bins are organised in ascending order by bin number:
	#Convert bins back to genomic positions
	
  counts[, start := bin * bin_size]
  counts[, end   := (bin + 1) * bin_size - 1]
  
  # Sort by (genomic) start position
  setorder(counts, start)
  
  #Add sequential bin number column to output (1 = first 50kb, 2 = second, etc.)
  counts[, bin_number := 1:.N]
  
  #Reorder columns (this is just formatting)
  counts <- counts[, .(bin_number, start, end, n_individuals, freq)]
  
  #Save output per chromosome
  out_file <- paste0(output_prefix, "_", chrom_name, ".txt")
  fwrite(counts, out_file)
}

cat("All chromosomes processed!\n")
