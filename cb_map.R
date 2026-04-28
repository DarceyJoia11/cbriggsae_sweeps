#This script was kindly supplied by my supervisor, Matthew Hartfield

# Load in tidyverse library
library(tidyverse)

# Defining selfing rate and inbreeding coefficient (for rescaled rec rate)
# This equation comes from pop gen theory, e.g., Nordborg 2000
# Here, F=/~ 0.98, meaning individuals are almost completely homozygous.
self <- 0.99
F <- self/(2-self)

# Read in genetic map from CaeNDR. Contains chrom number, markep pos, genetic distance
# And extract data from one chromosome ("I" by default)
dat <- read_delim("c_briggsae_genetic_map.tsv")
dat2 <- dat %>% filter(chrom=="I") %>% select(marker_pos,genetic)

# Define new array for centimorgan measurements
nr <- dim(dat2)[1]  #nr = number of rows (so number of SNPs)
cm <- rep(0,nr)     #cm will store rec rate (cM/Mb)
cmF <- rep(0,nr)    #rep(0, nr) creates a vector of 0s of length nr

# For each row (starting at 2) calculate difference and hence cM per Mb
#For each adjacent pair of markers, the recombination rate is...
#The change in cM/change in bp ...
# x 10^6 (this converts from cM per bp to cM per Mb)
#Start at i = 2 as you need a previous point to compute a difference
#cm[i] = recombination rate between marker i -1 and i
for(i in 2:nr){
	cm[i] <- (as.double(dat2[i,2]-dat2[i-1,2]))/(as.double(dat2[i,1]-dat2[i-1,1]))*1e6
}

# Secondary calculation, rescaling by 1-F (see Nordborg 2000)
# High selfing = fewer recombination events
# So recombination is (effectively) reduced
# F = 0.98 and 1-F = 0.02 ...
# Recombination is ~50x lower effectively
cmF <- cm*(1-F)

# Creating new map files
# Combines columns for: position (bp), rec rate (cM/Mb) and cumulative genetic position (cM)
# Then, replaces with selfing adjusted cM rate
dat_out <- cbind(dat2$marker_pos,cm,dat2$genetic)
dat_outF <- cbind(dat2$marker_pos,cmF,dat2$genetic*(1-F))

# Saving to disk
#N.B. - I did not use self maps in my project, only regular maps
write.table(dat_out,file="chr_I.map",sep="\t",quote=FALSE,row.names= FALSE,col.names= FALSE)
write.table(dat_outF,file="chr_I_self.map",quote= FALSE,row.names= FALSE,col.names= FALSE)
