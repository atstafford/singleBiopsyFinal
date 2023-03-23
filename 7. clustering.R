# PREREQUISITS: load section 1,6 data 

# hierarchical clustering based on diploid/aneuploid correlation matrix to cluster cooccuring CNAs
input <- cna.corr.car
input <- data.frame(t(input))
hiclust_diploid.aneu <- pvclust::pvclust(data = input, method.dist = 'euclidean', method.hclust = 'complete', nboot=1000, parallel=TRUE)

# select significant clusters
load("~/Documents/CNA/HPC/hiclust_diploid.aneu.rda")
hiclust <- hiclust_diploid.aneu
sig.hclust <- pvclust::pvpick(hiclust, alpha=0.95, pv="au", max.only = T)

# For each cluster extract the corresponding bins
list <- list()
for ( i in 1:length(sig.hclust$clusters) ) { 
  list[[i]] <- data.frame(bin = unlist(sig.hclust$clusters[[i]]), cluster = i)
}

# Create a dataframe with columns: bin | cluster
sig.hclust <- do.call('rbind', list)
sig.hclust$bin <- sub('.', '', sig.hclust$bin)

# Convert dataframe to numeric
sig.hclust <- data.frame(apply(sig.hclust, 2, function(x) as.numeric(as.character(x))))

# Add chromosome data based on bin number
sig.hclust$chr <- car.info$start.stop[match(sig.hclust$bin, car.info$start.stop$bin), 2]

# Save
saveRDS(sig.hclust, "~/Documents/CNA/Github/singleBiopsyITH/Data/sig.hclust.rds")

