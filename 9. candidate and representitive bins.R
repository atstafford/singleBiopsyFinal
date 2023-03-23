# PREREQUISITS: load section 1,7,8 data 

# Candidate bins

# Bind the three dataframes holding the univariate regression data
# pull bins that were individually significant (p<0.05) predictors of ITH
candidate.bins <- do.call('rbind',uniReg.out.list[c(2:3)])
candidate.bins <- candidate.bins[which(candidate.bins$pval <= 0.05),]

# add clonality
candidate.bins$pcSubclonal <- NA
for ( i in 1:nrow(candidate.bins) ) {
  bin <- candidate.bins$bin[i]
  CNA <- candidate.bins$CNA[i]
  
  if ( CNA=='gain') {
    candidate.bins$pcSubclonal[i] <- car.clonality$pcSubclonal$gain[which(car.clonality$pcSubclonal$bin==bin)]
  }
  if ( CNA=='loss') {
    candidate.bins$pcSubclonal[i] <- car.clonality$pcSubclonal$loss[which(car.clonality$pcSubclonal$bin==bin)]
  }
}

candidate.bins$cluster <- sig.hclust[match(candidate.bins$bin, sig.hclust$bin), 2]
rownames(candidate.bins) <- 1:nrow(candidate.bins)

# Assign unique clusters to bins with no cluster
new.clust <- max(candidate.bins$cluster[is.finite(candidate.bins$cluster)]) + 1
for (i in 1:nrow(candidate.bins) ) {
  if ( is.na(candidate.bins$cluster[i]) ) {
    candidate.bins$cluster[i] <- new.clust
    new.clust <- new.clust + 1
  }
}

# Pull representative bin: middle bin with the max coeff
list <- list()
c <- unique(candidate.bins$cluster)
for ( i in 1:length(c) ) {
  cluster <- c[i]
  data <- candidate.bins[which(candidate.bins$cluster == cluster),]
  data <- data[which(abs(data$coeff) == max(abs(data$coeff))),]
  if (nrow(data)==1) {
    list[[i]] <- data
  }
  else {
    middle <- round(nrow(data)/2)
    list[[i]] <- data[middle,]
  }
}
representitive.bins <- do.call('rbind',list)

# remove bins with overlap (>50%) with telo or centro
hg19.gaps <- read.csv("~/Documents/CNA/Github/singleBiopsyITH/Data/hg19gap.csv")
hg19.gaps <- hg19.gaps[hg19.gaps$type %in% c('telomere',"centromere"),]
hg19.gaps <- hg19.gaps[hg19.gaps$chrom %!in% c('chrY','chrX'),]
representitive.bins$teloCentro <- NA

for (i in 1:nrow(representitive.bins)) {
  bin <- representitive.bins$bin[i]
  chr <- paste('chr',car.info$start.stop$chr[which(car.info$start.stop$bin==bin)], sep = "")
  start <- car.info$start.stop$start[which(car.info$start.stop$bin==bin)]
  stop <- car.info$start.stop$stop[which(car.info$start.stop$bin==bin)]
  
  wd <- hg19.gaps[which(hg19.gaps$chrom==chr),]
  
  # check overlap
  overlap <- list()
  for (j in 1:nrow(wd)) {
    overlap[[j]] <- Overlap(c(start,stop), c(wd$chromStart[j], wd$chromEnd[j]))
  }
  representitive.bins$teloCentro[i] <- max(unlist(overlap))
  
  # as a percent of bin size
  representitive.bins$teloCentro[i] <- representitive.bins$teloCentro[i]/(stop-start)
}

representitive.bins <- representitive.bins[which(representitive.bins$teloCentro < 0.5),]

# save
saveRDS(candidate.bins, "~/Documents/CNA/Github/singleBiopsyITH/Data/candidate.bins.rds")
saveRDS(representitive.bins, "~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins.rds")