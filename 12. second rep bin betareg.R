# PREREQUISITS: load section 1,9,10 data 

hg19.gaps <- read.csv("~/Documents/CNA/Github/singleBiopsyITH/Data/hg19gap.csv")
hg19.gaps <- hg19.gaps[hg19.gaps$type %in% c('telomere',"centromere"),]
hg19.gaps <- hg19.gaps[hg19.gaps$chrom %!in% c('chrY','chrX'),]

# CHOOSE ANY BIN WITH MAX COEFF ----
representitive.bins_secondrep <- list()
for ( j in 1:500 ) {
  # run j iterations
  print(j)
  list <- list()
  c <- unique(candidate.bins$cluster)
  # for each cluster chose a rep bin that max coeff
  for ( i in 1:length(c) ) {
    cluster <- c[i]
    data <- candidate.bins[which(candidate.bins$cluster == cluster),]
    data <- data[which(abs(data$coeff) == max(abs(data$coeff))),]
    set.seed(j)
    list[[i]] <- sample_n(data,1)
  }
  representitive.bins_secondrep[[j]] <- do.call('rbind',list)
}

# remove bins with overlap (>50%) with telo or centro
for (k in 1:length(representitive.bins_secondrep)) {
  representitive.bins_secondrep[[k]]$teloCentro <- NA
  
  for (i in 1:nrow(representitive.bins_secondrep[[k]])) {
    bin <- representitive.bins_secondrep[[k]]$bin[i]
    chr <- paste('chr',car.info$start.stop$chr[which(car.info$start.stop$bin==bin)], sep = "")
    start <- car.info$start.stop$start[which(car.info$start.stop$bin==bin)]
    stop <- car.info$start.stop$stop[which(car.info$start.stop$bin==bin)]
    
    wd <- hg19.gaps[which(hg19.gaps$chrom==chr),]
    
    # check overlap
    overlap <- list()
    for (j in 1:nrow(wd)) {
      overlap[[j]] <- Overlap(c(start,stop), c(wd$chromStart[j], wd$chromEnd[j]))
    }
    representitive.bins_secondrep[[k]]$teloCentro[i] <- max(unlist(overlap))
    
    # as a percent of bin size
    representitive.bins_secondrep[[k]]$teloCentro[i] <- representitive.bins_secondrep[[k]]$teloCentro[i]/(stop-start)
  }
  
  representitive.bins_secondrep[[k]] <- representitive.bins_secondrep[[k]][which(representitive.bins_secondrep[[k]]$teloCentro < 0.5),]
}

# perform variable selection
bins_secondrep <- list()
clusters_secondrep <- list()
mods_secondrep <- list()
for ( j in 1:10 ) {
  print(j)
  reps <- representitive.bins_secondrep[[j]]
  
  multiReg.in <- t(car.raw[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
  multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
  multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
  rownames(multiReg.in) <- car.info$sampleIDs
  multiReg.in <- multiReg.in[,c(representitive.bins$bin)] # Keep only the candidate bins (columns)
  
  names(multiReg.in) = paste("bin_", names(multiReg.in), sep="")
  multiReg.in = multiReg.in %>%
    mutate(across(everything(), as.character))
  
  # Make dummy for beta
  multiReg.in <- dummy_cols(multiReg.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  multiReg.in$bin_ITH <- car.actualITH$actual
  remove <- attributes(alias(lm(bin_ITH ~ ., data = multiReg.in))$Complete)$dimnames[[1]]
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
  
  # Backwards AIC selection on 51 bin model
  x <- multiReg.in[ ,-which(colnames(multiReg.in) %in% c("bin_ITH"))]
  y <- car.actualITH$actual
  backsel <- betaselect(x, y, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
  keep <- backsel$variable
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %in% keep]
  multiReg.in$bin_ITH <- car.actualITH$actual
  mod <- betareg(bin_ITH ~., data = multiReg.in)
  x <- data.frame(vif(mod))
  
  while( max(x$vif.mod.) >=10 ) {
    remove <- rownames(x)[which(x$vif.mod.== max(x$vif.mod.))]
    multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
    mod <- betareg(bin_ITH ~., data = multiReg.in)
    x <- data.frame(vif(mod))
  }

  mods_secondrep[[j]] <- mod
  bins <- rownames(data.frame(summary(mod)$coefficients))[-1]
  bins <- data.frame(bin=as.numeric(str_extract_all(bins, "[0-9]+")), cna=str_sub(bins,-4,-1))
  bins$cluster <- candidate.bins[match(bins$bin, candidate.bins$bin), 11] 
  
  bins_secondrep[[j]] <- bins
  clusters_secondrep[[j]] <- unique(bins$cluster)
}

# are clusters the same across model options
for (i in 1:length(clusters_secondrep)) {
  if ( sum(clusters_secondrep[[i]] %!in% unique(hg19predictors$cluster)) >0 | sum(unique(hg19predictors$cluster) %!in% clusters_secondrep[[i]]) >0 ) {
    print(paste(i,"cluster mismatch"))
  }
}
for (i in 1:length(mods_secondrep)) {
  print(round(AIC(mods_secondrep[[i]]),2))
}
y <- list()
for (i in 1:length(bins_secondrep)) {
  y[[i]] <- bins_secondrep[[i]]$bin
}
y <- do.call('cbind', y)


# CHOOSE ANY BIN WITH NON-MAX COEFF ----
representitive.bins_nonmax <- list()
for ( j in 1:10 ) {
  # Select one representive bin per cluster
  print(j)
  list <- list()
  c <- unique(candidate.bins$cluster)
  for ( i in 1:length(c) ) {
    cluster <- c[i]
    data <- candidate.bins[which(candidate.bins$cluster == cluster),]
    if (nrow(data)==1) {
      list[[i]] <- data
    }
    else if (length(unique(data$coeff))==1) {
      list[[i]] <- sample_n(data,1)
    }
    else {
      data <- data[which(abs(data$coeff) != max(abs(data$coeff))),]
      set.seed(j)
      list[[i]] <- sample_n(data,1)
    }
  }
  representitive.bins_nonmax[[j]] <- do.call('rbind',list)
}

# remove bins with overlap (>50%) with telo or centro
for (k in 1:length(representitive.bins_nonmax)) {
  representitive.bins_nonmax[[k]]$teloCentro <- NA
  
  for (i in 1:nrow(representitive.bins_nonmax[[k]])) {
    bin <- representitive.bins_nonmax[[k]]$bin[i]
    chr <- paste('chr',car.info$start.stop$chr[which(car.info$start.stop$bin==bin)], sep = "")
    start <- car.info$start.stop$start[which(car.info$start.stop$bin==bin)]
    stop <- car.info$start.stop$stop[which(car.info$start.stop$bin==bin)]
    
    wd <- hg19.gaps[which(hg19.gaps$chrom==chr),]
    
    # check overlap
    overlap <- list()
    for (j in 1:nrow(wd)) {
      overlap[[j]] <- Overlap(c(start,stop), c(wd$chromStart[j], wd$chromEnd[j]))
    }
    representitive.bins_nonmax[[k]]$teloCentro[i] <- max(unlist(overlap))
    
    # as a percent of bin size
    representitive.bins_nonmax[[k]]$teloCentro[i] <- representitive.bins_nonmax[[k]]$teloCentro[i]/(stop-start)
  }
  
  representitive.bins_nonmax[[k]] <- representitive.bins_nonmax[[k]][which(representitive.bins_nonmax[[k]]$teloCentro < 0.5),]
}

# perform variable selection
bins_nonmax <- list()
clusters_nonmax <- list()
mods_nonmax <- list()
for ( j in 1:10 ) {
  print(j)
  reps <- representitive.bins_nonmax[[j]]
  
  multiReg.in <- t(car.raw[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
  multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
  multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
  rownames(multiReg.in) <- car.info$sampleIDs
  multiReg.in <- multiReg.in[,c(representitive.bins$bin)] # Keep only the candidate bins (columns)
  
  names(multiReg.in) = paste("bin_", names(multiReg.in), sep="")
  multiReg.in = multiReg.in %>%
    mutate(across(everything(), as.character))
  
  # Make dummy for beta
  multiReg.in <- dummy_cols(multiReg.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  multiReg.in$bin_ITH <- car.actualITH$actual
  remove <- attributes(alias(lm(bin_ITH ~ ., data = multiReg.in))$Complete)$dimnames[[1]]
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
  
  # Backwards AIC selection on 51 bin model
  x <- multiReg.in[ ,-which(colnames(multiReg.in) %in% c("bin_ITH"))]
  y <- car.actualITH$actual
  backsel <- betaselect(x, y, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
  keep <- backsel$variable
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %in% keep]
  multiReg.in$bin_ITH <- car.actualITH$actual
  mod <- betareg(bin_ITH ~., data = multiReg.in)
  x <- data.frame(vif(mod))
  
  while( max(x$vif.mod.) >=10 ) {
    remove <- rownames(x)[which(x$vif.mod.== max(x$vif.mod.))]
    multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
    mod <- betareg(bin_ITH ~., data = multiReg.in)
    x <- data.frame(vif(mod))
  }
  
  mods_nonmax[[j]] <- mod
  bins <- rownames(data.frame(summary(mod)$coefficients))[-1]
  bins <- data.frame(bin=as.numeric(str_extract_all(bins, "[0-9]+")), cna=str_sub(bins,-4,-1))
  bins$cluster <- candidate.bins[match(bins$bin, candidate.bins$bin), 11] 
  
  bins_nonmax[[j]] <- bins
  clusters_nonmax[[j]] <- unique(bins$cluster)
}


# are clusters the same across model options
for (i in 1:length(clusters_nonmax)) {
  if ( sum(clusters_nonmax[[i]] %!in% unique(hg19predictors$cluster)) >0 | sum(unique(hg19predictors$cluster) %!in% clusters_nonmax[[i]]) >0 ) {
    print(paste(i,"cluster mismatch"))
  }
}

for (i in 1:length(mods_nonmax)) {
  #print(round(mods_nonmax[[i]]$pseudo.r.squared,3))
  print(round(AIC(mods_nonmax[[i]]),2))
}

# examine bins
y <- list()
for (i in 1:length(bins_nonmax)) {
  y[[i]] <- bins_nonmax[[i]]$bin
}
y <- do.call('cbind', y)


# CHOOSE ANY RANDOM BIN ----
representitive.bins_rand <- list()
for ( j in 1:10 ) {
  # Select one representive bin per cluster
  print(j)
  list <- list()
  c <- unique(candidate.bins$cluster)
  for ( i in 1:length(c) ) {
    cluster <- c[i]
    data <- candidate.bins[which(candidate.bins$cluster == cluster),]
    set.seed(j)
    list[[i]] <- sample_n(data,1)
    }
  representitive.bins_rand[[j]] <- do.call('rbind',list)
}

# remove bins with overlap (>50%) with telo or centro
for (k in 1:length(representitive.bins_rand)) {
  representitive.bins_rand[[k]]$teloCentro <- NA
  
  for (i in 1:nrow(representitive.bins_rand[[k]])) {
    bin <- representitive.bins_rand[[k]]$bin[i]
    chr <- paste('chr',car.info$start.stop$chr[which(car.info$start.stop$bin==bin)], sep = "")
    start <- car.info$start.stop$start[which(car.info$start.stop$bin==bin)]
    stop <- car.info$start.stop$stop[which(car.info$start.stop$bin==bin)]
    
    wd <- hg19.gaps[which(hg19.gaps$chrom==chr),]
    
    # check overlap
    overlap <- list()
    for (j in 1:nrow(wd)) {
      overlap[[j]] <- Overlap(c(start,stop), c(wd$chromStart[j], wd$chromEnd[j]))
    }
    representitive.bins_rand[[k]]$teloCentro[i] <- max(unlist(overlap))
    
    # as a percent of bin size
    representitive.bins_rand[[k]]$teloCentro[i] <- representitive.bins_rand[[k]]$teloCentro[i]/(stop-start)
  }
  
  representitive.bins_rand[[k]] <- representitive.bins_rand[[k]][which(representitive.bins_rand[[k]]$teloCentro < 0.5),]
}

# perform variable selection
bins_rand <- list()
clusters_rand <- list()
mods_rand <- list()
for ( j in 1:5 ) {
  print(j)
  reps <- representitive.bins_rand[[j]]
  
  multiReg.in <- t(car.raw[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
  multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
  multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
  rownames(multiReg.in) <- car.info$sampleIDs
  multiReg.in <- multiReg.in[,c(representitive.bins$bin)] # Keep only the candidate bins (columns)
  
  names(multiReg.in) = paste("bin_", names(multiReg.in), sep="")
  multiReg.in = multiReg.in %>%
    mutate(across(everything(), as.character))
  
  # Make dummy for beta
  multiReg.in <- dummy_cols(multiReg.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  multiReg.in$bin_ITH <- car.actualITH$actual
  remove <- attributes(alias(lm(bin_ITH ~ ., data = multiReg.in))$Complete)$dimnames[[1]]
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
  
  # Backwards AIC selection on 51 bin model
  x <- multiReg.in[ ,-which(colnames(multiReg.in) %in% c("bin_ITH"))]
  y <- car.actualITH$actual
  backsel <- betaselect(x, y, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
  keep <- backsel$variable
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %in% keep]
  multiReg.in$bin_ITH <- car.actualITH$actual
  mod <- betareg(bin_ITH ~., data = multiReg.in)
  x <- data.frame(vif(mod))
  
  while( max(x$vif.mod.) >=10 ) {
    remove <- rownames(x)[which(x$vif.mod.== max(x$vif.mod.))]
    multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
    mod <- betareg(bin_ITH ~., data = multiReg.in)
    x <- data.frame(vif(mod))
  }
  
  mods_rand[[j]] <- mod
  bins <- rownames(data.frame(summary(mod)$coefficients))[-1]
  bins <- data.frame(bin=as.numeric(str_extract_all(bins, "[0-9]+")), cna=str_sub(bins,-4,-1))
  bins$cluster <- candidate.bins[match(bins$bin, candidate.bins$bin), 11] 
  
  bins_rand[[j]] <- bins
  clusters_rand[[j]] <- unique(bins$cluster)
}

# are clusters the same across model options
for (i in 1:length(clusters_rand)) {
  if ( sum(clusters_rand[[i]] %!in% unique(hg19predictors$cluster)) >0 | sum(unique(hg19predictors$cluster) %!in% clusters_rand[[i]]) >0 ) {
    print(paste(i,"cluster mismatch"))
  }
}

for (i in 1:length(mods_rand)) {
  #print(round(mods_nonmax[[i]]$pseudo.r.squared,3))
  print(round(AIC(mods_rand[[i]]),2))
}

# examine bins
y <- list()
for (i in 1:length(bins_rand)) {
  y[[i]] <- bins_rand[[i]]$bin
}
y <- do.call('cbind', y)


# save ----
saveRDS(mods_secondrep, "~/Documents/CNA/Github/singleBiopsyITH/Data/mods_secondrep.rds")
saveRDS(bins_secondrep, "~/Documents/CNA/Github/singleBiopsyITH/Data/bins_secondrep.rds")
saveRDS(clusters_secondrep, "~/Documents/CNA/Github/singleBiopsyITH/Data/clusters_secondrep.rds")
saveRDS(representitive.bins_secondrep, "~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins_secondrep.rds")

saveRDS(mods_nonmax, "~/Documents/CNA/Github/singleBiopsyITH/Data/mods_nonmax.rds")
saveRDS(bins_nonmax, "~/Documents/CNA/Github/singleBiopsyITH/Data/bins_nonmax.rds")
saveRDS(clusters_nonmax, "~/Documents/CNA/Github/singleBiopsyITH/Data/clusters_nonmax.rds")
saveRDS(representitive.bins_nonmax, "~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins_nonmax.rds")

saveRDS(mods_rand, "~/Documents/CNA/Github/singleBiopsyITH/Data/mods_rand.rds")
saveRDS(bins_rand, "~/Documents/CNA/Github/singleBiopsyITH/Data/bins_rand.rds")
saveRDS(clusters_rand, "~/Documents/CNA/Github/singleBiopsyITH/Data/clusters_rand.rds")
saveRDS(representitive.bins_rand, "~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins_rand.rds")
