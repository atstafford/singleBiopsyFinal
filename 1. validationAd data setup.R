# VALIDATION COHORT 

# In order to use upcoming functions, data must be in format of a list with a dataframe per patient, with columns: 
# chr | start | stop | sample1 | sample2...

# Read in the datasets
test_ad1 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.08.WGS.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad2 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad3 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad4 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.09.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad5 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

# Load patient datasets into a list
validationAd.list <- list(test_ad1,test_ad2,test_ad3,test_ad4,test_ad5)
rm(test_ad1,test_ad2,test_ad3,test_ad4,test_ad5)

# Convert data to numeric
validationAd.list <- lapply(validationAd.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})
x <- validationAd.list[[1]]

# Generate absolute copy number (minor + major)
for ( j in 1:length(validationAd.list)) {
  i <- 5
  l <- 1
  abscn <- list()
  
  while ( i < ncol(validationAd.list[[j]]) ) {
    
    for ( k in 1:nrow(validationAd.list[[j]]) ) {
      abscn[[l]] <- sum(validationAd.list[[j]][k,i], validationAd.list[[j]][k,i+1])
      l <- l + 1
    }
    i <- i + 2
  }
  
  validationAd.list[[j]] <- cbind(validationAd.list[[j]][,c(1,2,4)], matrix(unlist(abscn), nrow = nrow(validationAd.list[[j]])) )
  colnames(validationAd.list[[j]])[c(1:3)] <- c('chr','start','stop')
  colnames(validationAd.list[[j]])[c(4:ncol(validationAd.list[[j]]))] <- paste(j, seq(1, ncol(validationAd.list[[j]])-3, 1), sep = ".")
  validationAd.list[[j]] <- data.frame(validationAd.list[[j]])
}
x <- validationAd.list[[1]]

# Assess ploidy and recentre on a per sample basis
for ( i in 1:length(validationAd.list) ) {
  weights <- validationAd.list[[i]][,3] - validationAd.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(validationAd.list[[i]]) ) {
    ploidy <- round(sum(validationAd.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      print(paste(i, ploidy))
      validationAd.list[[i]][ ,j] <- validationAd.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}
x <- validationAd.list[[1]]

# Bin CN data to match training dataset (both hg19)
validationAdBinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = validationAd.list)

validationAdBinned.list <- lapply(validationAdBinned.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(validationAdBinned.list) ) {
  for ( k in 5:ncol(validationAdBinned.list[[i]]) ) {
    colnames(validationAdBinned.list[[i]])[k] <- paste(i, k-4, sep = '.')
  }
}

# make dataframe
validationAdBinned.list <- lapply(validationAdBinned.list, function(x) {
  x <- as.data.frame(x)
  x
})

# If absolute copy number >=3, its a gain, if <=2 its a loss
validationAdBinned.list <- lapply(validationAdBinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x <- as.data.frame(x)
  x
})

# Convert data to numeric (removes characters)
validationAdBinned.list <- lapply(validationAdBinned.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))), check.names = F)
  x
})

# remove bin column
for (i in 1:length(validationAdBinned.list)) {
  validationAdBinned.list[[i]] <- validationAdBinned.list[[i]][-1]
}

# collapse into one dataframe
validationAd.raw <- validationAdBinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop"))

# pull data info, clonality, diversity
validationAd.info <- PullDataInfo(rawdata = validationAd.raw) # hg19
validationAd.clonality <- PullDataClonality(rawdata = validationAd.raw, dataInfo = validationAd.info)
validationAd.diversity <- PullDataDiversity(rawdata = validationAd.raw, dataInfo = validationAd.info)
validationAd.actualITH <- data.frame(patient = lapply(data.frame(patient=validationAd.info$patientIDs), rep, validationAd.info$sampPerPatient),
                                     sample = validationAd.info$sampleIDs,
                                     actual = lapply(data.frame(actual=validationAd.diversity$pic.frac$pic.frac), rep, validationAd.info$sampPerPatient))

# Save
saveRDS(validationAd.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.raw.rds")
saveRDS(validationAd.info, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.info.rds")
saveRDS(validationAd.clonality, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.clonality.rds")
saveRDS(validationAd.diversity, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.diversity.rds")
saveRDS(validationAd.actualITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.actualITH.rds")

