# VALIDATION COHORT 

# In order to use upcoming functions, data must be in format of a list with a dataframe per patient, with columns: 
# chr | start | stop | sample1 | sample2...

# Read in the datasets
test_car1 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.01.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car2 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car3 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car4 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.04.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car5 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car6 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.06.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car7 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.07.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car8 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.08.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9p <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.09.Proximal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9d <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.09.Distal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car10 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.10.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

# Load patient datasets into a list
validationCar.list <- list(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)
rm(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)

# Convert data to numeric
validationCar.list <- lapply(validationCar.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})
x <- validationCar.list[[1]]

# Generate absolute copy number (minor + major)
for ( j in 1:length(validationCar.list)) {
  i <- 5
  l <- 1
  abscn <- list()
  
  while ( i < ncol(validationCar.list[[j]]) ) {
    
    for ( k in 1:nrow(validationCar.list[[j]]) ) {
      abscn[[l]] <- sum(validationCar.list[[j]][k,i], validationCar.list[[j]][k,i+1])
      l <- l + 1
    }
    i <- i + 2
  }
  
  validationCar.list[[j]] <- cbind(validationCar.list[[j]][,c(1,2,4)], matrix(unlist(abscn), nrow = nrow(validationCar.list[[j]])) )
  colnames(validationCar.list[[j]])[c(1:3)] <- c('chr','start','stop')
  colnames(validationCar.list[[j]])[c(4:ncol(validationCar.list[[j]]))] <- paste(j, seq(1, ncol(validationCar.list[[j]])-3, 1), sep = ".")
  validationCar.list[[j]] <- data.frame(validationCar.list[[j]])
}
x <- validationCar.list[[2]]

# Assess ploidy and recentre on a per sample basis
for ( i in 1:length(validationCar.list) ) {
  weights <- validationCar.list[[i]][,3] - validationCar.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(validationCar.list[[i]]) ) {
    ploidy <- round(sum(validationCar.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      print(paste(i, ploidy))
      validationCar.list[[i]][ ,j] <- validationCar.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}
x <- validationCar.list[[2]]

# Bin CN data to match training dataset (both hg19)
validationCarBinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = validationCar.list)

validationCarBinned.list <- lapply(validationCarBinned.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(validationCarBinned.list) ) {
  for ( k in 5:ncol(validationCarBinned.list[[i]]) ) {
    colnames(validationCarBinned.list[[i]])[k] <- paste(i, k-4, sep = '.')
  }
}

# make dataframe
validationCarBinned.list <- lapply(validationCarBinned.list, function(x) {
  x <- as.data.frame(x)
  x
})

# If absolute copy number >=3, its a gain, if <=2 its a loss
validationCarBinned.list <- lapply(validationCarBinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x <- as.data.frame(x)
  x
})

# Convert data to numeric (removes characters)
validationCarBinned.list <- lapply(validationCarBinned.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))), check.names = F)
  x
})

# remove bin column
for (i in 1:length(validationCarBinned.list)) {
  validationCarBinned.list[[i]] <- validationCarBinned.list[[i]][-1]
}

# collapse into one dataframe
validationCar.raw <- validationCarBinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop"))

# pull data info, clonality, diversirty
validationCar.info <- PullDataInfo(rawdata = validationCar.raw) # hg19
validationCar.clonality <- PullDataClonality(rawdata = validationCar.raw, dataInfo = validationCar.info)
validationCar.diversity <- PullDataDiversity(rawdata = validationCar.raw, dataInfo = validationCar.info)
validationCar.actualITH <- data.frame(patient = lapply(data.frame(patient=validationCar.info$patientIDs), rep, validationCar.info$sampPerPatient),
                                      sample = validationCar.info$sampleIDs,
                                      actual = lapply(data.frame(actual=validationCar.diversity$pic.frac$pic.frac), rep, validationCar.info$sampPerPatient))

# Save
saveRDS(validationCar.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.raw.rds")
saveRDS(validationCar.info, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.info.rds")
saveRDS(validationCar.clonality, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.clonality.rds")
saveRDS(validationCar.diversity, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.diversity.rds")
saveRDS(validationCar.actualITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.actualITH.rds")

