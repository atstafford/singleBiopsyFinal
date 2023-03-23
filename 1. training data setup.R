# Load data
rawdata <- read.table("~/Documents/CNA/Data/Training/sWGS Cross rename with type.txt", header=T, skip=2, sep="\t")

# Split into 2 separate raw data files (ad.raw and car.raw) 
ad.raw <- rawdata[,c(1:3,which(sub("^([[:alpha:]]*).*", "\\1", colnames(rawdata))=="A"))]
car.raw <- rawdata[,c(1:3,which(sub("^([[:alpha:]]*).*", "\\1", colnames(rawdata))=="C"))]

# pull data info, clonality, diversity
ad.info <- PullDataInfo(rawdata = ad.raw) # hg19
ad.clonality <- PullDataClonality(rawdata = ad.raw, dataInfo = ad.info)
ad.diversity <- PullDataDiversity(rawdata = ad.raw, dataInfo = ad.info)
ad.actualITH <- data.frame(patient = lapply(data.frame(patient=ad.info$patientIDs), rep, ad.info$sampPerPatient),
                           sample = ad.info$sampleIDs,
                           actual = lapply(data.frame(actual=ad.diversity$pic.frac$pic.frac), rep, ad.info$sampPerPatient))

car.info <- PullDataInfo(rawdata = car.raw) # hg19
car.clonality <- PullDataClonality(rawdata = car.raw, dataInfo = car.info)
car.diversity <- PullDataDiversity(rawdata = car.raw, dataInfo = car.info)
car.actualITH <- data.frame(patient = lapply(data.frame(patient=car.info$patientIDs), rep, car.info$sampPerPatient),
                            sample = car.info$sampleIDs,
                            actual = lapply(data.frame(actual=car.diversity$pic.frac$pic.frac), rep, car.info$sampPerPatient))

# set up data matrix
ad.matrices <- genMatrices(ad.raw)
car.matrices <- genMatrices(car.raw)

# save
saveRDS(ad.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.raw.rds")
saveRDS(ad.info, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.info.rds")
saveRDS(ad.clonality, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.clonality.rds")
saveRDS(ad.diversity, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.diversity.rds")
saveRDS(ad.actualITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.actualITH.rds")

saveRDS(car.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.raw.rds")
saveRDS(car.info, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.info.rds")
saveRDS(car.clonality, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.clonality.rds")
saveRDS(car.diversity, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.diversity.rds")
saveRDS(car.actualITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.actualITH.rds")

saveRDS(ad.matrices, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.matrices.rds")
saveRDS(car.matrices, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.matrices.rds")

