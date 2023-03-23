# TCGA PROCESSING ####

# Load clinical data for CRC (COAD)
COAD_clinical <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/TCGA/Clinical COAD.xlsx")
colnames(COAD_clinical)[c(1,2)] <- c("UUID_copyN","TCGA_barcode")

# Remove duplicates
COAD_clinical <- COAD_clinical[!duplicated(COAD_clinical$TCGA_barcode), ]

# Remove columns with no information
COAD_clinical <- COAD_clinical[,colSums(COAD_clinical=="'--") < nrow(COAD_clinical)]
COAD_clinical <- COAD_clinical[,-c(3,14,15,19,18,20,23,28,32,33)]

# Load additional data and merge
COAD_MSI <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/TCGA/TCGA_MSIstatus.xlsx", col_names = TRUE, skip =1)
COAD_MSI <- COAD_MSI[COAD_MSI$Organ=="COAD",]
COAD_clinical$MSI <- unlist(COAD_MSI[match(COAD_clinical$TCGA_barcode, COAD_MSI$`TCGA Participant Barcode`), 14])
COAD_clinical$CMS <- unlist(COAD_MSI[match(COAD_clinical$TCGA_barcode, COAD_MSI$`TCGA Participant Barcode`), 27])
COAD_clinical$MSI[is.na(COAD_clinical$MSI)] <- 'NA'

# Change stage names to remove A/B
colnames(COAD_clinical)[13] <- 'stage'
for ( i in 1:nrow(COAD_clinical) ) {
  COAD_clinical$stage[i] <- ifelse(COAD_clinical$stage[i]=="Stage IA", "Stage I", 
                                   ifelse(COAD_clinical$stage[i]=="Stage IIA", "Stage II",
                                          ifelse(COAD_clinical$stage[i]=="Stage IIB", "Stage II",
                                                 ifelse(COAD_clinical$stage[i]=="Stage IIC", "Stage II",
                                                        ifelse(COAD_clinical$stage[i]=="Stage IIIA", "Stage III",
                                                               ifelse(COAD_clinical$stage[i]=="Stage IIIB", "Stage III",
                                                                      ifelse(COAD_clinical$stage[i]=="Stage IIIC", "Stage III",
                                                                             ifelse(COAD_clinical$stage[i]=="Stage IVA", "Stage IV",
                                                                                    ifelse(COAD_clinical$stage[i]=="Stage IVB", "Stage IV",COAD_clinical$stage[i])))))))))
}

# Load TCGA copy number data, and keep only COAD
COAD_copyN <- read.table("~/Documents/CNA/Github/singleBiopsyITH/Data/TCGA/TCGAcn/TCGA_mastercalls.abs_segtabs.fixed.txt", header=T, skip=0, sep="\t")
COAD_copyN <- COAD_copyN[which((str_sub(COAD_copyN$Sample,1,str_length(COAD_copyN$Sample)-3)) %in% COAD_clinical$TCGA_barcode),c(1:4,9)]

# Remove patients missing CN data
COAD_copyN <- COAD_copyN[!is.na(COAD_copyN$Modal_Total_CN),]

# Rename columns
colnames(COAD_copyN) <- c('sample','chr','start','stop','cn')

# Pull cn data into element corresponding to sample ID
names <- unique(COAD_copyN$sample)
COAD.list <- list()
for ( i in 1:length(names) ) {
  COAD.list[[i]] <- COAD_copyN[which(COAD_copyN$sample == names[i]),]
  COAD.list[[i]] <- COAD.list[[i]][-1] #remove col holding sample as its sample per list element
}

# Assess ploidy and recentre 
for ( i in 1:length(COAD.list) ) {
  weights <- COAD.list[[i]][,3] - COAD.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(COAD.list[[i]]) ) {
    ploidy <- round(sum(COAD.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      print(paste(i, ploidy))
      COAD.list[[i]][ ,j] <- COAD.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}

# Bin CN data to match training dataset
COADbinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = COAD.list)
COADbinned.list <- lapply(COADbinned.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(COADbinned.list) ) {
  for ( k in 5:ncol(COADbinned.list[[i]]) ) {
    colnames(COADbinned.list[[i]])[k] <- names[i]
  }
}

# make dataframe and remove bin column
COADbinned.list <- lapply(COADbinned.list, function(x) {
  x <- as.data.frame(x)
  x <- x[-1]
  x
})

# If absolute copy number >=3, its a gain, if <=2 its a loss
COADbinned.list <- lapply(COADbinned.list, function(x) {
  x[,-c(1:3)] <- ifelse(x[,-c(1:3)] >= 3, 3, ifelse(x[,-c(1:3)] < 2, 1, 2))
  x <- as.data.frame(x)
  x
})

# Convert data to numeric (removes characters)
COADbinned.list <- lapply(COADbinned.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))), check.names = F)
  x
})

# collapse into one dataframe
coad.raw <- COADbinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop"))
colnames(coad.raw) <- paste(colnames(coad.raw), ".", sep = "")

# pull data info, clonality, diversity
coad.info <- PullDataInfo(rawdata = coad.raw) # hg19

# Save
saveRDS(coad.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/coad.raw.rds")
saveRDS(COAD_clinical, "~/Documents/CNA/Github/singleBiopsyITH/Data/COAD_clinical.rds")
saveRDS(coad.info, "~/Documents/CNA/Github/singleBiopsyITH/Data/coad.info.rds")

