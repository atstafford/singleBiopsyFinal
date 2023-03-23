# TCGA PROCESSING ####

# Load clinical data for CRC (COAD)
READ_clinical <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/TCGA/Clinical COAD.xlsx")
colnames(READ_clinical)[c(1,2)] <- c("UUID_copyN","TCGA_barcode")

# Remove duplicates
READ_clinical <- READ_clinical[!duplicated(READ_clinical$TCGA_barcode), ]

# Remove columns with no information
READ_clinical <- READ_clinical[,colSums(READ_clinical=="'--") < nrow(READ_clinical)]
READ_clinical <- READ_clinical[,-c(3,14,15,19,18,20,23,28,32,33)]

# Load additional data and merge
READ_MSI <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/TCGA/TCGA_MSIstatus.xlsx", col_names = TRUE, skip =1)
READ_MSI <- READ_MSI[READ_MSI$Organ=="READ",]
READ_clinical$MSI <- unlist(READ_MSI[match(READ_clinical$TCGA_barcode, READ_MSI$`TCGA Participant Barcode`), 14])
READ_clinical$CMS <- unlist(READ_MSI[match(READ_clinical$TCGA_barcode, READ_MSI$`TCGA Participant Barcode`), 27])
READ_clinical$MSI[is.na(READ_clinical$MSI)] <- 'NA'

# Change stage names to remove A/B
colnames(READ_clinical)[13] <- 'stage'
for ( i in 1:nrow(READ_clinical) ) {
  READ_clinical$stage[i] <- ifelse(READ_clinical$stage[i]=="Stage IA", "Stage I", 
                                   ifelse(READ_clinical$stage[i]=="Stage IIA", "Stage II",
                                          ifelse(READ_clinical$stage[i]=="Stage IIB", "Stage II",
                                                 ifelse(READ_clinical$stage[i]=="Stage IIC", "Stage II",
                                                        ifelse(READ_clinical$stage[i]=="Stage IIIA", "Stage III",
                                                               ifelse(READ_clinical$stage[i]=="Stage IIIB", "Stage III",
                                                                      ifelse(READ_clinical$stage[i]=="Stage IIIC", "Stage III",
                                                                             ifelse(READ_clinical$stage[i]=="Stage IVA", "Stage IV",
                                                                                    ifelse(READ_clinical$stage[i]=="Stage IVB", "Stage IV",READ_clinical$stage[i])))))))))
}

# Load TCGA copy number data, and keep only COAD
READ_copyN <- read.table("~/Documents/CNA/Github/singleBiopsyITH/Data/TCGA/TCGAcn/TCGA_mastercalls.abs_segtabs.fixed.txt", header=T, skip=0, sep="\t")
READ_copyN <- READ_copyN[which((str_sub(READ_copyN$Sample,1,str_length(READ_copyN$Sample)-3)) %in% READ_clinical$TCGA_barcode),c(1:4,9)]

# Remove patients missing CN data
READ_copyN <- READ_copyN[!is.na(READ_copyN$Modal_Total_CN),]

# Rename columns
colnames(READ_copyN) <- c('sample','chr','start','stop','cn')

# Pull cn data into element corresponding to sample ID
names <- unique(READ_copyN$sample)
READ.list <- list()
for ( i in 1:length(names) ) {
  READ.list[[i]] <- READ_copyN[which(READ_copyN$sample == names[i]),]
  READ.list[[i]] <- READ.list[[i]][-1] #remove col holding sample as its sample per list element
}

# Assess ploidy and recentre 
for ( i in 1:length(READ.list) ) {
  weights <- READ.list[[i]][,3] - READ.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(READ.list[[i]]) ) {
    ploidy <- round(sum(READ.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      print(paste(i, ploidy))
      READ.list[[i]][ ,j] <- READ.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}

# Bin CN data to match training dataset
READbinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = READ.list)
READbinned.list <- lapply(READbinned.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(READbinned.list) ) {
  for ( k in 5:ncol(READbinned.list[[i]]) ) {
    colnames(READbinned.list[[i]])[k] <- names[i]
  }
}

# make dataframe and remove bin column
READbinned.list <- lapply(READbinned.list, function(x) {
  x <- as.data.frame(x)
  x <- x[-1]
  x
})


# If absolute copy number >=3, its a gain, if <=2 its a loss
READbinned.list <- lapply(READbinned.list, function(x) {
  x[,-c(1:3)] <- ifelse(x[,-c(1:3)] >= 3, 3, ifelse(x[,-c(1:3)] < 2, 1, 2))
  x <- as.data.frame(x)
  x
})

# Convert data to numeric (removes characters)
READbinned.list <- lapply(READbinned.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))), check.names = F)
  x
})

# collapse into one dataframe
READbinned.raw <- READbinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop"))
colnames(READbinned.raw) <- paste(colnames(READbinned.raw), ".", sep = "")

# Save
saveRDS(READbinned.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/COADbinned.raw.rds")
saveRDS(READ_clinical, "~/Documents/CNA/Github/singleBiopsyITH/Data/COAD_clinical.rds")

