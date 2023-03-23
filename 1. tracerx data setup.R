# VALIDATION COHORT ####

# In order to use upcoming functions, data must be in format of a list with a dataframe per patient, with columns: 
# chr | start | stop | sample1 | sample2...

# Read in the datasets
tracerx <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/TracerX/cn_data.xlsx", sheet = 5)

# Pull cn data into element corresponding to sample ID
names <- unique(tracerx$sample)
patient <- unique(sub("\\-.*", "", names))

tracerx.list <- list()
for ( i in 1:length(names) ) {
  tracerx.list[[i]] <- tracerx[which(tracerx$sample == names[i]),]
  tracerx.list[[i]] <- tracerx.list[[i]][ ,c(2,3,4,6)]
  colnames(tracerx.list[[i]]) <- c("chr", "start", "stop", names[i])
  
  # make numeric
  tracerx.list[[i]] <- data.frame(apply(tracerx.list[[i]], 2, function(x) as.numeric(as.character(x))))
  
  # autosomes only
  tracerx.list[[i]] <- tracerx.list[[i]][which(tracerx.list[[i]]$chr <= 22), ]
}

# Assess ploidy and recentre on a per sample basis
for ( i in 1:length(tracerx.list) ) {
  weights <- tracerx.list[[i]][,3] - tracerx.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(tracerx.list[[i]]) ) {
    ploidy <- round(sum(tracerx.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      tracerx.list[[i]][ ,j] <- tracerx.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}

# liftover to hg19
tracerx.hg19 <- list()
for (k in 1:length(tracerx.list)) {
  print(paste(k,"/",length(tracerx.list)))
  iddf <- tracerx.list[[k]][ ,c(1:3)]
  iddf$chr <- paste("chr",iddf$chr,sep="")
  idrange <- makeGRangesFromDataFrame(iddf)
  chain <- import.chain("~/Documents/CNA/Github/singleBiopsyITH/Data/hg38Tohg19.over.chain")
  hg19ids <- liftOver(idrange, chain)
  hg19ids <- as.data.frame(hg19ids)
  
  tracerx.hg19.df <- cbind(ID=1:nrow(iddf), iddf, data.frame(chr.hg19=NA, start.hg19=NA, stop.hg19=NA))
  
  for ( i in 1:nrow(tracerx.hg19.df) ) {
    wd <- hg19ids[which(hg19ids$group==tracerx.hg19.df$ID[i]),] # pull potential hg19 id to match the hg38 ID
    wd <- wd[which(wd$seqnames == tracerx.hg19.df$chr[i]), ] # keep only if on same chromosome
    
    if (nrow(wd) == 0) {
      tracerx.hg19.df$chr.hg19[i] <- NA
      tracerx.hg19.df$start.hg19[i] <- NA
      tracerx.hg19.df$stop.hg19[i] <- NA
    } else {
      tracerx.hg19.df$chr.hg19[i] <- paste("chr",sub('...','',wd$seqnames[1]),sep="")
      tracerx.hg19.df$start.hg19[i] <- min(wd$start)
      tracerx.hg19.df$stop.hg19[i] <- max(wd$end)
    }
  }
  tracerx.hg19[[k]] <- tracerx.list[[k]] %>% inner_join(tracerx.hg19.df, by = c("start","stop"))
}

# keep only hg19 coordinates
for (i in 1:length(tracerx.hg19)) {
  tracerx.hg19[[i]] <- tracerx.hg19[[i]][, c(7,8,9,4)]
  colnames(tracerx.hg19[[i]])[c(1:3)] <- c("chr","start","stop")
  tracerx.hg19[[i]] <- na.omit(tracerx.hg19[[i]])
  tracerx.hg19[[i]]$chr <- substr(tracerx.hg19[[i]]$chr, 4, nchar(tracerx.hg19[[i]]$chr))
}

# Bin CN data to match training dataset
tracerxBinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = tracerx.hg19)
x <- tracerxBinned.list[[1]]

# Pull samples from same patient into dataframe per patient
tracerxcn.list <- list()
for ( i in 1:length(patient) ) {
  x <- which(sub("\\-.*","",names) == patient[i])
  tracerxcn.list[[i]] <- tracerxBinned.list[x] %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
  tracerxcn.list[[i]] <- tracerxcn.list[[i]][-1]
}

# Convert data to numeric (removes characters)
tracerxcn.list <- lapply(tracerxcn.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))))
  x
})

# Drop patient without MR
drop <- which(sapply(tracerxcn.list, function(x) ncol(x) < 5) == TRUE)
tracerxcn.list <- tracerxcn.list[-drop]

# collapse into one dataframe
tracerx.raw <- tracerxcn.list %>% purrr::reduce(full_join, by = c("chr","start","stop"))

# If absolute copy number >=3, its a gain, if <=2 its a loss
tracerx.raw[ ,-c(1:3)] <- ifelse(tracerx.raw[ ,-c(1:3)] >= 3, 3, ifelse(tracerx.raw[ ,-c(1:3)] < 2, 1, 2))

# pull data info, clonality, diversity
tracerx.info <- PullDataInfo(rawdata = tracerx.raw) # hg38
tracerx.clonality <- PullDataClonality(rawdata = tracerx.raw, dataInfo = tracerx.info)
tracerx.diversity <- PullDataDiversity(rawdata = tracerx.raw, dataInfo = tracerx.info)
tracerx.actualITH <- data.frame(patient = lapply(data.frame(patient=tracerx.info$patientIDs), rep, tracerx.info$sampPerPatient),
                                sample = tracerx.info$sampleIDs,
                                actual = lapply(data.frame(actual=tracerx.diversity$pic.frac$pic.frac), rep, tracerx.info$sampPerPatient))

# pull matrix
tracerx.matrices <- genMatrices(tracerx.raw)

# Save
saveRDS(tracerx.raw, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.raw.rds")
saveRDS(tracerx.info, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.info.rds")
saveRDS(tracerx.clonality, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.clonality.rds")
saveRDS(tracerx.diversity, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.diversity.rds")
saveRDS(tracerx.actualITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.actualITH.rds")
saveRDS(tracerx.matrices, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.matrices.rds")
















# old
# create 1mb bins for hg38
hg38.length <- SeqinfoForUCSCGenome("hg38")
hg38.length <- data.frame("chr"=substr(hg38.length@seqnames,4,nchar(hg38.length@seqnames)), "length"=hg38.length@seqlengths)
hg38.length <- hg38.length[c(1:22),]
hg38.length$max <- plyr::round_any(hg38.length$length, 1000000, f = ceiling)

chr <- list()
start <- list()
stop <- list()
for ( i in 1:nrow(hg38.length) ) {
  max <- plyr::round_any(hg38.length$length[i], 1000000, f = ceiling)
  chr[[i]] <- rep(hg38.length$chr[i], max/1000000)
  start[[i]] <- seq(0 , max-1000000, 1000000)
  stop[[i]] <- seq(1000000 , max, 1000000)
}

tracerx.hg38binned <- cbind(bin=1:length(unlist(chr)), chr=unlist(chr), start=unlist(start), stop=unlist(stop))
tracerx.hg38binned <- data.frame(apply(tracerx.hg38binned, 2, function(x) as.numeric(as.character(x))))
