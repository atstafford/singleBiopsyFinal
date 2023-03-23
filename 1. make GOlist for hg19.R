# make GO list for hg19
GOlist <- read.delim("~/Documents/CNA/Data/hg38.1_all_gene_GO_annotations.txt")
GOlist <- GOlist[GOlist$Chromosome.scaffold.name %in% paste(seq(1,22,1)), ]
colnames(GOlist)[c(5,3,4)] <- c("chr","start","stop")
GOlist$chr <- as.character(paste('chr',GOlist$chr,sep = "", collapse = NULL))

iddf <- GOlist[,c(5,3,4)]
iddf <- iddf[!duplicated(iddf),]

idrange <- GenomicRanges::makeGRangesFromDataFrame(iddf)
chain <- rtracklayer::import.chain("~/Documents/CNA/Github/singleBiopsyITH/Data/hg38Tohg19.over.chain")
hg19ids <- rtracklayer::liftOver(idrange, chain)
hg19ids <- as.data.frame(hg19ids)

GOlist.hg19 <- cbind(ID=1:nrow(iddf), iddf, data.frame(chr.hg19=NA, start.hg19=NA, stop.hg19=NA))

for ( i in 1:2 ) {
  print(paste(i, "/", nrow(GOlist.hg19)))
  wd <- hg19ids[which(hg19ids$group==GOlist.hg19$ID[i]),] # pull potential hg19 id to match the hg38 ID
  wd <- wd[which(wd$seqnames == GOlist.hg19$chr[i]), ] # keep only if on same chromosome
  
  if (nrow(wd) == 0) {
    GOlist.hg19$chr.hg19[i] <- NA
    GOlist.hg19$start.hg19[i] <- NA
    GOlist.hg19$stop.hg19[i] <- NA
  } else {
    GOlist.hg19$chr.hg19[i] <- sub('...','',wd$seqnames[1])
    GOlist.hg19$start.hg19[i] <- min(wd$start)
    GOlist.hg19$stop.hg19[i] <- max(wd$end)
  }
}

GOlist.hg19 <- GOlist %>% inner_join(GOlist.hg19)

#save
saveRDS(GOlist.hg19, "~/Documents/CNA/Github/singleBiopsyITH/Data/hg19.1_all_gene_GO_annotations.rds")