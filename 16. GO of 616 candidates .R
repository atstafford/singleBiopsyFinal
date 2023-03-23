# PREREQUISITS: load section 1,8 data 

# Pull bins that are individually predictive
GO <- rbind(uniReg.out.list[[2]][which(uniReg.out.list[[2]]$sig=='sig'),],
            uniReg.out.list[[3]][which(uniReg.out.list[[3]]$sig=='sig'),])

# Gene anno function needs to have cols: bin | chr | start | stop
GO <- merge(GO, car.info$start.stop)
GO <- GO[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]

# load gene list
GOlist.hg19 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/hg19.1_all_gene_GO_annotations.rds")
GOlist.hg19 <- GOlist.hg19[c(8,11,13,14,15,1)]
GOlist.hg19 <- GOlist.hg19[!duplicated(GOlist.hg19),]

# annotate
betaregCandGenes <- gene_anno(GO, GOlist.hg19)
betaregCandGenes$start <- car.info$start.stop[match(betaregCandGenes$bin, car.info$start.stop$bin), 3]
betaregCandGenes$stop <- car.info$start.stop[match(betaregCandGenes$bin, car.info$start.stop$bin), 4]
betaregCandGenes <- betaregCandGenes[c(1,3,7,8,2,6,4,5)]

# Add back in the gain/loss status to the data. 
betaregCandGenes$CNA <- GO[match(betaregCandGenes$bin, GO$bin),5]
betaregCandGenes$CNA.freq <- GO[match(betaregCandGenes$bin, GO$bin),6]
betaregCandGenes$symbol <- mapIds(org.Hs.eg.db, keys = betaregCandGenes$geneID,
                                 column = c('SYMBOL'), keytype = 'ENSEMBL')

# count number of genes gained and lost
GOgain <- betaregCandGenes$entrezgene[which(betaregCandGenes$CNA=="gain")]
GOgain <- GOgain[!duplicated(GOgain)]
length(GOgain)

GOloss <- betaregCandGenes$entrezgene[which(betaregCandGenes$CNA=="loss")]
GOloss <- GOloss[!duplicated(GOloss)]
length(GOloss)

# GO analysis
GO_cand <- limma::goana(list(Up=GOgain,Down=GOloss))
GO_cand$p.Up.adj = p.adjust(GO_cand$P.Up, method = "fdr")
GO_cand$p.Down.adj = p.adjust(GO_cand$P.Down, method = "fdr")
GO_cand$p.abs <- apply(GO_cand[,c(8:9)], 1 , min)
GO_cand$adjsig <- 'insig'
GO_cand$adjsig[which(GO_cand$p.abs<=0.05)] <- 'sig'

# collapse into parent terms
data <- GO_cand[which(GO_cand$p.abs <= 0.05),]
data$ID <- rownames(data)
simMatrix_cand <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.abs), data$ID)
reducedTerms_cand <- rrvgo::reduceSimMatrix(simMatrix_cand, scores = scores, threshold=0.6, orgdb="org.Hs.eg.db")
parentGO_cand <- data.frame(value=tapply(reducedTerms_cand$score, reducedTerms_cand$parentTerm, FUN=sum))
parentGO_cand$term <- rownames(parentGO_cand)
parentGO_cand$fraction <- ( (parentGO_cand$value/sum(parentGO_cand$value))*100 )
parentGO_cand <- parentGO_cand %>% add_row(term = "other", value=NA, fraction = sum(parentGO_cand$fraction[which(parentGO_cand$fraction<=4)]))
parentGO_cand <- parentGO_cand[which(parentGO_cand$fraction>4),]

plot <- ggplot(parentGO_cand, aes(x=reorder(term, -fraction), y=fraction, fill=term)) +
  geom_bar(stat="identity", width=0.75) +
  scale_fill_brewer(palette="Dark2") +
  ylab("% of GO parent terms") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_custom() +
  theme(axis.text.x = element_text(size=28, angle=90, hjust=1, vjust=0.4),
        axis.title.x = element_blank(),
        legend.position = "none")

jpeg('tempfig.jpeg', width = (3*37.795*4.75), height = (3*37.795*7.7))
plot
dev.off()

#bootstrapping 90% of 616
GO <- rbind(uniReg.out.list[[2]][which(uniReg.out.list[[2]]$sig=='sig'),],
            uniReg.out.list[[3]][which(uniReg.out.list[[3]]$sig=='sig'),])
candGOboot1 <- list()

for (i in 1:100) {
  print(i)
  x <- sample_n(GO, 585)
  x <- merge(x, car.info$start.stop)
  x <- x[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]
  
  gene.anno <- gene_anno(x, GOlist.hg19)
  gene.anno$start <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 3]
  gene.anno$stop <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 4]
  gene.anno <- gene.anno[c(1,3,7,8,2,6,4,5)]
  
  gene.anno$CNA <- x[match(gene.anno$bin, x$bin),5]
  gene.anno$CNA.freq <- x[match(gene.anno$bin, x$bin),6]
  
  GOgain <- gene.anno$entrezgene[which(gene.anno$CNA=="gain")]
  GOgain <- GOgain[!duplicated(GOgain)]
  GOloss <- gene.anno$entrezgene[which(gene.anno$CNA=="loss")]
  GOloss <- GOloss[!duplicated(GOloss)]
  
  GOanno <- limma::goana(list(Up=GOgain,Down=GOloss))
  GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
  GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")
  GOanno$p.abs <- apply(GOanno[,c(8:9)], 1 , min)
  GOanno$adjsig <- 'insig'
  GOanno$adjsig[which(GOanno$p.abs<=0.05)] <- 'sig'
  GOanno <- GOanno[which(GOanno$adjsig=="sig"),]
  
  if (nrow(GOanno)==0) {
    next
  }
  
  x <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]))
  simMatrix_all <- rrvgo::calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
  scores <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]),  scores = GOanno$p.abs[which(GOanno$p.abs <= 0.05)])
  scores <- setNames(-log10(scores$scores), scores$ID)
  reducedTerms_all <- rrvgo::reduceSimMatrix(simMatrix_all, threshold=0.6, orgdb="org.Hs.eg.db", scores = scores)
  candGOboot1[[i]] <- aggregate(reducedTerms_all$score, by=list(Category=reducedTerms_all$parentTerm), FUN=sum)
}

# candGOboot1 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/candGOboot1.rds")
terms <- c("DNA packaging", "nucleosome assembly","DNA replication-dependent nucleosome organization","chromatin assembly",
           "chromatin organization involved in negative regulation of transcription","chromatin organization involved in regulation of transcription")
df <- data.frame(i=1:100, chromatin=0)
for (i in 1:length(candGOboot1)) {
  wd <- candGOboot1[[i]]
  wd <- slice_max(wd, order_by = wd$x, n=1)
  if (wd$Category[1] %in% terms) {
    df$chromatin[i] <- 1
  }
}
length(df$chromatin[which(df$chromatin==1)])/nrow(df)

#bootstrapping 616 random bins
GO <- rbind(uniReg.out.list[[2]])
GO <- merge(GO, car.info$start.stop)
GO <- GO[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]
gene.anno <- gene_anno(GO, GOlist.hg19)
gene.anno2 <- gene.anno
gene.anno$start <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 3]
gene.anno$stop <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 4]
gene.anno <- gene.anno[c(1,3,7,8,2,6,4,5)]
gene.anno$CNA <- GO[match(gene.anno$bin, GO$bin),5]
gene.anno$CNA.freq <- GO[match(gene.anno$bin, GO$bin),6]

candGOboot2 <- list()

for (i in 1:1000) {
  print(i)
  sample <- sample(1:2694, 616)
  x <- gene.anno[which(gene.anno$bin %in% sample),]
  
  GOgain <- gene.anno$entrezgene
  GOgain <- GOgain[!duplicated(GOgain)]
  
  GOanno <- limma::goana(list(Up=GOgain,Down=GOloss))
  GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
  GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")
  GOanno$p.abs <- apply(GOanno[,c(8:9)], 1 , min)
  GOanno$adjsig <- 'insig'
  GOanno$adjsig[which(GOanno$p.abs<=0.05)] <- 'sig'
  GOanno <- GOanno[which(GOanno$adjsig=="sig"),]
  
  if (nrow(GOanno)==0) {
    next
  }
  
  x <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]))
  simMatrix_all <- rrvgo::calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
  scores <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]),  scores = GOanno$p.abs[which(GOanno$p.abs <= 0.05)])
  scores <- setNames(-log10(scores$scores), scores$ID)
  reducedTerms_all <- rrvgo::reduceSimMatrix(simMatrix_all, threshold=0.6, orgdb="org.Hs.eg.db", scores = scores)
  
  candGOboot2[[i]] <- aggregate(reducedTerms_all$score, by=list(Category=reducedTerms_all$parentTerm), FUN=sum)
}


# candGOboot2 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/candGOboot2.rds")
terms <- c("DNA packaging", "nucleosome assembly","DNA replication-dependent nucleosome organization","chromatin assembly",
           "chromatin organization involved in negative regulation of transcription","chromatin organization involved in regulation of transcription")
df <- data.frame(i=1:786, chromatin=0)
for (i in 1:length(candGOboot2)) {
  wd <- candGOboot2[[i]]
  wd <- slice_max(wd, order_by = wd$x, n=5)
  wd <- wd$Category %in% terms
  if (TRUE %in% wd) { # does one for the chromatin terms appear in the top 5 of random 616 bins?
    df$chromatin[i] <- 1
  }
}

length(df$chromatin[which(df$chromatin==1)])/nrow(df)

# save
saveRDS(betaregCandGenes, "~/Documents/CNA/Github/singleBiopsyITH/Data/betaregCandGenes.rds")
saveRDS(candGOboot1, "~/Documents/CNA/Github/singleBiopsyITH/Data/candGOboot1.rds")
saveRDS(candGOboot2, "~/Documents/CNA/Github/singleBiopsyITH/Data/candGOboot2.rds")
