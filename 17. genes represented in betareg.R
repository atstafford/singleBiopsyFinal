# PREREQUISITS: load section 1,7,8,9,11 data 

# load gene list
GOlist.hg19 <- readRDS("~/Documents/CNA/Github/Data/hg19.1_all_gene_GO_annotations.rds")
GOlist.hg19 <- GOlist.hg19[c(8,11,13,14,15,1)]
GOlist.hg19 <- GOlist.hg19[!duplicated(GOlist.hg19),]

# pull bins in all clusters that are represented in final model
clusters <- unique(hg19predictors$cluster)
bins <- c(sig.hclust$bin[which(sig.hclust$cluster %in% clusters)], #need to add in the ones missing from sig.hclust
          candidate.bins$bin[which(candidate.bins$cluster %in% (clusters[which(clusters %!in% sig.hclust$cluster)]))])

GO <- rbind(uniReg.out.list[[2]][which(uniReg.out.list[[2]]$sig=='sig'),],
            uniReg.out.list[[3]][which(uniReg.out.list[[3]]$sig=='sig'),])
GO <- GO[which(GO$bin %in% bins),]
GO <- merge(GO, car.info$start.stop)
GO <- GO[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]
GO$cluster <- candidate.bins[match(GO$bin, candidate.bins$bin), 11]

# annotate per cluster
gene.anno_perCluster <- list()
for (i in 1:length(unique(GO$cluster))) {
  print(i)
  cluster <- unique(GO$cluster)[i]
  wd <- GO[which(GO$cluster==cluster),]
  gene.anno_perCluster[[i]] <- gene_anno(wd, GOlist.hg19)
  gene.anno_perCluster[[i]]$start <- car.info$start.stop[match(gene.anno_perCluster[[i]]$bin, car.info$start.stop$bin), 3]
  gene.anno_perCluster[[i]]$stop <- car.info$start.stop[match(gene.anno_perCluster[[i]]$bin, car.info$start.stop$bin), 4]
  gene.anno_perCluster[[i]]$direction <- wd[match(gene.anno_perCluster[[i]]$bin, wd$bin), 7]
  gene.anno_perCluster[[i]] <- gene.anno_perCluster[[i]][c(1,3,7,8,2,6,4,5,9)]
  gene.anno_perCluster[[i]]$CNA <- wd[match(gene.anno_perCluster[[i]]$bin, wd$bin),5]
  gene.anno_perCluster[[i]]$CNA.freq <- wd[match(gene.anno_perCluster[[i]]$bin, wd$bin),6]
  gene.anno_perCluster[[i]]$cluster <- candidate.bins[match(gene.anno_perCluster[[i]]$bin, candidate.bins$bin), 11]
}

betaregRepGenes <- do.call("rbind",gene.anno_perCluster)
betaregRepGenes$symbol <- mapIds(org.Hs.eg.db, keys = betaregRepGenes$geneID,
                                 column = c('SYMBOL'), keytype = 'ENSEMBL')

#RBM43: RNA binding motif protein, modulates cyclin B1 exp.
#NMI: interacts with oncogenes N-myc and C-myc (and STAT)
#TNFAIP6: hyaluronan binding protein. Role in inflammation. induces by proinflam cytokines. role in tumour dev.
#MIR4773-2 / MIR4773-1: microRNAs
#RIF1: localise to aberrant telomers, maybe involved in DNA repair. promote tumour growth maybe by wnt
#NEB: encode component of cytoskeletal matrix 

# find no genes/bins per cluster
cluster <- list()
chr <- list()
no.bins <- list()
no.genes <- list()
for (i in 1:length(gene.anno_perCluster)) {
  wd <- gene.anno_perCluster[[i]]
  cluster[[i]] <- unique(wd$cluster)
  chr[[i]] <- unique(wd$chr)
  no.bins[[i]] <- length(unique(wd$bin))
  no.genes[[i]] <- length(unique(wd$entrezgene))
}
x <- data.frame(cluster=unlist(cluster), chr=unlist(chr), no.bins=unlist(no.bins), no.genes=unlist(no.genes))

sum(x$no.bins) #bins represented in model
sum(x$no.bins)/car.info$noBins # as a %

sum(x$no.genes) #genes represented in model
sum(x$no.genes)/20000 # as a %

mean(x$no.bins) #bins per cluster
min(x$no.bins)
max(x$no.bins)

mean(x$no.genes) #genes per cluster
min(x$no.genes)
max(x$no.genes)

# load gene lists
hallmarkGenes <- readRDS('~/Documents/CNA/Github//Data/GeneLists/hallmarkGenes.rds')
cosmic <- read_csv('~/Documents/CNA/Github//Data/GeneLists/COSMIC 11_12_03 2022.csv')
chromatinGenes <- read.table('~/Documents/CNA/Github//Data/GeneLists/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Github//Data/GeneLists/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Github//Data/GeneLists/histoneGenes', header = T, sep = "\t")

# for each cluster, pull genes occurring from gene lists
GOI.percluster <- list()
for (i in 1:length(gene.anno_perCluster)) { 
  x <- GOlist.hg19$HGNC.symbol[which(GOlist.hg19$NCBI.gene..formerly.Entrezgene..ID %in% gene.anno_perCluster[[i]]$entrezgene)]
  
  GOI.percluster[[i]] <- list(hallmark=x[x %in% hallmarkGenes$gene_symbol],
                       cosmic=x[x %in% cosmic$`Gene Symbol`],
                       chromatin=x[x %in% chromatinGenes$Symbol],
                       DNAmod=x[x %in% DNAmodGenes$Symbol],
                       histone=x[x %in% histoneGenes$Symbol])
}

# across all clusters, pull genes occurring from gene lists
GOI.acrossclusters <- GOlist.hg19$HGNC.symbol[which(GOlist.hg19$NCBI.gene..formerly.Entrezgene..ID %in% betaregRepGenes$entrezgene)]
GOI.acrossclusters <- list(hallmark = GOI.acrossclusters[GOI.acrossclusters %in% hallmarkGenes$gene_symbol], 
                    cosmicTier1 = GOI.acrossclusters[GOI.acrossclusters %in% cosmic$`Gene Symbol`[which(cosmic$Tier==1)]],
                    cosmicTier2 = GOI.acrossclusters[GOI.acrossclusters %in% cosmic$`Gene Symbol`[which(cosmic$Tier==2)]],
                    chromatin = GOI.acrossclusters[GOI.acrossclusters %in% chromatinGenes$Symbol],
                    DNAmod = GOI.acrossclusters[GOI.acrossclusters %in% DNAmodGenes$Symbol],
                    histone = GOI.acrossclusters[GOI.acrossclusters %in% histoneGenes$Symbol])

# print list of represented genes 
paste(sort(unique(c(GOI.acrossclusters$hallmark))), sep = ',', collapse = ", ")
paste(sort(unique(c(GOI.acrossclusters$histone, GOI.acrossclusters$DNAmod, GOI.acrossclusters$chromatin))), sep = ',', collapse = ", ")
paste(sort(unique(c(GOI.acrossclusters$cosmicTier1))), sep = ',', collapse = ", ")
paste(sort(unique(c(GOI.acrossclusters$cosmicTier2))), sep = ',', collapse = ", ")

# GO of bins directly in model
GO <- hg19predictors
betaregGenes <- gene_anno(GO, GOlist.hg19)
betaregGenes$start <- car.info$start.stop[match(betaregGenes$bin, car.info$start.stop$bin), 3]
betaregGenes$stop <- car.info$start.stop[match(betaregGenes$bin, car.info$start.stop$bin), 4]
betaregGenes <- betaregGenes[c(1,3,7,8,2,6,4,5)]
betaregGenes$CNA <- GO[match(betaregGenes$bin, GO$bin),7]
betaregGenes$symbol <- mapIds(org.Hs.eg.db, keys = betaregGenes$geneID,
                                 column = c('SYMBOL'), keytype = 'ENSEMBL')

# save
saveRDS(betaregRepGenes, "~/Documents/CNA/Github/singleBiopsyITH/Data/betaregRepGenes.rds")
saveRDS(betaregGenes, "~/Documents/CNA/Github/singleBiopsyITH/Data/betaregGenes.rds")
saveRDS(GOI.acrossclusters, "~/Documents/CNA/Github/singleBiopsyITH/Data/GOI.acrossclusters.rds")
