# PREREQUISITS: load section 10,11 data 

# download RNA-seq data
expquery <- TCGAbiolinks::GDCquery(project = "TCGA-COAD",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts")

TCGAbiolinks::GDCdownload(expquery, method = "api")
rnaseq <- TCGAbiolinks::GDCprepare(expquery)

# Set directory and pull relevent files
directory <- "~/Documents/CNA/Github/singleBiopsyITH/GDCdata/TCGA-COAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification"
sampleFiles <- grep("\\counts.tsv$",list.files(path = directory, recursive = TRUE),value=TRUE)

# Collect UUIDs and match to barcodes
UUID <- basename(dirname(sampleFiles)) # Extract UUID (the folder name for each patient)
UUID <- TCGAutils::UUIDtoBarcode(UUID, from_type = "file_id", legacy = F) # Collect patient barcode to align with TCGA clinical data
UUID$TCGA_barcode <- substr(UUID$associated_entities.entity_submitter_id, 1, 12) # Trim barcode to match with clinical
UUID <- UUID[-2]
colnames(UUID)[1] <- "UUID_RNA"

# pull RNA seq files
setwd("~/Documents/CNA/Github/singleBiopsyITH/GDCdata/TCGA-COAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification")
RNAfiles.list <- lst()
for (i in 1:length(sampleFiles)) {
  print(i)
  RNAfiles.list[[i]] <- read.table(file = sampleFiles[i], header = T, sep = "\t")[-c(1:4),]
}

# only use unstranded
unstrand_RNAfiles.list <- lst()
for (i in 1:length(RNAfiles.list)) {
  unstrand_RNAfiles.list[[i]] <- RNAfiles.list[[i]][ ,c(1,4)]
}

# add sample name
barcodes <- UUID[match(basename(dirname(sampleFiles)), UUID$UUID_RNA), 2] 
for (i in 1:length(unstrand_RNAfiles.list)) {
  colnames(unstrand_RNAfiles.list[[i]])[2] <- barcodes[i]
}

# merge
countData <- powerjoin::power_full_join(unstrand_RNAfiles.list, by = "gene_id")

# remove second RNAseq for duolicate patients
countData <- countData[ ,c(1, which(nchar(colnames(countData)) == 12))]




# load clinical data
COAD_clinical <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/COAD_clinical.rds")
metaData <- COAD_clinical[ ,which(colnames(COAD_clinical) %in% c("TCGA_barcode","stage","site_of_resection_or_biopsy","tissue_or_organ_of_origin","MSI","CMS"))]

# Add predicted ITH back into COAD_clinical
metaData$ITH <- coad.predictedITH[match(metaData$TCGA_barcode, coad.predictedITH$TCGA_barcode), 3]

# Keep only those samples whose ITH could be predicted
metaData <- metaData[!is.na(metaData$ITH),]

# Mark those patients that didnt have a CNA in any of the predictive bins
metaData$intercept <- "no"
metaData$intercept[which(metaData$ITH==coef(betareg_39)["(Intercept)"])] <- "yes"

# Exclude certain patients based on control requirements in DE analysis
#metaData <- metaData[which(metaData$MSI=="MSS"),]

# Add UUID_RNA
metaData$UUID_RNA <- UUID[match(metaData$TCGA_barcode, UUID$TCGA_barcode),1]  

# remove patients with no RNA data
metaData <- metaData[!is.na(metaData$UUID_RNA),]
metaData <- metaData[which(metaData$TCGA_barcode %in% colnames(countData)), ]

# take top25% and bottom 25%
metaData$divGroup4 <- ntile(metaData$ITH, 4)
metaData$divGroup3 <- "med"
metaData$divGroup3[which(metaData$divGroup4==4)] <- "top25"
metaData$divGroup3[which(metaData$divGroup4==1)] <- "bottom25"
#metaData_sub <- metaData[which(metaData$divGroup3!="med"),]

# make column order same as metadata
countData_sub <- countData[ ,c("gene_id", metaData$TCGA_barcode)]

# build dds objects
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData_sub, 
                                     colData = metaData, 
                                     design = ~ divGroup3, tidy = TRUE)

# run deseq
dde <- DESeq2::DESeq(dds)

# extract results
res <- DESeq2::results(dde, contrast=c("divGroup3","top25","bottom25"))
res <- res[order(res$padj),]
head(res)

# Make a basic volcano plot
jpeg('tempfig.jpeg', width = 400, height = 400)
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=NULL, xlim=c(-3,3), cex=1, cex.axis=2, cex.lab=2, cex.main=2))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# pull sig DE genes
# some of my diversity bins are rare implying many different mechanism for driving high ITH
# get changes in chromatin genes when keeping all MSI/MSS, doing top25 v bottom25 (remove med). dont control for anything
# up reg parent GOs are meaningless, down reg are immune
# lots of histone/chromatin genes DE
# may be able to include med
res.sig <- as.data.frame(res)
res.sig <- res.sig[which(res.sig$padj<0.05),]

res.sig$ensemble <- gsub("\\..*","", rownames(res.sig))
res.sig$gene_id <- mapIds(org.Hs.eg.db, keys = res.sig$ensemble,
       column = c('SYMBOL'), keytype = 'ENSEMBL')
res.sig$entrez <- mapIds(org.Hs.eg.db, keys = res.sig$ensemble,
                          column = c('ENTREZID'), keytype = 'ENSEMBL')
# annotate if in genelist
hallmarkGenes <- readRDS('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/hallmarkGenes.rds')
cosmic <- read_csv('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/COSMIC 11_12_03 2022.csv')
chromatinGenes <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/histoneGenes', header = T, sep = "\t")

res.sig$hallmark <- ifelse(res.sig$gene_id %in% hallmarkGenes$gene_symbol, TRUE, FALSE)
res.sig$cosmic <- ifelse(res.sig$gene_id %in% cosmic$`Gene Symbol`, TRUE, FALSE)
res.sig$chromatin <- ifelse(res.sig$gene_id %in% chromatinGenes$Symbol, TRUE, FALSE)
res.sig$DNAmodGenes <- ifelse(res.sig$gene_id %in% DNAmodGenes$Symbol, TRUE, FALSE)
res.sig$histoneGenes <- ifelse(res.sig$gene_id %in% histoneGenes$Symbol, TRUE, FALSE)

res.sig$gene_id[res.sig$cosmic==TRUE & res.sig$log2FoldChange>0] %in% x$`Gene Symbol`
res.sig$gene_id[res.sig$cosmic==TRUE & res.sig$log2FoldChange<0]

x <- res.sig$gene_id[res.sig$chromatin==TRUE & res.sig$log2FoldChange>0]
x <- res.sig$gene_id[res.sig$chromatin==TRUE & res.sig$log2FoldChange<0]

x <- res.sig$gene_id[res.sig$DNAmodGenes==TRUE & res.sig$log2FoldChange>0]
x <- res.sig$gene_id[res.sig$DNAmodGenes==TRUE & res.sig$log2FoldChange<0]

x <- res.sig$gene_id[res.sig$histoneGenes==TRUE & res.sig$log2FoldChange>0]
x <- res.sig$gene_id[res.sig$histoneGenes==TRUE & res.sig$log2FoldChange<0]

x <- cosmic[cosmic$Hallmark=="Yes",c(1,6,15,16)]
x <- x[!is.na(x$Hallmark),]

# overlap of DE genes and genes represented in model
overlap <- data.frame(gene_id=res.sig$gene_id[res.sig$ensemble %in% betaregRepGenes$geneID])
overlap$cosmic <- ifelse(overlap$gene_id %in% cosmic$`Gene Symbol`, TRUE, FALSE)
overlap$chromatin <- ifelse(overlap$gene_id %in% chromatinGenes$Symbol, TRUE, FALSE)
overlap$DNAmodGenes <- ifelse(overlap$gene_id %in% DNAmodGenes$Symbol, TRUE, FALSE)
overlap$histoneGenes <- ifelse(overlap$gene_id %in% histoneGenes$Symbol, TRUE, FALSE)

# GO of DE genes
Up <- res.sig[which(res.sig$log2FoldChange>0),]
Down <- res.sig[which(res.sig$log2FoldChange<0),]
GOanno <- limma::goana(list(Up=Up$entrez,Down=Down$entrez))

GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")

data <- GOanno[which(GOanno$p.Up.adj <= 0.05),]
data$ID <- rownames(data)
simMatrix_up <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.Up.adj), data$ID)
reducedTerms_up <- rrvgo::reduceSimMatrix(simMatrix_up, scores = scores, threshold=0.6, orgdb="org.Hs.eg.db")
parentGO_DEup <- data.frame(tapply(reducedTerms_up$score, reducedTerms_up$parentTerm, FUN=sum))

data <- GOanno[which(GOanno$p.Down.adj <= 0.05),]
data$ID <- rownames(data)
simMatrix_down <- rrvgo::calculateSimMatrix(data$ID, orgdb="org.Hs.eg.db", ont = c("BP"), method="Rel")
scores <- setNames(-log10(data$p.Down.adj), data$ID)
reducedTerms_down <- rrvgo::reduceSimMatrix(simMatrix_down, scores = scores, threshold=0.6, orgdb="org.Hs.eg.db")
parentGO_DEdown <- data.frame(tapply(reducedTerms_down$score, reducedTerms_down$parentTerm, FUN=sum))

# save
saveRDS(countData, "~/Documents/CNA/Github/singleBiopsyITH/Data/countData.rds")
saveRDS(metaData, "~/Documents/CNA/Github/singleBiopsyITH/Data/metaData.rds")
saveRDS(dds, "~/Documents/CNA/Github/singleBiopsyITH/Data/dds.rds")
saveRDS(dde, "~/Documents/CNA/Github/singleBiopsyITH/Data/dde.rds")
saveRDS(res, "~/Documents/CNA/Github/singleBiopsyITH/Data/res.rds")
saveRDS(res.sig, "~/Documents/CNA/Github/singleBiopsyITH/Data/res.sig.rds")
saveRDS(overlap, "~/Documents/CNA/Github/singleBiopsyITH/Data/overlap.rds")

saveRDS(countData_betaregGenes, "~/Documents/CNA/Github/singleBiopsyITH/Data/countData_betaregGenes.rds")
saveRDS(countData_betaregRepGenes, "~/Documents/CNA/Github/singleBiopsyITH/Data/countData_betaregRepGenes.rds")
saveRDS(countData_betaregCandGenes, "~/Documents/CNA/Github/singleBiopsyITH/Data/countData_betaregCandGenes.rds")






















# WIP
# calc heterogeneity of betaregGenes in high vs low ITH. 
# you can have high div because your pheno is controlled -> low div across patients because converging on fitness peak
# genes in mod/rep/cand/all have higher variability in expression in low diversity 
# high diversity tumours are more similar across patients with regards to RNA
x <- assay(vsdata)
x <- x[which(gsub("\\..*","",rownames(x)) %in% betaregRepGenes$geneID),]

stdev.per.gene_low <- data.frame(st.dev=rowSds(x[ ,which(colnames(x) %in% metaData$TCGA_barcode[metaData$divGroup3=="bottom25"])]))
stdev.per.gene_med <- data.frame(st.dev=rowSds(x[ ,which(colnames(x) %in% metaData$TCGA_barcode[metaData$divGroup3=="med"])]))
stdev.per.gene_high <- data.frame(st.dev=rowSds(x[ ,which(colnames(x) %in% metaData$TCGA_barcode[metaData$divGroup3=="top25"])]))
#rownames(stdev.per.gene_low) <- rownames(assay(vsdata))[which(gsub("\\..*","",rownames(x)) %in% betaregRepGenes$geneID)]
#rownames(stdev.per.gene_high) <- rownames(assay(vsdata))[which(gsub("\\..*","",rownames(x)) %in% betaregRepGenes$geneID)]
#rownames(stdev.per.gene_med) <- rownames(assay(vsdata))[which(gsub("\\..*","",rownames(x)) %in% betaregCandGenes$geneID)]
rownames(stdev.per.gene_low) <- rownames(assay(vsdata))
rownames(stdev.per.gene_med) <- rownames(assay(vsdata))
rownames(stdev.per.gene_high) <- rownames(assay(vsdata))

high <- cbind(stdev.per.gene_high, group="high")
med <- cbind(stdev.per.gene_med, group="med")
low <- cbind(stdev.per.gene_low, group="low")
x <- rbind(low, high)
x$group <- factor(x$group, levels = c("low","med", "high"))

ggplot(x, aes(x=group, y=st.dev, colour=group)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 0.5, position = position_jitter(0.15)) +
  geom_boxplot(width=1) +
  #stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 

wilcox.test(high$st.dev, med$st.dev, paired = TRUE, alternative = "two.sided")
wilcox.test(high$st.dev, low$st.dev, paired = TRUE, alternative = "two.sided")
wilcox.test(med$st.dev, low$st.dev, paired = TRUE, alternative = "two.sided")

z <- data.frame(cbind(gene=1:nrow(high), high=high$st.dev, med=med$st.dev, low=low$st.dev))
df_long <- z %>% gather(key = "diversity", value = "st.dev", high,med,low)
df_long %>% group_by(diversity) %>%  summarise(n = n(), mean = mean(st.dev), 
                                               sd = sd(st.dev))
df_long$diversity <- factor(df_long$diversity, levels = c("low","med","high"))
ggplot(df_long, aes(x = diversity, y = st.dev)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + theme(legend.position="top")

friedman.test(y = df_long$st.dev, groups = df_long$diversity, blocks = df_long$gene)

# bivariate plot of high vs low
# majority of points are under x=y meaning per gene, st dev tends to be higher in low div than high div
x <- assay(vsdata)
stdev.per.gene_low <- data.frame(st.dev=rowSds(x[ ,which(colnames(x) %in% metaData$TCGA_barcode[metaData$divGroup3=="bottom25"])]))
stdev.per.gene_high <- data.frame(st.dev=rowSds(x[ ,which(colnames(x) %in% metaData$TCGA_barcode[metaData$divGroup3=="top25"])]))
z <- data.frame(cbind(low=stdev.per.gene_low$st.dev, high=stdev.per.gene_high$st.dev))
z$color <- "A"
z$color[which(z$low>z$high)] <- "B"

sum(z$color=="A")
sum(z$color=="B")

wilcox.test(z$low, z$high, paired = TRUE, alternative = "two.sided")

plot <- ggplot(data = z, aes(x = low, y = high, color=as.factor(color))) +
  geom_point(size = 2, alpha=.2) +
  geom_abline(linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") +
  xlab("St.dev (low diversity patients)") +
  ylab("St.dev (high diversity patients)") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "none")

jpeg('tempfig.jpeg', width = 400, height = 375)
plot
dev.off()



# Remove rows with few fragments in total
nrow(ddsHTSeq)
keep <- rowSums(counts(ddsHTSeq)) > 1
dds <- ddsHTSeq[keep,]
nrow(dds)

# Set ITH low as the first level in condition factor 
dds$divGroup2 <- relevel(dds$divGroup2, 'low')
dds$divGroup4 <- relevel(dds$divGroup4, '1')

# Run differential expression pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Transform raw count to stabilize variance 
vsd <- vst(dds, blind=FALSE)

# PCA
plotPCA(vsd, intgroup = c("divGroup4")) + stat_ellipse()

# Volcano plot and add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
par(mfrow=c(1,1), mar=c(1,1,1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Subset results to padj<0.1 and sort by log2 fold to get sig gene with up/down reg
resSig <- res[ which(res$padj < 0.1 ), ]
head( resSig[ order( resSig$log2FoldChange ), ] )
tail( resSig[ order( resSig$log2FoldChange ), ] )
summary(res)

# Plot difference for a specific gene
plotCounts(dds, gene="ENSG00000120211.4", intgroup = "MSI")

# MA plot (blue show sig?)
# Weakly expressedd have no chance of differential expression because low counts suffer from high poisson noise
plotMA( res, ylim = c(-1, 1), colSig = "red" )

# Visualize dispersion estimates (within-group variability)
# Black points are dispersion estimates, and will flux around true values. 
# A trend line is fitter and each gene estimate is shrunk towards it. (blue points)
# Blue circles above main cloud have high gene-wise dispersion so werent shrunk
plotDispEsts( dds, ylim = c(1e-6, 1e1) )



# Add gene names
res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
res$ensembl <- sapply(strsplit(res$ensembl, ".", fixed=T), function(x) x[1]) # Remove the dot and numbers after in order to match
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl )
idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res,4)

# Visualize distance in heatmap
library("pheatmap")
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))

mat_col <- data.frame(condition = vsd$condition)
rownames(mat_col) <- colnames(sampleDistMatrix)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

rownames(sampleDistMatrix) <- rownames(mat_col)
rownames(sampleDistMatrix) <- rownames(sampleDistMatrix)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         annotation_row = mat_col,show_rownames = FALSE)
