
m_df = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
m_df$gs_name <- gsub('HALLMARK_(\\S+)','\\1',m_df$gs_name)
m_df <- m_df[,c(3:4)]
isc <- read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/merlossuarez_2011_ISC_signature.txt')[,1] # Include ISC signature
geneinfo <- read.delim('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/complete_gene_info.txt')
isc_entrez <- geneinfo[which(geneinfo$Symbol %in% isc),]
isc_entrez$gs_name <- "ISC"
isc_entrez <- isc_entrez[,c(12,9)]
colnames(isc_entrez) <- colnames(m_df)
wnt_entrez <- data.frame("WNT", read.table('~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/wnt_signalling_entrez.txt')[,1]) # Include extra WNT_signalling
colnames(wnt_entrez) <- colnames(m_df)
hallmarkGenes <- rbind(m_df, isc_entrez, wnt_entrez)

# save
saveRDS(hallmarkGenes, '~/Documents/CNA/Github/singleBiopsyITH/Data/GeneLists/hallmarkGenes.rds')