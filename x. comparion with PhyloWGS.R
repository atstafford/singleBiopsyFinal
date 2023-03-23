# load data
phyloWGS <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/Raynaud2018.xlsx", col_names = TRUE, skip =0)

# keep only CRC
phyloWGS <- phyloWGS[which(phyloWGS$`tumor type`=="CRC"),]
phyloWGS <- phyloWGS[which(phyloWGS$`tumor subtype`=="CIN"),]
phyloWGS$predictedITH <- coad.predictedITH[match(substr(phyloWGS$sample_name,1,nchar(phyloWGS$sample_name)-3), coad.predictedITH$TCGA_barcode), 3]


# plot
ggplot(phyloWGS, aes(x=predictedITH, y=`number of clones`)) +
  geom_smooth(method = "lm", formula = y ~ x, color="black", fill="black") +
  geom_point(shape=21, fill = "black", size = 3) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'left', label.y = 'top', size=10) +
  theme_custom()

cor.test(phyloWGS$`number of clones`, phyloWGS$predictedITH)
cor.test(phyloWGS$`ABSOLUTE purity`, phyloWGS$predictedITH)
cor.test(phyloWGS$`TCGA purity (curated)`, phyloWGS$predictedITH)
cor.test(phyloWGS$`ABSOLUTE purity`, phyloWGS$`number of clones`)