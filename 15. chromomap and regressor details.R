# PREREQUISITS: load section 1,10,11 data 

hg19 <- makeHg19()
hg19$start <- 1
hg19 <- hg19[c(1,6,2,3)]
hg19$chrom <- as.character(as.numeric(hg19$chrom))
hg19 <- hg19[which(hg19$chrom %in% hg19predictors$chr),]

# make annotation file
annotation <- hg19predictors[,c(1:4,7)]
annotation$chr <- as.character((annotation$chr))
annotation$cna <- as.character((annotation$cna))
chromoMap::chromoMap(list(hg19), list(annotation),
          data_based_color_map = T, 
          data_type = "categorical",
          data_colors = list(c("#81A88D","#B40F20","#046C9A")),
          canvas_width = 600, canvas_height = 600,
          chr_color = c("#CCCCCC"),
          #chr_width = 13, chr_length = 6, ch_gap = -6, 
          chr_width = 15, chr_length = 4, ch_gap = 10, 
          #top_margin = 0,
          legend = T, lg_x = 50, lg_y = 150, text_font_size = 17 )

# plot coefficients with clonality
data <- data.frame(coeff=summary(betareg_39, type = "pearson")$coefficients[[1]][-1,1])
data$regressor <- substr(rownames(data),5,nchar(rownames(data)))
data$bin <- as.numeric(sub("\\_.*", "", data$regressor))
data$cna <- sub('.*_', '', data$regressor)

x <- car.clonality$pcSubclonal
x <- x[which(x$bin %in% hg19predictors$bin),]

list <- list()
for (i in 1:nrow(data)) {
  bin <- data$bin[i]
  if (data$cna[i] == "gain") {
    list[[i]] <- c(x$gain[which(x$bin==bin)])
  }
  if (data$cna[i] == "loss") {
    list[[i]] <- c(x$loss[which(x$bin==bin)])
  }
}

data$pcSubclonal <- unlist(list)

plot1 <- ggplot(data=data, aes(x=fct_reorder(regressor, coeff) , y=coeff, fill=coeff)) +
  geom_bar(stat="identity") +
  scale_fill_gradientn(colours = c("#000066","#3300CC","#6600FF","#CC66FF","#FF6699","#FF0066","#CC0033","#990000"), limits = c(-2,2)) +
  scale_y_continuous(name="Coefficient", breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5,2)) +
  coord_flip() +
  theme_custom() +
  theme(axis.title.y = element_blank(), 
        legend.position = "none", 
        plot.margin=margin(t=0.2,r=0,b=0.5,l=0.5,"cm"))

plot2 <- ggplot(data=data, aes(x=fct_reorder(regressor, coeff) , y=pcSubclonal)) +
  geom_bar(stat="identity", width = 0.5, fill=alpha("#330033",0.5)) +
  geom_text(aes(label=round(pcSubclonal,1)), nudge_y = 0.15, size=8) +
  scale_y_continuous(name="% subclonal", breaks=c(0,1)) +
  coord_flip() +
  theme_custom() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        plot.margin=margin(t=0.2,r=0.8,b=0.5,l=0.5,"cm"))

plot <- cowplot::plot_grid(plot1, plot2, ncol = 2, align = 'h', rel_widths =c(1, 0.3))

#jpeg('tempfig.jpeg', width = (3*37.795*5.92), height = (3*37.795*8))
jpeg('tempfig.jpeg', width = (33), height = (35), units = 'cm', res = 300)
plot
dev.off()

# clonality histogram
plot <- ggplot(z, aes(round(pcSubclonal,2))) + 
  geom_histogram(binwidth = 0.05, color="black", fill=alpha("#3399FF",0.5)) +
  xlab("% subclonal") +
  scale_x_continuous(breaks=seq(0,1,0.2)) +
  theme_custom()

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
plot
dev.off()

# what genes are in the clonal bin?
clonalBin <- data[data$pcSubclonal==0,3]
clonalBin <- car.info$start.stop[car.info$start.stop$bin==as.numeric(clonalBin),]

GOlist.hg19 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/hg19.1_all_gene_GO_annotations.rds")
GOlist.hg19 <- GOlist.hg19[c(8,11,13,14,15,1)]
GOlist.hg19 <- GOlist.hg19[!duplicated(GOlist.hg19),]

clonalBin <- gene_anno(clonalBin, GOlist.hg19)
clonalBin$symbol <- mapIds(org.Hs.eg.db, keys = clonalBin$geneID,
                                 column = c('SYMBOL'), keytype = 'ENSEMBL')





# tile of cna per predictors per patient in training
x <- car.raw[-c(1:3)]
rownames(x) <- 1:nrow(x)
x <- x[which(rownames(x) %in% hg19predictors$bin), ]
x$bin <- rownames(x)
x$chr <- hg19predictors[match(x$bin, hg19predictors$bin),2]

x$bin <- as.factor(x$bin)
y <- gather(x, key="Sample", value="cna", -bin, -chr)

plot <- ggplot(y, aes(x=bin, y=Sample, fill = as.factor(cna))) +
  geom_tile(color = "black", lwd = 0.1, linetype = 1) +
 scale_fill_manual(values=c("#046C9A","white","#B40F20")) +
  facet_grid(. ~ chr, scales = "free", space = "free") +
  theme(
    axis.text.x = element_text(size=14, angle=90, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.title = element_text(size=24),strip.text = element_text(size=24),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "none") +
  guides()

jpeg('tempfig.jpeg', width = 750, height = 750)
plot
dev.off()


