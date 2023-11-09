# PREREQUISITS: load section 1 data 

# Prepare data for heatmap showing genome-wide correlations between diploid/aneuploid
# Having two multiregion samples per patient will drive up the correlation. 
# Therefore, we will consider the bin-bin correlations using the average copy number per patient

# for crc
list <- list()
l <- 1
for ( k in 1:nrow(car.matrices$diploid.aneu) ) { 
    i <- 1
    
    while ( i < ncol(car.matrices$diploid.aneu) ) { 
      list[[l]] <- ((car.matrices$diploid.aneu[k,i]) + (car.matrices$diploid.aneu[k,i+1])) / 2 
      i <- i + 2
      l <- l + 1
    }
  }
cna.corr.car <- data.frame(t(matrix(unlist(list), ncol=car.info$noBins)), check.names = FALSE)

cna.corr.car <- t(cna.corr.car)
cna.corr.car <- cor(cna.corr.car, method = 'pearson', use = 'pairwise.complete.obs')
colnames(cna.corr.car) <- 1:ncol(cna.corr.car)
rownames(cna.corr.car) <- 1:nrow(cna.corr.car)

# plot crc
hm.car <- cnaCoocHM(corrMatrix=cna.corr.car, dataInfo=car.info)
hm1 <- cowplot::plot_grid(hm.car$legend,
                          annotate_figure(plot_grid(hm.car$plot, ncol = 1),
                                          left = text_grob('Chromosome', size = 28, rot = 90, vjust = 2),
                                          bottom = text_grob('Chromosome', size = 28, vjust = -1)),
                          ncol = 1, align = 'v', rel_heights = c(0.08, 1))
#jpeg('tempfig.jpeg', width = 1028.409, height = 1028.409)
jpeg('tempfig.jpeg', width = (40), height = (40), units = 'cm', res = 300)
hm1
dev.off()

# for lung
list <- list()
l <- 1
for ( k in 1:nrow(tracerx.matrices$diploid.aneu) ) { 
  i <- 1
  
  while ( i < ncol(tracerx.matrices$diploid.aneu) ) { 
    list[[l]] <- ((tracerx.matrices$diploid.aneu[k,i]) + (tracerx.matrices$diploid.aneu[k,i+1])) / 2 
    i <- i + 2
    l <- l + 1
  }
}
cna.corr.tracerx <- data.frame(t(matrix(unlist(list), ncol=tracerx.info$noBins)), check.names = FALSE)

cna.corr.tracerx <- t(cna.corr.tracerx)
cna.corr.tracerx <- cor(cna.corr.tracerx, method = 'pearson', use = 'pairwise.complete.obs')
colnames(cna.corr.tracerx) <- 1:ncol(cna.corr.tracerx)
rownames(cna.corr.tracerx) <- 1:nrow(cna.corr.tracerx)

# plot lung
hm.tracerx <- cnaCoocHM(corrMatrix=cna.corr.tracerx, dataInfo=tracerx.info)
hm2 <- cowplot::plot_grid(hm.tracerx$legend,
                          annotate_figure(plot_grid(hm.tracerx$plot, ncol = 1),
                                          left = text_grob('Chromosome', size = 28, rot = 90, vjust = 2),
                                          bottom = text_grob('Chromosome', size = 28, vjust = -1)),
                          ncol = 1, align = 'v', rel_heights = c(0.08, 1))
#jpeg('tempfig.jpeg', width = 1028.409, height = 1028.409)
jpeg('tempfig.jpeg', width = (40), height = (40), units = 'cm', res = 300)
hm2
dev.off()

# save
saveRDS(cna.corr.car, '~/Documents/CNA/Github/singleBiopsyITH/Data/cna.corr.car.rds') # needed for clustering
saveRDS(cna.corr.tracerx, '~/Documents/CNA/Github/singleBiopsyITH/Data/cna.corr.tracerx.rds') # needed for clustering



