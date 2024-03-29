# PREREQUISITS: load section 1 data 

# First charcterise correlations between aberrations in individual bins and genome CNA diversity
# We will use 3 of the matrices: diploid.aneu (dip=0, aneu=1), diploid.gain (dip=0 gain=1), diploid.loss (dip=0 loss=1)
uniReg.in.list <- car.matrices[c(1,2,3)]

uniReg.in.list <- lapply(uniReg.in.list, function(x) {
  # Transpose to make a bin per column
  x <- t(x)
  
  # Convert to character and then a factor
  x <- data.frame(apply(x, 2, as.character), row.names = rownames(x), check.names = FALSE)
  x <- data.frame(lapply(x, factor, levels=c(0,1), labels=c('diploid','CNA')), row.names = rownames(x), check.names = FALSE)
  
  # # Add PIC.frac (stored in col 5 of ith df) to be predicted in first column
  x <- cbind(PIC.frac = car.actualITH$actual, x)
  x
})

# Run regression for each bin, and collect the coefficients for a gain, a loss, or noCNA at that bin

# Create list of 3 dataframes to store univariate regression output
uniReg.out.list <- list(data.frame(CNA='aneuploid', chr=car.info$start.stop$chr, bin=car.info$start.stop$bin, coeff=NA, pval=NA),
                        data.frame(CNA='gain', chr=car.info$start.stop$chr, bin=car.info$start.stop$bin, coeff=NA, pval=NA),
                        data.frame(CNA='loss', chr=car.info$start.stop$chr, bin=car.info$start.stop$bin, coeff=NA, pval=NA))
names(uniReg.out.list) = c('diploid.aneu', 'diploid.gain', 'diploid.loss')

# For each input matrix
for ( i in 1:length(uniReg.in.list) ) { 
  
  # Skipt first column as it holds ITH
  for ( k in 2:ncol(uniReg.in.list[[i]]) ) {
    
    # Creating a working dataframe holding ITH and the predictor bin
    data <- uniReg.in.list[[i]][,c(1,k)]
    data <- na.omit(data)
    colnames(data) <- c('ITH','predictor')
    
    # Set diploid as baseline
    data <- data %>%
      mutate(predictor = relevel(predictor, ref = 'diploid')) 
    
    # Regression can only be run when there are >1 factors present
    if ( length(unique(as.character(data$predictor))) == 1 ) { 
      next
    }
    
    # Run regression
    reg <- lm(ITH ~ ., data)
    
    # Extract the coefficient and pvalue
    coeffs <- data.frame(t(summary(reg)$coefficients[,1])) 
    names(coeffs) <- unlist(reg$xlevels)
    pvals <- data.frame(t(summary(reg)$coefficients[,4]))
    names(pvals) <- unlist(reg$xlevels)
    
    if ('CNA' %in% names(coeffs)) {
      uniReg.out.list[[i]][(uniReg.out.list[[i]]$bin==(k-1)),c(4:5)] <- c(coeffs$CNA,pvals$CNA)
    }
  }
}


# Add frequency of diploid, gain, loss

# For each regression output matrix
for ( i in 1:length(uniReg.out.list) ) {
  
  freq <- list()
  l <- 1
  
  # Pull the frequency of a diploid, gain, or loss occuring in at least one sample and divde by noPatient to get fraction
  for ( k in 1:nrow(uniReg.out.list[[i]]) ) {
    bin <- uniReg.out.list[[i]]$bin[k]
    status <- uniReg.out.list[[i]]$CNA[k]
    
    if ( status == 'aneuploid') {
      freq[[l]] <- (car.clonality$CNA.clo.counts$gain[bin] + car.clonality$CNA.clo.counts$loss[bin]) / car.info$noPatients
      l <- l + 1
    }
    else if ( status == 'gain') {
      freq[[l]] <- car.clonality$CNA.clo.counts$gain[bin] / car.info$noPatients
      l <- l + 1
    }
    else if ( status == 'loss') {
      freq[[l]] <- car.clonality$CNA.clo.counts$loss[bin] / car.info$noPatients
      l <- l + 1
    }
  }
  
  uniReg.out.list[[i]]$CNA.freq <- unlist(freq)
}


# Add adjusted p-value
uniReg.out.list <- lapply(uniReg.out.list, function(x) {
  x$sig <- 'insig'
  x$sig[which(x$pval<=0.05)] <- 'sig'
  x$pcorrect = p.adjust(x$pval, method = "BH")
  x$adjsig <- 'insig'
  x$adjsig[which(x$pcorrect<=0.05)] <- 'sig'
  x
})


# PLOT

# Store the 3 figures in a list (sig)
plot.list <- list()

for ( i in 1:length(uniReg.out.list) ) {
  # Set colours for: dark background, insig bars, sig bars, light background
  if ( i == 1 ) { col <- c("white",alpha("#81A88D", 0.2),"#81A88D","#CCCCCC"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Aneuploidy' } 
  if ( i == 2 ) { col <- c("white",alpha("#B40F20", 0.2),"#B40F20","#CCCCCC"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Gains only' } 
  if ( i == 3 ) { col <- c("white",alpha("#046C9A", 0.2),"#046C9A","#CCCCCC"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Losses only' }
  
  plot.list[[i]] <- ggplot(uniReg.out.list[[i]]) +
    geom_bar(aes(x=bin, y=coeff, fill=factor(sig)), stat="identity", width = 1) + 
    geom_rect(data=car.info$chr.end, aes(NULL,NULL,xmin=start,xmax=end,fill=col),ymin=-Inf,ymax=Inf,alpha=0.1) + 
    scale_fill_manual(values = col) + 
    
    geom_line(aes(x=bin, y=CNA.freq/3.33), group=1, color='#000033', linetype=1, size=0.7) + 
    geom_hline(yintercept = 0, size=0.3, color="black") + 
    
    ggtitle(tit) +
    scale_x_continuous(expand = c(0,0), name="chromosome", breaks=car.info$chr.mid, labels = labs) +
    scale_y_continuous(breaks = c(seq(-0.1,0.3,0.1)), limits = c(-0.14,0.34), 
                       sec.axis = sec_axis(trans = ~.*3.33)) +
    
    theme_custom() +
    theme(legend.position = "none", axis.title = element_blank()) 
}

# Store the 3 figures in a list (adjsig)
plot.list <- list()

for ( i in 1:length(uniReg.out.list) ) {
  # Set colours for: dark background, insig bars, sig bars, light background
  if ( i == 1 ) { col <- c("white",alpha("#81A88D", 0.2),"#81A88D","#CCCCCC"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Aneuploidy' } 
  if ( i == 2 ) { col <- c("white",alpha("#B40F20", 0.2),"#B40F20","#CCCCCC"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Gains only' } 
  if ( i == 3 ) { col <- c("white",alpha("#046C9A", 0.2),"#046C9A","#CCCCCC"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Losses only' }
  
  plot.list[[i]] <- ggplot(uniReg.out.list[[i]]) +
    geom_bar(aes(x=bin, y=coeff, fill=factor(adjsig)), stat="identity", width = 1) + 
    geom_rect(data=car.info$chr.end, aes(NULL,NULL,xmin=start,xmax=end,fill=col),ymin=-Inf,ymax=Inf,alpha=0.1) + 
    scale_fill_manual(values = col) + 
    
    geom_line(aes(x=bin, y=CNA.freq/3.33), group=1, color='#000033', linetype=1, size=0.7) + 
    geom_hline(yintercept = 0, size=0.3, color="black") + 
    
    ggtitle(tit) +
    scale_x_continuous(expand = c(0,0), name="chromosome", breaks=car.info$chr.mid, labels = labs) +
    scale_y_continuous(breaks = c(seq(-0.1,0.3,0.1)), limits = c(-0.14,0.34), 
                       sec.axis = sec_axis(trans = ~.*3.33)) +
    
    theme_custom() +
    theme(legend.position = "none", axis.title = element_blank()) 
}

# Prepare common legend 
coeff.legend <- as_ggplot(cowplot::get_legend(ggplot(uniReg.out.list$diploid.aneu) +
                                                geom_bar(aes(x=bin, y=coeff, fill=factor(sig)), stat="identity") + 
                                                geom_line(aes(x=bin, y=500*CNA.freq, colour='line')) + 
                                                scale_fill_manual(values = c(alpha("#81A88D", 0.2),"#81A88D"), labels = c('p>0.05', 'p\u22640.05')) +
                                                labs(colour = NULL, fill = NULL) +
                                                scale_color_manual(values = '#000033', labels = 'Freq.') + 
                                                
                                                guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1),
                                                       colour = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 4)) +
                                                
                                                theme_legend())) #or 0.8 for square


# Assemble plots
coeff.plot <- cowplot::plot_grid(coeff.legend,
                                 plot_grid(plot.list[[1]], plot.list[[2]], plot.list[[3]], ncol = 1, align = "v", rel_heights = c(1,1,1)),
                                 ncol = 1, align = 'v', rel_heights = c(0.09, 1))


coeff.plot <- cowplot::plot_grid(coeff.legend,
                                    annotate_figure(plot_grid(plot.list[[1]], plot.list[[2]], plot.list[[3]], nrow = 3, align = "v"),
                                                    left = text_grob('Coefficient', size = 28, rot = 90),
                                                    right = text_grob('CNA Frequency', size = 28, rot = 270),
                                                    bottom = text_grob('Chromosome', size = 28, vjust = -1)),
                                    ncol = 1, align = 'v', rel_heights = c(0.09, 1))


#jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*7.87))
jpeg('tempfig.jpeg', width = (40), height = (40), units = 'cm', res = 300)
coeff.plot
dev.off()

# Save
saveRDS(uniReg.out.list, "~/Documents/CNA/Github/singleBiopsyITH/Data/uniReg.out.list.rds")



