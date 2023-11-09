# PREREQUISITS: load section 1 data 

# create annotation layer of common CRC genes
annotation <- data.frame( 
  x = c(960,1253,1311,1479,1808,2035,2360,2387,2480),
  label = c('APC','EGFR','MET','MYC','CCND1','RB1','TP53','ERBB2','SMAD4'),
  y = NA)

# create plots
clofreq.ad <- cloFreqPlot(clonalityData = ad.clonality, dataInfo = ad.info, annotation = annotation, title="Adenoma", ylab = NULL, xlab = NULL)
clofreq.car <- cloFreqPlot(clonalityData = car.clonality, dataInfo = car.info, annotation = annotation, title="Carcinoma", ylab = NULL, xlab = NULL)
clofreq.tracerx <- cloFreqPlot(clonalityData = tracerx.clonality, dataInfo = tracerx.info, annotation = annotation, title="TracerX", ylab = NULL, xlab = NULL)

# Test if clonality is correlated with frequency
cor.test(clofreq.ad$CloFreqDF$freq, clofreq.ad$CloFreqDF$pcSubclonal, type = 'pearson')
cor.test(clofreq.car$CloFreqDF$freq, clofreq.car$CloFreqDF$pcSubclonal, type = 'pearson')
cor.test(clofreq.tracerx$CloFreqDF$freq, clofreq.tracerx$CloFreqDF$pcSubclonal, type = 'pearson')


# plot adenoma and carincoma together
ClonalityPlot1 <- cowplot::plot_grid(clofreq.ad$legend,
                                    annotate_figure(plot_grid(clofreq.ad$plot, clofreq.car$plot, ncol = 1),
                                                    left = text_grob('Fraction of patients with CNA', size = 28, rot = 90),
                                                    bottom = text_grob('Chromosome', size = 28, vjust = -1)),
                                    ncol = 1, align = 'v', rel_heights = c(0.08, 1))
#jpeg('tempfig.jpeg', width = 1575, height = 892.34)
jpeg('tempfig.jpeg', width = (72), height = (35), units = 'cm', res = 300)
ClonalityPlot1
dev.off()

# plot tracerx
ClonalityPlot2 <- cowplot::plot_grid(clofreq.tracerx$legend,
                                     annotate_figure(plot_grid(clofreq.tracerx$plot, ncol = 1),
                                                     left = text_grob('Fraction of patients with CNA', size = 28, rot = 90, hjust=0.4),
                                                     bottom = text_grob('Chromosome', size = 28, vjust = -1)),
                                     ncol = 1, align = 'v', rel_heights = c(0.16, 1))
#jpeg('tempfig.jpeg', width = 1575, height = 460.34)
jpeg('tempfig.jpeg', width = (72), height = (20), units = 'cm', res = 300)
ClonalityPlot2
dev.off()


# t-test to compare the number of clonal CNAs between adenomas and carcinomas, as a proportion of PGA
t.test(car.clonality$patientClo$clonal/car.clonality$patientClo$CNA, ad.clonality$patientClo$clonal/ad.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of subclonal CNAs between adenomas and carcinomas, as a proportion of PGA
t.test(car.clonality$patientClo$subclonal/car.clonality$patientClo$CNA, ad.clonality$patientClo$subclonal/ad.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of clonal and subclonal CNAs in adenoma patients
t.test(ad.clonality$patientClo$subclonal/ad.clonality$patientClo$CNA,ad.clonality$patientClo$clonal/ad.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of clonal and subclonal CNAs in carcinoma patients
t.test(car.clonality$patientClo$subclonal/car.clonality$patientClo$CNA,car.clonality$patientClo$clonal/car.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)


