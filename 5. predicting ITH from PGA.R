# PREREQUISITS: load section 1 data 

# create plots
pgapred.ad <- predictITHfromPGA(diversityData = ad.diversity$pga,title="Adenoma", ylab = NULL, xlab = NULL, toplim = 0.6)
pgapred.car <- predictITHfromPGA(diversityData = car.diversity$pga, title="Carcinoma",  ylab = NULL, xlab = NULL, toplim = 0.6)
remove <- rownames(tracerx.diversity$pga)[tracerx.diversity$pga$prop.aneu>0.8]
pgapred.tracerx <- predictITHfromPGA(diversityData = tracerx.diversity$pga[rownames(tracerx.diversity$pga) %!in% remove,], ylab = NULL, xlab = NULL, toplim = 0.8)

summary(pgapred.ad$r2)
summary(pgapred.car$r2)
summary(pgapred.tracerx$r2)

# plot adenoma and carcinoma together
pgaPredictPlot1 <- annotate_figure(plot_grid(pgapred.ad$plot, pgapred.car$plot, ncol=1), 
                                   left = text_grob('CNA diversity', rot = 90, size = 28, vjust = 2),
                                   bottom = text_grob('PGA', size = 28, vjust = 0))
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*8.38))
pgaPredictPlot1
dev.off()

# plot tracerx (remove outliers PGA>0.8)
pgaPredictPlot2 <- annotate_figure(plot_grid(pgapred.tracerx$plot, ncol=1), 
                                   left = text_grob('CNA diversity', rot = 90, size = 28, vjust = 2),
                                   bottom = text_grob('PGA', size = 28, vjust = 0))
jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
pgaPredictPlot2
dev.off()



