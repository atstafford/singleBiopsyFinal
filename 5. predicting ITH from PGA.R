# PREREQUISITS: load section 1 data 

# create plots
pgapred.ad <- predictITHfromPGA(diversityData = ad.diversity$pga,title="Adenoma", ylab = NULL, xlab = NULL, toplim = 0.6)
pgapred.car <- predictITHfromPGA(diversityData = car.diversity$pga, title="Carcinoma",  ylab = NULL, xlab = NULL, toplim = 0.6)
#remove <- rownames(tracerx.diversity$pga)[tracerx.diversity$pga$prop.aneu>0.8]
#pgapred.tracerx <- predictITHfromPGA(diversityData = tracerx.diversity$pga[rownames(tracerx.diversity$pga) %!in% remove,], ylab = NULL, xlab = NULL, toplim = 0.8)
pgapred.tracerx <- predictITHfromPGA(diversityData = tracerx.diversity$pga, ylab = NULL, xlab = NULL, toplim = 1)

summary(pgapred.ad$r2)
summary(pgapred.car$r2)
summary(pgapred.tracerx$r2)

# plot adenoma and carcinoma together
#pgaPredictPlot1 <- annotate_figure(plot_grid(pgapred.ad$plot, pgapred.car$plot, ncol=1), 
#                                   left = text_grob('CNA diversity', rot = 90, size = 28, vjust = 2),
#                                   bottom = text_grob('PGA', size = 28, vjust = 0))
pgaPredictPlot1 <- annotate_figure(plot_grid(pgapred.ad$plot, pgapred.car$plot, ncol=2), 
                                   left = text_grob('CNA diversity (CNAdf)', rot = 90, size = 28, vjust = 0.5),
                                   bottom = text_grob('PGA', size = 28, vjust = 0))
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*8.38))
jpeg('tempfig.jpeg', width = (40), height = (20), units = 'cm', res = 300)
pgaPredictPlot1
dev.off()

# plot tracerx 
pgaPredictPlot2 <- annotate_figure(plot_grid(pgapred.tracerx$plot, ncol=1), 
                                   left = text_grob('CNAdf', rot = 90, size = 28, vjust = 1),
                                   bottom = text_grob('PGA', size = 28, vjust = 0))
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
pgaPredictPlot2
dev.off()


# COAD and PREDICTED ITH
x <- cbind(coad.pga, coad.predictedITH)
ggplot(data=x, aes(y=prop.aneu, x=predicted)) +
  geom_smooth(method = "lm", formula = y ~ x, color="black", fill="#273046") +
  #geom_point(size=3, shape=21, color="black", fill=colourChoice) +
  geom_point(size=5, shape=21, color="black", fill=alpha("#273046", 0.9)) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(after_stat(adj.rr.label), "*\", \"*", after_stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 0, label.y = 'top', size=10, ) +
  ggtitle('TCGA (COAD)') +
  coord_fixed(ratio=1) +
  #scale_y_continuous(ylab, expand = c(0.01,0), limits = c(0,1)) +
  #scale_x_continuous(xlab, expand = c(0.01,0), limits = c(0,1), breaks = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2))) +
  theme_custom()


