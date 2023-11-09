# PREREQUISITS: load section 1,11 data 

# Are the beta bins spread evenly across both samples? What is the difference when using sample 1 vs 2
data <- gather(car.predictedITH, 'CNA diversity', value, -sample, -patient)

data$range <- NA
for (i in 1:length(unique(data$patient))) {
  wd <- data[which(data$patient == unique(data$patient)[i]), ]
  data$range[which(data$patient == unique(data$patient)[i])] <- max(wd$value) - min(wd$value)
}
data$range <- as.factor(data$range)

plot <- ggplot(data, aes(x=fct_reorder(patient, value, mean) , y=value)) +
  geom_line(aes(group = patient, alpha=fct_rev(range)), size=5) +
  geom_point(aes(size =`CNA diversity`, shape =`CNA diversity`, fill=`CNA diversity`)) +
  scale_shape_manual(values=c(23, 21)) +
  scale_size_manual(values=c(6,4)) +
  scale_fill_manual(values=c("#273046","#9999CC")) +
  ylab("CNAdf - predicted") +
  xlab("Patient") +
  theme_custom() +
  #coord_flip() +
  theme(legend.position = "top", legend.text = element_text(size=28), legend.title = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, size=20)) +
        #axis.text.x = element_text(size=28), 
        #axis.text.y = element_text(size=20)) +
  guides(color = FALSE, size = FALSE, alpha = FALSE, fill = guide_legend(override.aes = list(size = 8)))

#jpeg('tempfig.jpeg', width = (3*37.795*13.89), height = (3*37.795*4.39))
jpeg('tempfig.jpeg', width = (40), height = (20), units = 'cm', res = 300)
plot
dev.off()

# difference between two samples
diff <- list()
error <- list()
for ( i in 1:length(unique(car.predictedITH$patient)) ) {
  patient <- unique(car.predictedITH$patient)[i]
  wd <- car.predictedITH[which(car.predictedITH$patient ==  patient),]
  diff[[i]] <- max(wd$predicted)-min(wd$predicted)
}
error <- car.predictedITH
error$betweenSampleDiff <- unlist(diff)
mean(unlist(diff))
median(unlist(diff))
max(unlist(diff))
min(unlist(diff))

# difference between actual and predicted
diff <- list()
for ( i in 1:nrow(error) ) {
  if ( error$actual[i] >= error$predicted[i] ) {
    diff[[i]] <- error$actual[i] - error$predicted[i]
  }
  else if ( error$actual[i] < error$predicted[i] ) {
    diff[[i]] <- error$predicted[i] - error$actual[i]
  }
}
error$sampleError <- unlist(diff)
mean(unlist(diff))
median(unlist(diff))

# difference between actual and predicted ave per patient
diff <- list()
meanSampleError <- list()
for ( k in 1:length(unique(error$patient)) ) {
  wd2 <- error[which(error$patient == unique(error$patient)[k]),]
  for ( i in 1:nrow(wd2) ) {
    if ( wd2$actual[i] >= wd2$predicted[i] ) {
      diff[[i]] <- wd2$actual[i] - wd2$predicted[i]
    }
    else if ( wd2$actual[i] < wd2$predicted[i] ) {
      diff[[i]] <- wd2$predicted[i] - wd2$actual[i]
    }
  }
  meanSampleError[[k]] <- mean(unlist(diff))
}
error$meanSampleError <- rep(unlist(meanSampleError), each=2)

# correlation between error and PGA
error$pga <- car.diversity$pga$prop.aneu
summary(lm(betweenSampleDiff ~ pga, data=error))
summary(lm(sampleError ~ pga, data=error))
summary(lm(meanSampleError ~ pga, data=error))

# correlation between error and actual CNA diversity
summary(lm(betweenSampleDiff ~ actual, data=error))
summary(lm(sampleError ~ actual, data=error))
summary(lm(meanSampleError ~ actual, data=error))

# potential cut off
beta <- betareg_39
cutoff <- list()
aic <- list()

wd <- car.predictedITH
for ( i in 1:70 ) {
  remove <- grep(pattern=max(wd$actual),x = wd$actual)
  wd <- wd[-remove,]
  cutoff[[i]] <- max(wd$actual)
  mod <- lm(actual ~ predicted, data = wd)
  aic[[i]] <- round(AIC(mod),2)
}

x <- data.frame(cutoff=unlist(cutoff), aic=unlist(aic))

plot <- ggplot(x, aes(x=cutoff, y=aic)) +
  geom_point(aes(fill = aic), shape=21, size=5, alpha=0.9) +
  scale_fill_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 0.24, linetype="longdash", color="blue") +
  xlab("CNAdf cutoff") +
  ylab("AIC") +
  theme(plot.margin=margin(t=0,r=0.5,b=0,l=0,"cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=28, colour='black'),
        axis.text = element_text(size=28, colour='black'),
        axis.ticks.length=unit(0.2, "cm"),
        legend.position = "none")

#jpeg('tempfig.jpeg', width = 400, height = 375)
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
plot
dev.off()




