# PREREQUISITS: load section 1 data 

# Define the maximum pic depending on the number of samples
max.pics <- list()
for ( d in 1:50 ) {
  if ( d %% 3 == 0 ) {
    n1 <- n2 <- n3 <- d/3
  }
  else if ( d %% 3 == 1 ) {
    n1 <- ((d-1)/3) + 1
    n2 <- ((d-1)/3)
    n3 <- ((d-1)/3)
  }
  else if ( d %% 3 == 2 ) {
    n1 <- ((d-2)/3) + 1
    n2 <- ((d-2)/3) + 1
    n3 <- ((d-2)/3)
  }
  max.pics[[d]] <- pic.score <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
}

# Define the maximum shannon depending on the number of samples
max.shan <- list()
for ( d in 1:50 ) {
  if ( d %% 3 == 0 ) {
    n1 <- n2 <- n3 <- d/3
  }
  else if ( d %% 3 == 1 ) {
    n1 <- ((d-1)/3) + 1
    n2 <- ((d-1)/3)
    n3 <- ((d-1)/3)
  }
  else if ( d %% 3 == 2 ) {
    n1 <- ((d-2)/3) + 1
    n2 <- ((d-2)/3) + 1
    n3 <- ((d-2)/3)
  }
  max.shan[[d]] <- -1 * ( ((n1/d)*SciViews::ln(n1/d)) + ((n2/d)*SciViews::ln(n2/d)) + ((n3/d)*SciViews::ln(n3/d)) )
}

# Calculate shannon using all samples for each patient
shan <- list()
shan.frac <- list()
for ( i in 1:validationCar.info$noPatients) {
  # Set working data as the cols holding the samples
  wd <- validationCar.raw[ ,which(sub("\\..*", "", colnames(validationCar.raw))==validationCar.info$patientIDs[i])]
  
  # Calculate freq of loss/dip/gain per bin
  x <- data.frame(loss=rowSums(wd == 1), diploid=rowSums(wd == 2), gain=rowSums(wd == 3))

  # Calculate Shannon div
  shan[[i]] <- vegan::diversity(x, index="shannon")
  shan[[i]] <- na.omit(shan[[i]])
  
  # Record the number of samples
  upto <- ncol(wd)
  
  # Store as a dataframe
  shan[[i]] <- data.frame(patient = i, shan = sum(shan[[i]], na.rm = TRUE), sampleNo = upto)
  
  # Define the max possible diversity given the number of sample, as maxPIC*number of bins
  max.ITH <- max.shan[[upto]] * nrow(x)
  
  # Store as a dataframe
  shan.frac[[i]] <- data.frame(patient = i, shan.frac = sum(shan[[i]], na.rm = TRUE)/max.ITH, sampleNo = upto)
  

  
}
testShan <- do.call('rbind',shan)
testShan.frac <- do.call('rbind',shan.frac)


# Calculate PIC using all samples for each patient
pic <- list()
pic.frac <- list()
for ( i in 1:validationCar.info$noPatients) {
  # Set working data as the cols holding the samples
  wd <- validationCar.raw[ ,which(sub("\\..*", "", colnames(validationCar.raw))==validationCar.info$patientIDs[i])]
  
  # Record the number of samples
  upto <- ncol(wd)
  
  # Use PIC function on wd
  pic[[i]] <- PIC(wd, upto, c(1:upto))
  pic[[i]] <- na.omit(pic[[i]])
  
  # Store as a dataframe
  pic[[i]] <- data.frame(patient = i, PIC = sum(pic[[i]], na.rm = TRUE), sampleNo = upto)
  
  # Define the max possible diversity given the number of sample, as maxPIC*number of bins
  max.ITH <- max.pics[[upto]] * nrow(x)
  
  # Store as a dataframe
  pic.frac[[i]] <- data.frame(patient = i, PIC.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH, sampleNo = upto)
  
}

testPIC <- do.call('rbind',pic)
testPIC.frac <- do.call('rbind',pic.frac)

scores <- cbind(testShan, testShan.frac[2], testPIC[2], testPIC.frac[2])
scores$PIC.frac<- scores$PIC.frac*800
scores$shan.frac<- scores$shan.frac*800

scores <- gather(scores, Metric, count, -patient, -sampleNo)
scores$Value <- "Raw"
scores$Value[which(scores$Metric=="shan.frac" | scores$Metric=="PIC.frac")] <- "As fraction" 
scores$Value <- as.factor(scores$Value)
scores$Value <- factor(scores$Value, levels = c("Raw","As fraction"))
scores$Method <- "PIC"
scores$Method[which(scores$Metric=="shan.frac" | scores$Metric=="shan")] <- "Shannon" 

labels <- list()
for (i in 1:validationCar.info$noPatients) {
  s <- unique(scores$sampleNo[which(scores$patient==validationCar.info$patientIDs[i])])
  labels[[i]] <- paste(i," (n=",s,")", sep = "", collapse = "")
}
labels <- unlist(labels)

plot <- ggplot(data = scores, aes(x=reorder(patient, count) , y=count, group=Metric)) +
  geom_line() +
  geom_point(size=5, aes(color=Method, shape=Value)) +
  scale_color_manual(values=c("#F98400", "#660066")) +
  scale_y_continuous(name="Raw diversity values, summed over bins", 
                     sec.axis = sec_axis(trans = ~./800, name="Values as fraction of maximum diversity")) +
  scale_x_discrete(name="Patient (sample number)", labels=labels)+
  scale_linetype_manual(values=c("solid","dotted")) +
  theme_custom() +
  theme(axis.text.y = element_text(size=28, colour='black', hjust=1),
        axis.title.y.right = element_text(size=28, colour='black', vjust=1, hjust=0.2),
        axis.title.y.left = element_text(size=28, colour='black', vjust=1, hjust=0.8),
        axis.text.x = element_text(size=28, colour='black', angle = 90, hjust=1),
        legend.position = "top",
        legend.text = element_text(size=28, colour='black'),
        legend.title = element_blank())  

#jpeg('tempfig.jpeg', width = (3*37.795*5.94), height = (3*37.795*5.64))
jpeg('tempfig.jpeg', width = (28), height = (20), units = 'cm', res = 300)
plot
dev.off()



