# PREREQUISITS: load section 1 data 

patients <- list()
number <- list()
pic.frac <- list()
l <- 1

# for a given patient p...
for (p in 1:validationCar.info$noPatients) {
  patient <- validationCar.info$patientIDs[p]
  wd <- validationCar.raw[ , which(sub("\\..*", "",colnames(validationCar.raw)) == patient)]
  
  # define combinations using d samples...
  for (d in 2:ncol(wd)) {
    combo <- combn(ncol(wd), d, simplify = T)
    
    # and use to calculate pic.frac
    for (c in 1:ncol(combo)) {
      patients[[l]] <- patient
      number[[l]] <- d
      cols <- combo[,c]
      pic <- 1 - ((rowSums(wd[,cols]==1, na.rm = TRUE)/d)^2 + 
                    (rowSums(wd[,cols]==2, na.rm = TRUE)/d)^2 + 
                    (rowSums(wd[,cols]==3, na.rm = TRUE)/d)^2)
      
      pic <- sum(pic, na.rm = TRUE)
      
      if ( d %% 3 == 0 ) {
        n1 <- n2 <- n3 <- d/3
      } else if ( d %% 3 == 1 ) {
        n1 <- ((d-1)/3) + 1
        n2 <- ((d-1)/3)
        n3 <- ((d-1)/3)
      } else if ( d %% 3 == 2 ) {
        n1 <- ((d-2)/3) + 1
        n2 <- ((d-2)/3) + 1
        n3 <- ((d-2)/3)
      }
      
      max.pic <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
      
      pic.frac[[l]] <- as.numeric(pic) / (max.pic*unlist(validationCar.info$binsPerPatient[p]))
      l <- l + 1
    }
  }
}

data <- data.frame(patients=unlist(patients), number=unlist(number), pic.frac=unlist(pic.frac))
data$patients <- factor(data$patients, levels = c(1,2,3,4,5,6,7,8,9,10,11))


plot <- ggplot(data = data, aes(x = number, y = pic.frac)) +
  geom_jitter(size=2.5, shape=21, width=0.2, height = 0.05,fill=alpha("#273046",0.8), color="black") +
  facet_grid(patients ~ ., space = 'free_y') +
  scale_x_continuous(name = "Number of samples used to calculate CNAdf", breaks = seq(1,13,1)) +
  scale_y_continuous(name = "CNAdf",  n.breaks = 3, limits = c(0,0.5)) +
  theme(
    plot.margin=margin(t=0.2,r=0.5,b=0.5,l=0.5,"cm"),
    plot.title = element_text(size=28, colour='black',face='bold', hjust=0),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size=28, colour='black'),
    axis.text.x = element_text(size=28, colour='black'),
    axis.text.y = element_text(size=20, colour='black'),
    strip.text.y = element_text(size=28),
    strip.background.y = element_rect(colour = "black", fill=alpha("#6699CC",0.2)),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = 'white', colour = NA, size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.ticks.length=unit(0.2, "cm"))

#jpeg('tempfig.jpeg', width = (3*37.795*8.78), height = (3*37.795*8.9))
jpeg('tempfig.jpeg', width = (40), height = (35), units = 'cm', res = 300)
plot
dev.off()
