# calculate
# 0=notCNA, 1=subclonalCNA, 2=clonalCNA
# keep only those bins which can be subclonal
# there are NAs in the tracerx.raw because they were binned into even bins (still hg38)
x <- car.clonality$clonal.data
y <- tracerx.clonality$clonal.data

x <- cbind(data.frame(bin=1:nrow(x)), x)
y <- cbind(data.frame(bin=1:nrow(y)), y)

#x <- x[rowSums(x[-1]==1, na.rm = T)>0, ]
#y <- y[rowSums(y[-1]==1, na.rm = T)>0, ]

a <- data.frame(group = "CRC", pcSubclonal = rowSums(x==1)/ncol(x))
b <- data.frame(group = "Lung", pcSubclonal = rowSums(y==1)/ncol(y))

z <- rbind(a,b)

# violin plot (point per bin), 
# y axis=number of patients with subclonal alteration as fraction of total number of patients
plot <- ggplot(z, aes(x=group, y=pcSubclonal, colour=group)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 1, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 

wilcox.test(a$pcSubclonal, b$pcSubclonal, paired = TRUE, alternative = "two.sided")



# Distribution of subclonal CNAs
# leave as per bin
IPHc <- function(dataClonality, dataInfo) { # need matrix with CN per patient (col) per bin (row)
  # two species: 1=subclonalCNA, 2clonalCNA
  # calculate pic per bin and sum across bins
  iph <- list()
  for (i in 1:nrow(dataClonality)) {
    d <- length(dataClonality[i,which(dataClonality[i,]!=0)])
    iph[i] <- 1 - ( ( length(dataClonality[i,which(dataClonality[i,]==1)]) / d ) ^2 + 
                      ( length(dataClonality[i,which(dataClonality[i,]==2)]) / d )^2  )
    
    if ( d %% 2 == 0 ) {
      n1 <- n2 <- d/2
    }
    if ( d %% 2 == 1 ) {
      n1 <- ((d-1)/2) + 1
      n2 <- ((d-1)/2)
    }
    
    max.pic <- 1 - ( (n1/d)^2 + (n2/d)^2  )
    iph[i] <- as.numeric(iph[i]) / max.pic
  }
  iph <- unlist(iph)
  return(iph)
}
x <- IPHc(car.clonality$clonal.data, car.info)
mean(x)

IPHc <- function(dataClonality) { # need matrix with CN per patient (col) per bin (row)
  #  i've figured out the calculation of how diverse the subclonal regions are - for the cohort you count 
  # the number of subclonal bins (your denominator and the only 'species') and per bin you count the number 
  # of patients that are subclonal at that bin. then it's the sum of those fractions squared that will give 
  # you the diversity of subclonal bins along the genome. a score of 1 means subclonal bins are localised to 
  # only one bin and a score close to 0 means they are very spread along the genome
  # two species: 1=subclonalCNA, 2clonalCNA
  # calculate pic per bin and sum across bins
  
  d <- length(rowSums(dataClonality==1) > 0)
  n <- rowSums(dataClonality==1, na.rm = T)
  
  iphc <- n/d
  iphc <- 1-sum(iphc^2)
  
  return(iphc)
}
a <- IPHc(car.clonality$clonal.data)
b <- IPHc(ad.clonality$clonal.data)
c <- IPHc(tracerx.clonality$clonal.data)



