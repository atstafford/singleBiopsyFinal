# PREREQUISITS: load section 1 data 

# colours
palette.training <- c("#660066")
palette.validation <- c("#E6AB02")
palette.lung <- c("#F98400")

# violin plots showing PGA and pic.frac between training and validation
data <- data.frame(rbind(cbind(value=validationCar.diversity$pga$prop.aneu, cohort="Validation", variable="PGA"),
                      cbind(value=validationCar.diversity$pga$pic.frac, cohort="Validation", variable="CNAdiversity"),
                      cbind(value=car.diversity$pga$prop.aneu, cohort="Training", variable="PGA"),
                      cbind(value=car.diversity$pga$pic.frac, cohort="Training", variable="CNAdiversity")))

pvPGA <- round(wilcox.test(as.numeric(car.diversity$pga$prop.aneu), as.numeric(validationCar.diversity$pga$prop.aneu), alternative = "two.sided")$p.value,3)
pvPIC <- round(wilcox.test(as.numeric(car.diversity$pga$pic.frac), as.numeric(validationCar.diversity$pga$pic.frac), alternative = "two.sided")$p.value,3)

trainVval <- ggplot(data, aes(x = variable, y = as.numeric(value), fill = cohort)) +
  geom_split_violin(width=1.2, alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE) +
  annotate('text', y=c(0.7,1.1), x=c(0.70,1.70), 
           label=c(paste("p", "==", pvPIC, sep = ""), paste("p", "==", pvPGA, sep = "")), size=10, parse=TRUE, hjust = 0) +
  annotate('segment', x = c(0.8,1.8), xend = c(1.2,2.1), y = c(0.65,1.05), yend = c(0.65,1.05)) +
  scale_x_discrete(name = NULL, labels = c("CNA diversity", "PGA")) +
  scale_y_continuous(name = "Proportion", breaks = seq(0, 1, 0.2), limits = c(0, 1.1)) +
  scale_fill_manual(values= c(palette.training, palette.validation), name = "Cohort") +
  theme_custom() +
  theme(legend.text = element_text(size=28),
        legend.title = element_blank(),
        legend.position = "top")
  
jpeg('tempfig.jpeg', width = (3*37.795*5.33), height = (3*37.795*4.29))
trainVval
dev.off()

t.test(car.diversity$pga$prop.aneu, validationCar.diversity$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)
t.test(car.diversity$pga$pic.frac, validationCar.diversity$pga$pic.frac,
       alternative = 'two.sided', var.equal = FALSE)


# violin plots showing PGA and pic.frac between CRC and lung
data <- data.frame(rbind(cbind(value=tracerx.diversity$pga$prop.aneu, cohort="Lung", variable="PGA"),
                         cbind(value=tracerx.diversity$pga$pic.frac, cohort="Lung", variable="CNAdiversity"),
                         cbind(value=car.diversity$pga$prop.aneu, cohort="CRC", variable="PGA"),
                         cbind(value=car.diversity$pga$pic.frac, cohort="CRC", variable="CNAdiversity"),
                         cbind(value=validationCar.diversity$pga$prop.aneu, cohort="CRC", variable="PGA"),
                         cbind(value=validationCar.diversity$pga$pic.frac, cohort="CRC", variable="CNAdiversity")))

pvPGA <- format(wilcox.test(as.numeric(car.diversity$pga$prop.aneu), as.numeric(tracerx.diversity$pga$prop.aneu), alternative = "two.sided")$p.value, scientific = TRUE, digits = 3)
pvPIC <- format(wilcox.test(as.numeric(car.diversity$pga$pic.frac), as.numeric(tracerx.diversity$pga$pic.frac), alternative = "two.sided")$p.value, scientific = TRUE, digits = 3)

trainVlung <- ggplot(data, aes(x = variable, y = as.numeric(value), fill = cohort)) +
  geom_split_violin(width=1.2, alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .8,  show.legend = FALSE) +
  annotate('text', y=c(0.85,1.2), x=c(0.70,1.70), 
           label=c(paste("p", "==", pvPIC, sep = ""), paste("p", "==", pvPGA, sep = "")), size=10, parse=TRUE, hjust = 0) +
  annotate('segment', x = c(0.8,1.8), xend = c(1.2,2.1), y = c(0.8,1.15), yend = c(0.8,1.15)) +
  scale_x_discrete(name = NULL, labels = c("CNA diversity", "PGA")) +
  scale_y_continuous(name = "Proportion", breaks = seq(0, 1, 0.2), limits = c(0, 1.21)) +
  scale_fill_manual(values= c(palette.training, palette.lung), name = "Cohort") +
  theme_custom() +
  theme(legend.text = element_text(size=28),
        legend.title = element_blank(),
        legend.position = "top")

t.test(car.diversity$pga$prop.aneu, tracerx.diversity$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)
t.test(car.diversity$pga$pic.frac, tracerx.diversity$pga$pic.frac,
       alternative = 'two.sided', var.equal = FALSE)

jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.29))
trainVlung
dev.off()

# compare adenoma and carcinoma PGA
t.test(car.diversity$pga$prop.aneu, ad.diversity$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)

# compare adenoma and carcinoma pic.frac
t.test(car.diversity$pic.frac$pic.frac, ad.diversity$pic.frac$pic.frac,
       alternative = 'two.sided', var.equal = FALSE)

# compare IPH
carIPH <- IPH(dataMatrix = car.raw[ ,-c(1:3)], dataInfo = car.info, nboot = 1000)
adIPH <- IPH(dataMatrix = ad.raw[ ,-c(1:3)], dataInfo = ad.info, nboot = 1000)
tracerxIPH <- IPH(dataMatrix = tracerx.raw[ ,-c(1:3)], dataInfo = tracerx.info, nboot = 1000)

# compare distribution of subclonal CNA
times <- c(rowSums(car.clonality$clonal.data==1, na.rm = T))
bins <- c(as.numeric(rownames(car.clonality$clonal.data)))
x <- rep(bins, times = times)
hist(x)
ks.test(x, "punif",0, max(bins))

times <- c(rowSums(ad.clonality$clonal.data==1, na.rm = T))
bins <- c(as.numeric(rownames(ad.clonality$clonal.data)))
x <- rep(bins, times = times)
hist(x)
ks.test(x, "punif",0, max(bins))

times <- c(rowSums(tracerx.clonality$clonal.data==1, na.rm = T))
bins <- c(as.numeric(rownames(tracerx.clonality$clonal.data)))
x <- rep(bins, times = times)
hist(x)
ks.test(x, "punif",0, max(bins))
