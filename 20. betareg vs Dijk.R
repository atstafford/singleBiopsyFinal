# Dijk used absolute meanCN calls, the val cohort used sequenza for CN calls
# Working with absolute CN data, without ploidy recentering
# Take average across MR samples per bin to generate absolute meanCN 

# Read in the datasets
test_car1 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.01.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car2 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car3 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car4 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.04.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car5 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car6 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.06.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car7 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.07.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car8 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.08.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9p <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.09.Proximal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9d <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.09.Distal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car10 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.10.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

# Load patient datasets into a list
validationCar.list <- list(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)
rm(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)

# Convert data to numeric
validationCar.list <- lapply(validationCar.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})
x <- validationCar.list[[1]]

# Generate absolute copy number (minor + major)
for ( j in 1:length(validationCar.list)) {
  i <- 5
  l <- 1
  abscn <- list()
  
  while ( i < ncol(validationCar.list[[j]]) ) {
    
    for ( k in 1:nrow(validationCar.list[[j]]) ) {
      abscn[[l]] <- sum(validationCar.list[[j]][k,i], validationCar.list[[j]][k,i+1])
      l <- l + 1
    }
    i <- i + 2
  }
  
  validationCar.list[[j]] <- cbind(validationCar.list[[j]][,c(1,2,4)], matrix(unlist(abscn), nrow = nrow(validationCar.list[[j]])) )
  colnames(validationCar.list[[j]])[c(1:3)] <- c('chr','start','stop')
  colnames(validationCar.list[[j]])[c(4:ncol(validationCar.list[[j]]))] <- paste(j, seq(1, ncol(validationCar.list[[j]])-3, 1), sep = ".")
  validationCar.list[[j]] <- data.frame(validationCar.list[[j]])
}

validation_dijk.list <- lapply(validationCar.list, function(x) {
  abs.meanCN <- rowMeans(x[,-c(1:3)], na.rm = TRUE)
  ploidy <- mean(colMeans(x[,-c(1:3)], na.rm = TRUE))
  x <- cbind(x[,c(1:3)], ploidy, abs.meanCN)
})

dijk.output <- list()
for (i in 1:length(validation_dijk.list)) {
  seg_val <- validation_dijk.list[[i]]$abs.meanCN
  seg_len <- validation_dijk.list[[i]]$stop-validation_dijk.list[[i]]$start
  ploidy <- validation_dijk.list[[i]]$ploidy[1]
  dijk.output[[i]] <- dijk(seg_val,seg_len, ploidy, purity = NULL)
}

validation_dijk <- data.frame(matrix(unlist(dijk.output), nrow=length(dijk.output), byrow=TRUE))
colnames(validation_dijk) <- c('CNH', 'purity', 'ploidy')
validation_dijk$patient <- validationCar.info$patientIDs

# And for training
chr <- sub("\\:.*", "", rownames(car.raw))
start <- sub(".*[:]([^.]+)[-].*", "\\1", rownames(car.raw))
stop <- sub(".*-", "", rownames(car.raw))
names <- unique(colnames(car.raw))
patient <- sub("\\_.*", "", names)
sample <- sub("_", ".", names)
type <- as.factor(sub("^([[:alpha:]]*).*", "\\1", sub(".*_", "", names)))

training_dijk.list <- list()
for ( i in 1:car.info$noPatients ) {
  training_dijk.list[[i]] <- data.frame(car.info$start.stop, car.raw[,which(sub("\\..*", "", colnames(car.raw))==car.info$patientIDs[i])])
}

training_dijk.list <- lapply(training_dijk.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x <- as.data.frame(x)
  x
})

training_dijk.list <- lapply(training_dijk.list, function(x) {
  abs.meanCN <- rowMeans(x[,-c(1:4)], na.rm = TRUE)
  ploidy <- mean(colMeans(x[,-c(1:4)], na.rm = TRUE))
  x <- cbind(x[,c(2:4)], ploidy, abs.meanCN)
})

dijk.output <- list()
for (i in 1:length(training_dijk.list)) {
  seg_val <- training_dijk.list[[i]]$abs.meanCN
  seg_len <- training_dijk.list[[i]]$stop-training_dijk.list[[i]]$start
  ploidy <- training_dijk.list[[i]]$ploidy[1]
  dijk.output[[i]] <- dijk(seg_val,seg_len, ploidy, purity=NULL)
}

training_dijk <- data.frame(matrix(unlist(dijk.output), nrow=length(dijk.output), byrow=TRUE))
colnames(training_dijk) <- c('CNH', 'purity', 'ploidy')
training_dijk$patient <- car.info$patientIDs 

# Does Dijk predict our actual (PIC.frac)
DIJKvPICFRAC <- rbind(data.frame(pic.frac=validationCar.diversity$pic.frac$pic.frac, CNH=validation_dijk$CNH, patient=validation_dijk$patient, type="test"),
                      data.frame(pic.frac=car.diversity$pic.frac$pic.frac, CNH=training_dijk$CNH, patient=training_dijk$patient, type="train"))

plot <- ggplot(data = DIJKvPICFRAC[which(DIJKvPICFRAC$type=="train"),], aes(y = pic.frac, x = CNH)) +
  geom_smooth(method = 'lm', fill="#273046", color="black") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=10) +
  #geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(size = 5, fill=alpha("#273046",0.9), color="black", shape=21) +
  scale_y_continuous(limits = c(0,0.55), breaks = seq(0,0.4,0.1)) +
  #geom_text_repel(aes(label=patient, size=5), min.segment.length = 0, box.padding = 0.8, nudge_y = 0.001) +
  #scale_color_manual(values = c('#993366',"#333366","#663366")) + 
  ylab("CNA diversity (CNAdf)") +
  xlab("CNA diversity (Dijk's CNH)") +
  theme_custom()

#jpeg('tempfig.jpeg', width = (3*37.795*5.94), height = (3*37.795*5.64))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
plot
dev.off()

# compare Dijk methods with my prediction
MINEvDIJK <- rbind(car.predictedITH, validationCar.predictedITH)
MINEvDIJK$CNH <- "NA"
x <- rbind(training_dijk, validation_dijk)
MINEvDIJK$CNH <- x[match(MINEvDIJK$patient, x$patient), 1]

plot <- ggplot(data = MINEvDIJK, aes(x = CNH, y = predicted, color=type)) +
  geom_smooth(method = 'lm') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=8) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(size = 6) +
  scale_color_manual(values = c('#993366',"#333366","#663366")) + 
  ylab("Predicted using my method") +
  xlab("Dijk's CNH") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "top",
        legend.text = element_text(size=24, colour='black'),
        legend.title = element_blank() ) 

jpeg('tempfig.jpeg', width = 1000, height = 1000, pointsize = 12)
plot
dev.off()

# compare with my prediction averaged per sample
meanITH <- list()
CNH <- list()
type <- list()
for (i in 1:length(unique(MINEvDIJK$patient)) ) {
  patient <- unique(MINEvDIJK$patient)[i]
  wd <- MINEvDIJK$predicted[which(MINEvDIJK$patient==patient)]
  meanITH[[i]] <-  mean(wd, na.rm = TRUE)
  CNH[[i]] <-  unique(MINEvDIJK$CNH[which(MINEvDIJK$patient==patient)])
  type[[i]] <- unique(MINEvDIJK$type[which(MINEvDIJK$patient==patient)])
}

meanMINEvDIJK <- data.frame("patient"=unique(MINEvDIJK$patient), "CNH"=unlist(CNH), 'meanITH'=unlist(meanITH), "type"=unlist(type))

plot <- ggplot(data = meanMINEvDIJK, aes(y = CNH, x = meanITH, color=type)) +
  geom_smooth(method = 'lm') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=8) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(size = 6) +
  #geom_text_repel(aes(label=patient, size=5), min.segment.length = 0, box.padding = 0.8, nudge_y = 0.001) +
  scale_color_manual(values = c('#993366',"#333366","#663366")) + 
  ylab("Predicted using my method, averaged across samples") +
  xlab("Dijk's CNH") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "top",
        legend.text = element_text(size=24, colour='black'),
        legend.title = element_blank() ) 

jpeg('tempfig.jpeg', width = 1000, height = 1000, pointsize = 12)
plot
dev.off()

