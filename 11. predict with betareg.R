# PREREQUISITS: load section 1,9,10 data 

# choose model
beta <- betareg_39
summary(beta, type = "pearson")
round(AIC(beta),2)

# pull predictors and locatinos for hg19
predictors <- rownames(data.frame(summary(beta)$coefficients))[-1]
hg19predictors <- data.frame(bin=as.numeric(str_extract_all(predictors, "[0-9]+")), cna=str_sub(predictors,-4,-1))
hg19predictors <- merge(unique(hg19predictors), car.info$start.stop)
hg19predictors$cluster <- candidate.bins[match(hg19predictors$bin, candidate.bins$bin), 11] 
hg19predictors$pcSubclonal <- NA
both <- hg19predictors$bin[which(duplicated(hg19predictors$bin)==T)]
hg19predictors$cna[which(hg19predictors$bin %in% hg19predictors$bin[which(duplicated(hg19predictors$bin)==T)] )] <- 'both'

hg19predictors <- hg19predictors[c(1,3:7,2)]
hg19predictors <- hg19predictors[!duplicated(hg19predictors),]

# predict training
car.predict <- predictITHfromBeta(rawData = car.raw, dataInfo = car.info, predictors = hg19predictors, title="Training", 
                        actualITH = car.actualITH, beta = beta)
car.predictedITH <- car.predict$predictedITH
car.predict$plot
summary(lm(actual ~ predicted, data=car.predictedITH))
betaPredictPlot1 <- annotate_figure(plot_grid(car.predict$plot, ncol=1), 
                                    left = text_grob('CNAdf - predicted', rot = 90, size = 28, vjust = 1),
                                    bottom = text_grob('CNAdf - actual', size = 28, vjust = 0))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
betaPredictPlot1
dev.off()

# predict validation cohort
validationCar.predict <- predictITHfromBeta(rawData = validationCar.raw, dataInfo = validationCar.info, predictors = hg19predictors,title="Validation",  
                                  actualITH = validationCar.actualITH, beta = beta)
validationCar.predictedITH <- validationCar.predict$predictedITH
validationCar.predict$plot
summary(lm(actual ~ predicted, data=validationCar.predictedITH))
betaPredictPlot1 <- annotate_figure(plot_grid(validationCar.predict$plot, ncol=1), 
                                    left = text_grob('CNAdf - predicted', rot = 90, size = 28, vjust = 1),
                                    bottom = text_grob('CNAdf - actual', size = 28, vjust = 0))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
betaPredictPlot1
dev.off()

# plot training and validation together  
betaPredictPlot1 <- annotate_figure(plot_grid(car.predict$plot, validationCar.predict$plot, ncol=2), 
                                     left = text_grob('CNA diversity', rot = 90, size = 28, vjust = 1),
                                     bottom = text_grob('Predicted CNA diversity', size = 28, vjust = 0))
  
jpeg('tempfig.jpeg', height = (3*37.795*4.30),  width = (3*37.795*8.38))
betaPredictPlot1
dev.off()

# predict tracerx cohort
tracerx.predict <- predictITHfromBeta(rawData = tracerx.raw, dataInfo = tracerx.info, predictors = hg19predictors,title="NSCLC",  
                                            actualITH = tracerx.actualITH, beta = beta)
tracerx.predictedITH <- tracerx.predict$predictedITH
tracerx.predict$plot

betaPredictPlot2 <- annotate_figure(plot_grid(tracerx.predict$plot, ncol=1), 
                                    left = text_grob('CNAdf - predicted', rot = 90, size = 28, vjust = 1),
                                    bottom = text_grob('CNAdf - actual', size = 28, vjust = 0))

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
betaPredictPlot2
dev.off()

# predict tracerx cohort no staqge 1
stage <- read_excel("~/Documents/CNA/Github//Data/TracerX/cn_data.xlsx", sheet = 1, skip = 1)
stage <- stage[,c(1,2)]
stage$Stage <- substr(stage$Stage, 1, 1)
stage23NSCLC <- stage$TRACERxID[which(stage$Stage!=1)]
stage23NSCLC.raw <- tracerx.raw[,which(substr(colnames(tracerx.raw), 1, nchar(colnames(tracerx.raw))-3) %in% stage23NSCLC)]
stage23NSCLC.raw <- cbind(tracerx.raw[,c(1:3)], stage23NSCLC.raw)
stage23NSCLC.info <- PullDataInfo(stage23NSCLC.raw)
stage23NSCLC.diversity <- PullDataDiversity(rawdata = stage23NSCLC.raw, dataInfo = stage23NSCLC.info)
stage23NSCLC.actualITH <- data.frame(patient = lapply(data.frame(patient=stage23NSCLC.info$patientIDs), rep, stage23NSCLC.info$sampPerPatient),
                                sample = stage23NSCLC.info$sampleIDs,
                                actual = lapply(data.frame(actual=stage23NSCLC.diversity$pic.frac$pic.frac), rep, stage23NSCLC.info$sampPerPatient))


stage23NSCLC.predict <- predictITHfromBeta(rawData = stage23NSCLC.raw, dataInfo = stage23NSCLC.info, predictors = hg19predictors,title="NSCLC (stage 2/3)",  
                                      actualITH = stage23NSCLC.actualITH, beta = beta)
stage23NSCLC.predictedITH <- stage23NSCLC.predict$predictedITH
stage23NSCLC.predict$plot

betaPredictPlot2 <- annotate_figure(plot_grid(stage23NSCLC.predict$plot, ncol=1), 
                                    left = text_grob('CNAdf - predicted', rot = 90, size = 28, vjust = 1),
                                    bottom = text_grob('CNAdf - actual', size = 28, vjust = 0))

#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
betaPredictPlot2
dev.off()

# predict adenomas cohort
ad.predict <- predictITHfromBeta(rawData = ad.raw, dataInfo = ad.info, predictors = hg19predictors,title=NULL,  
                                      actualITH = ad.actualITH, beta = beta)
ad.predictedITH <- ad.predict$predictedITH
ad.predict$plot
summary(lm(actual ~ predicted, data=ad.predictedITH))

betaPredictPlot3 <- annotate_figure(plot_grid(ad.predict$plot, ncol=1), 
                                    left = text_grob('CNAdf - predicted', rot = 90, size = 28, vjust = 1),
                                    bottom = text_grob('CNAdf - actual', size = 28, vjust = 0))
#jpeg('tempfig.jpeg', width = (3*37.795*4.8), height = (3*37.795*4.19))
jpeg('tempfig.jpeg', width = (20), height = (20), units = 'cm', res = 300)
betaPredictPlot3
dev.off()

# predict COAD
coad.predict <- predictITHfromBeta(rawData = coad.raw, dataInfo = coad.info, predictors = hg19predictors, title=NULL,  
                                 actualITH = NA, beta = beta)
coad.predictedITH <- coad.predict$predictedITH
coad.predictedITH <- cbind(TCGA_barcode=substr(coad.info$patientIDs,1,12), coad.predictedITH)

# predict READ
read.predict <- predictITHfromBeta(rawData = read.raw, dataInfo = read.info, predictors = hg19predictors, title=NULL,  
                                   actualITH = NA, beta = beta)
read.predictedITH <- read.predict$predictedITH
read.predictedITH <- cbind(TCGA_barcode=substr(read.info$patientIDs,1,12), read.predictedITH)


# freq of predictors in ad and car
x <- ad.raw
x <- car.raw
rownames(x) <- 1:nrow(x)
x <- x[hg19predictors$bin,-c(1:3)]
colnames(x)[which(colnames(x)=="A12.1.1")] <- "A12.2"
x$MCR <- rownames(x)
x <- melt(x, id='MCR')
colnames(x) <- c('MCR',"Sample","CNA")
x$patient <- sub("\\..*", "",x$Sample)
x$CNA <- ifelse(x$CNA==1, "loss", ifelse(x$CNA==2,"diploid","gain"))
x$chr <- hg19predictors[match(x$MCR, hg19predictors$bin), 2]
y <- ggplot(x, aes(x=Sample, y=MCR, fill = as.factor(CNA))) +
  geom_tile(color="black") +
  scale_fill_manual(name = "CNA",values = c(alpha("#CCCCCC",0.2),"#B40F20","#046C9A")) +
  facet_grid(cols=vars(patient), rows=vars(chr), scales = 'free', space = 'free') +
  theme_custom() +
  theme(
    legend.position = "top",
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_blank(),
    panel.border=element_blank(),
    #axis.text.x = element_text(angle=90, vjust=0.5, hjust=),
    axis.text.x = element_blank(),
    #strip.text.x = element_blank(),
    strip.text.x = element_text(size=28, angle=90),
    strip.text.y = element_text(size=28),
    strip.background.y = element_rect(colour = 'white', fill=NULL),
    strip.background.x = element_rect(colour = 'white', fill=NULL),
    legend.text = element_text(size=28),
    legend.title = element_text(size=28)) 
#jpeg('tempfig.jpeg', width = (35), height = (35), units = 'cm', res = 300)
jpeg('tempfig.jpeg', width = (65), height = (35), units = 'cm', res = 300)
y
dev.off()

# save
saveRDS(hg19predictors, "~/Documents/CNA/Github/singleBiopsyITH/Data/hg19predictors.rds")
saveRDS(car.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.predictedITH.rds")
saveRDS(validationCar.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.predictedITH.rds")
saveRDS(tracerx.predictedITH,"~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.predictedITH.rds")
saveRDS(ad.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.predictedITH.rds")
saveRDS(coad.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/coad.predictedITH.rds")
























# predict EPICC
multiReg.in <- t(epicc.raw[,-c(1:3)]) 
multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(multiReg.in) <- epicc.info$sampleIDs
multiReg.in <- multiReg.in[,c(hg19predictors$bin)]
names(multiReg.in) = paste("bin_", names(multiReg.in), sep="")
multiReg.in = multiReg.in %>%
  mutate(across(everything(), as.character))
multiReg.in <- dummy_cols(multiReg.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
rownames(multiReg.in) <- epicc.info$sampleIDs
multiReg.in$bin_407_gain <- 0
multiReg.in$bin_606_gain <- 0
multiReg.in$bin_898_gain <- 0
multiReg.in$bin_1074_gain <- 0
multiReg.in$bin_1139_gain <- 0
multiReg.in$bin_1158_gain <- 0
multiReg.in$bin_1173_gain <- 0
multiReg.in$bin_1193_gain <- 0
multiReg.in$bin_1698_loss <- 0
multiReg.in$bin_1816_gain <- 0
multiReg.in$bin_2306_gain <- 0
multiReg.in$bin_2326_gain <- 0
multiReg.in$bin_2336_gain <- 0
multiReg.in$bin_2385_gain <- 0
multiReg.in$bin_1455_gain <- 0
multiReg.in$bin_1470_gain <- 0
multiReg.in$bin_1488_gain <- 0
multiReg.in$bin_1524_loss <- 0
multiReg.in$bin_2083_gain <- 0
multiReg.in$bin_2371_gain <- 0
epicc.predictedITH <- data.frame(predicted = predict(beta, multiReg.in))

epicc.predictedITH <- data.frame(cbind(epicc.actualITH, epicc.predictedITH))
epicc.predictedITH$type <- 'epicc'

# Save
saveRDS(hg19predictors, "~/Documents/CNA/Github/singleBiopsyITH/Data/hg19predictors.rds")
saveRDS(car.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/car.predictedITH.rds")
saveRDS(validationCar.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.predictedITH.rds")
saveRDS(tracerx.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.predictedITH.rds")
saveRDS(tcga.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/tcga.predictedITH.rds")
saveRDS(ad.predictedITH, "~/Documents/CNA/Github/singleBiopsyITH/Data/ad.predictedITH.rds")






# old
# convert hg19 predictors into hg38 for tracerx predictions
iddf <- hg19predictors[,c(2:4)]
iddf$chr <- as.character(paste('chr',iddf$chr,sep = "", collapse = NULL))
idrange <- makeGRangesFromDataFrame(iddf)
chain <- import.chain("~/Documents/CNA/Data/hg19Tohg38.over.chain")
hg38ids <- liftOver(idrange, chain)
numfound <- unlist(lapply(hg38ids,length))
hg38ids <- as.data.frame(hg38ids)
hg38ids$group_name <- rep(hg19predictors$bin, numfound)

hg38predictors <- data.frame(bin=hg19predictors$bin, chr=NA, start=NA, stop=NA)

for ( i in 1:nrow(hg19predictors) ) {
  wd <- hg38ids[which(hg38ids$group_name==hg19predictors$bin[i]),]
  wd <- wd[which(wd$seqnames == paste("chr",hg19predictors$chr[i],sep="")), ]
  
  hg38predictors$chr[i] <- sub('...','',wd$seqnames[1])
  hg38predictors$start[i] <- min(wd$start)
  hg38predictors$stop[i] <- max(wd$end)
}

# pull copy number for hg38 predictor bins
x <- newAlignBins(bins = hg38predictors, cn.list = list(tracerx.raw))
x <- x[[1]]
x[ ,-c(1:4)] <- ifelse(x[ ,-c(1:4)] >= 3, 3, ifelse(x[ ,-c(1:4)] < 2, 1, 2))
rownames(x) <- x$bin

saveRDS(hg38predictors, "~/Documents/CNA/Github/singleBiopsyITH/Data/hg38predictors.rds")