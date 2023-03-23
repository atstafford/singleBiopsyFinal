# PREREQUISITS: load section 11 data 

car_survival <- read_excel("~/Documents/CNA/Github/singleBiopsyITH/Data/clinical_data_training_cohort.xlsx", col_names = TRUE, skip =0)

# add in ITH
car_survival$actualITH <- car.predictedITH[match(car_survival$sampleID, car.predictedITH$patient), 3]
car_survival$predictedITH <- car.predictedITH[match(car_survival$sampleID, car.predictedITH$patient), 4]
car_survival <- car_survival[!is.na(car_survival$predictedITH),]

# split
car_survival$divGroup4 <- ntile(car_survival$predictedITH, 4)
car_survival$divGroup3 <- "med"
car_survival$divGroup3[which(car_survival$divGroup4==4)] <- "top25"
car_survival$divGroup3[which(car_survival$divGroup4==1)] <- "bottom25"

# pan stage
z <- car_survival
z <- z[which(z$divGroup3 != "med"),]

title <- paste("Pan-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(RFS_time, RFS_cens) ~ z$divGroup3, data = z)
panOS.MSS <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                        font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                        ylab="Survival Probability", xlab="Time (years)")
panOS.MSS