READ_survival <- READ_clinical

# Add predicted ITH back into COAD_clinical
READ_survival$ITH <- tcga.predictedITH[match(READ_survival$TCGA_barcode, read.predictedITH$TCGA_barcode), 2]

# Keep only those samples whose ITH could be predicted
READ_survival <- READ_survival[!is.na(READ_survival$ITH),]

# Mark those patients that didnt have a CNA in any of the predictive bins
READ_survival$intercept <- "no"
READ_survival$intercept[which(READ_survival$ITH==coef(betareg_33bin)["(Intercept)"])] <- "yes"

# Set up KM 
READ_survival$OS_cens <- 0
READ_survival$OS_cens[which(READ_survival$vital_status=="Dead")] <- 1
READ_survival <- READ_survival[-which(READ_survival$days_to_death==1 | READ_survival$days_to_death==0),]

READ_survival$days_to_last_follow_up <- as.numeric(as.character(READ_survival$days_to_last_follow_up))
READ_survival$days_to_death <- as.numeric(as.character(READ_survival$days_to_death))

READ_survival$OS_time <- NA

for (i in 1:nrow(READ_survival)) {
  if (READ_survival$OS_cens[i] == 0) {
    READ_survival$OS_time[i] <- READ_survival$days_to_last_follow_up[i]/365
  }
  else if (READ_survival$OS_cens[i] == 1) {
    READ_survival$OS_time[i] <- READ_survival$days_to_death[i]/365
  }
}

READ_survival$stage2 <- NA
READ_survival$stage2[which(READ_survival$stage == "Stage I" | READ_survival$stage == "Stage II")] <- "early"
READ_survival$stage2[which(READ_survival$stage == "Stage III" | READ_survival$stage == "Stage IV")] <- "late"

READ_survival$divGroup4 <- ntile(READ_survival$ITH, 4)
READ_survival$divGroup3 <- "med"
READ_survival$divGroup3[which(READ_survival$divGroup4==4)] <- "top25"
READ_survival$divGroup3[which(READ_survival$divGroup4==1)] <- "bottom25"

READ_survival$MSI <- substr(READ_survival$MSI, 1, 3)

# hist of MSI MSS
dist <- ggplot(READ_survival, aes(x=ITH)) + 
  geom_histogram(data=subset(READ_survival, MSI == 'MSS'), fill = "#1B9E77", alpha = 0.3, binwidth = 0.02) +
  geom_histogram(data=subset(READ_survival, MSI == 'MSI'), fill = "#E6AB02", alpha = 0.3, binwidth = 0.02) +
  geom_density(data=subset(READ_survival, MSI == 'MSS'), aes(fill=MSI), color = "black", alpha = 0.5) +
  geom_density(data=subset(READ_survival, MSI == 'MSI'), color = 'black', aes(fill = MSI), alpha = 0.5) +
  scale_fill_manual(name=NULL, breaks=c('MSS','MSI'),
                     values=c('MSS'='#1B9E77', 'MSI'='#E6AB02')) +
  theme_minimal() +
  theme(legend.title=element_text(size=24),
      legend.text=element_text(size=24),
      legend.position = "top",
      axis.title = element_text(size=24),
      axis.text = element_text(size=24))

ks.test(COAD_survival$ITH[which(COAD_survival$MSI=="MSS")], COAD_survival$ITH[which(COAD_survival$MSI=="MSI")])

# hist of cut
hist <- ggplot(READ_survival, aes(round(ITH,2))) + 
  geom_histogram(binwidth = 0.01, color="black", aes(fill = divGroup3)) +
  xlab("Predicted CNA diversity") +
  scale_fill_manual(name=NULL, labels = c("bottom25%", "mid", "top25%"), values = c("bottom25"="#7570B3", "med"="#CCCCCC", "top25"="#D95F02")) +
  theme_minimal() +
  theme(legend.title=element_text(size=24),
        legend.text=element_text(size=24),
        legend.position = "top",
        axis.title = element_text(size=24),
        axis.text = element_text(size=24))

# pan stage MSS
z <- READ_survival
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSS Pan-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
panOS <- ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#7570B3","#D95F02"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
panOS

# early 
z <- READ_survival[which(READ_survival$stage2=="early"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("Early-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
earlyOS <- ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=7, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
earlyOS

# late 
z <- READ_survival[which(READ_survival$stage2=="late"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSS late-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
lateOS <- ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#7570B3","#D95F02"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
lateOS



