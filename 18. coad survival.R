# PREREQUISITS: load section 11 data 

COAD_survival <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/COAD_clinical.rds")

# Add predicted ITH back into COAD_clinical
COAD_survival$ITH <- coad.predictedITH[match(COAD_survival$TCGA_barcode, coad.predictedITH$TCGA_barcode), 3]

# Keep only those samples whose ITH could be predicted
COAD_survival <- COAD_survival[!is.na(COAD_survival$ITH),]

# Mark those patients that didnt have a CNA in any of the predictive bins
COAD_survival$intercept <- "no"
COAD_survival$intercept[which(COAD_survival$ITH==coef(betareg_39)["(Intercept)"])] <- "yes"

# Set up KM 
COAD_survival$OS_cens <- 0
COAD_survival$OS_cens[which(COAD_survival$vital_status=="Dead")] <- 1
COAD_survival <- COAD_survival[-which(COAD_survival$days_to_death==1 | COAD_survival$days_to_death==0),]

COAD_survival$days_to_last_follow_up <- as.numeric(as.character(COAD_survival$days_to_last_follow_up))
COAD_survival$days_to_death <- as.numeric(as.character(COAD_survival$days_to_death))

COAD_survival$OS_time <- NA

for (i in 1:nrow(COAD_survival)) {
  if (COAD_survival$OS_cens[i] == 0) {
    COAD_survival$OS_time[i] <- COAD_survival$days_to_last_follow_up[i]/365
  }
  else if (COAD_survival$OS_cens[i] == 1) {
    COAD_survival$OS_time[i] <- COAD_survival$days_to_death[i]/365
  }
}

COAD_survival$stage2 <- NA
COAD_survival$stage2[which(COAD_survival$stage == "Stage I" | COAD_survival$stage == "Stage II")] <- "early"
COAD_survival$stage2[which(COAD_survival$stage == "Stage III" | COAD_survival$stage == "Stage IV")] <- "late"

COAD_survival$divGroup4 <- ntile(COAD_survival$ITH, 4)
COAD_survival$divGroup3 <- "med"
COAD_survival$divGroup3[which(COAD_survival$divGroup4==4)] <- "top25"
COAD_survival$divGroup3[which(COAD_survival$divGroup4==1)] <- "bottom25"

COAD_survival$MSI <- substr(COAD_survival$MSI, 1, 3)

# hist of MSI MSS
dist <- ggplot(COAD_survival, aes(x=ITH)) + 
  geom_histogram(data=subset(COAD_survival, MSI == 'MSS'), fill = "#1B9E77", alpha = 0.3, binwidth = 0.02) +
  geom_histogram(data=subset(COAD_survival, MSI == 'MSI'), fill = "#E6AB02", alpha = 0.3, binwidth = 0.02) +
  geom_density(data=subset(COAD_survival, MSI == 'MSS'), aes(fill=MSI), color = "black", alpha = 0.5) +
  geom_density(data=subset(COAD_survival, MSI == 'MSI'), color = 'black', aes(fill = MSI), alpha = 0.5) +
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
hist <- ggplot(COAD_survival, aes(round(ITH,2))) + 
  geom_histogram(binwidth = 0.01, color="black", aes(fill = divGroup3)) +
  xlab("Predicted CNA diversity") +
  scale_fill_manual(name=NULL, labels = c("bottom25%", "mid", "top25%"), values = c("bottom25"="#000000", "med"="#CCCCCC", "top25"="#CC0033")) +
  theme_minimal() +
  theme(legend.title=element_text(size=24),
        legend.text=element_text(size=24),
        legend.position = "top",
        axis.title = element_text(size=24),
        axis.text = element_text(size=24))

# pan stage MSS
z <- COAD_survival[which(COAD_survival$MSI=="MSS"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSS Pan-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
panOS.MSS <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
panOS.MSS

# pan stage MSI
z <- COAD_survival[which(COAD_survival$MSI=="MSI"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSI Pan-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
panOS.MSI <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
panOS.MSI

# early MSS
z <- COAD_survival[which(COAD_survival$MSI=="MSS" & COAD_survival$stage2=="early"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSS early-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
earlyOS.MSS <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=7, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
#jpeg('tempfig.jpeg', width = 400, height = 375)
#earlyOS.MSS
#dev.off()

# late MSS
z <- COAD_survival[which(COAD_survival$MSI=="MSS" & COAD_survival$stage2=="late"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSS late-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
lateOS.MSS <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
lateOS.MSS

# early MSI
z <- COAD_survival[which(COAD_survival$MSI=="MSI"),]
z <- z[which(z$stage2=="early"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSI early-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
earlyOS.MSI <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
earlyOS.MSI

# late MSI
z <- COAD_survival[which(COAD_survival$MSI=="MSI"),]
z <- z[which(z$stage2=="late"),]
z <- z[which(z$divGroup3 != "med"),]

title <- paste("MSI late-stage (n=", nrow(z), ")", "\nhigh=", nrow(z[which(z$divGroup3=="top25"),]), " low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = "")
fit <- survival::survfit(Surv(OS_time, OS_cens) ~ z$divGroup3, data = z)
lateOS.MSI <- survminer::ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#000000","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")
lateOS.MSI

plot <- plot_grid(dist, hist, panOS.MSS$plot, panOS.MSI$plot, earlyOS.MSS$plot, earlyOS.MSI$plot, lateOS.MSS$plot, lateOS.MSI$plot, ncol = 2)

jpeg('tempfig.jpeg', width = 888, height = 1800)
plot
dev.off()


# ITH per stage 
# Perform chi squared test
chisq <- chisq.test(table(COAD_survival$divGroup3, COAD_survival$stage))


