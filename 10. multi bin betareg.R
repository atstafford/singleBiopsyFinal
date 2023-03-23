# PREREQUISITS: load section 1,9 data 

# Build input for 51-bin model
multiReg.in <- t(car.raw[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(multiReg.in) <- car.info$sampleIDs
multiReg.in <- multiReg.in[,c(representitive.bins$bin)] # Keep only the candidate bins (columns)

names(multiReg.in) = paste("bin_", names(multiReg.in), sep="")
multiReg.in = multiReg.in %>%
  mutate(across(everything(), as.character))

# Make dummy for beta
multiReg.in <- dummy_cols(multiReg.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
multiReg.in$bin_ITH <- car.actualITH$actual
remove <- attributes(alias(lm(bin_ITH ~ ., data = multiReg.in))$Complete)$dimnames[[1]]
multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]

# 51-bin betareg
betareg_51 <- betareg:betareg(bin_ITH ~., data = multiReg.in)
summary(betareg_51, type = "pearson")
round(AIC(betareg_51),2)
vif <- data.frame(car::vif(betareg_51))

# Backwards AIC selection on 51 bin model
variables <- multiReg.in[ ,-which(colnames(multiReg.in) %in% c("bin_ITH"))]
selection_back <- frmselection::betaselect(variables, car.actualITH$actual, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
keep <- selection_back$variable
multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %in% keep]
multiReg.in$bin_ITH <- car.actualITH$actual
betareg_40 <- betareg::betareg(bin_ITH ~., data = multiReg.in)
summary(betareg_40, type = "pearson")
round(AIC(betareg_40),2)

# sequential removal of VIF>10 variables
new.model <- betareg_40
vif <- data.frame(vif(new.model))
while( max(vif$vif.new.model.) >=10 ) {
  remove <- rownames(vif)[which(vif$vif.new.model.== max(vif$vif.new.model.))]
  multiReg.in <- multiReg.in[ ,colnames(multiReg.in) %!in% remove]
  new.model <- betareg::betareg(bin_ITH ~., data = multiReg.in)
  vif <- data.frame(car::vif(new.model))
}

betareg_39 <- new.model
summary(betareg_39, type = "pearson")
round(AIC(betareg_39),2)
#lrtest(betareg_33bin)

# Save
saveRDS(betareg_51, "~/Documents/CNA/Github/singleBiopsyITH/Data/betareg_51.rds")
saveRDS(betareg_40, "~/Documents/CNA/Github/singleBiopsyITH/Data/betareg_40.rds")
saveRDS(betareg_39, "~/Documents/CNA/Github/singleBiopsyITH/Data/betareg_39.rds")


