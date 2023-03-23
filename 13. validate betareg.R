# PREREQUISITS: load section 1,10 data 
beta <- betareg_39

# bootstrap resampling of random 30 bins to generate models to compare AIC
rand30 <- t(car.raw[,-c(1:3)])
rand30 <- data.frame(apply(rand30, 2, as.character), check.names = FALSE)
rand30 <- data.frame(lapply(rand30, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(rand30) <- car.info$sampleIDs

rand30 <- cbind(ITH=car.actualITH$actual, rand30)
names(rand30) = paste("bin_", names(rand30), sep="")
rand30 = rand30 %>%
  mutate(across(everything(), as.character))

rand30 <- dummy_cols(rand30[-1], remove_selected_columns = TRUE, remove_first_dummy = TRUE) 
rand30$bin_ITH <- car.actualITH$actual

set.seed(12345)
pR2 <- list()
ll <- list()
tb <- list()
aic2 <- list()
for (i in 1:100000) {
  print(i)
  testbins <- sample(c(1:car.info$noBins), 30)
  tb[[i]] <- testbins
  input <- rand30[ , which(as.numeric(gsub("[^\\d]+", "", colnames(rand30), perl=TRUE)) %in% testbins)]
  input$bin_ITH <- rand30$bin_ITH
  try({
    beta <- betareg(bin_ITH ~., data = input)
    pR2[[i]] <- round(beta$pseudo.r.squared,3)
    ll[[i]] <- round(beta$loglik,2)
    aic2[[i]] <- round(AIC(beta),2)
  }, silent = TRUE)
}

length(unlist(aic2))
z <- data.frame(unlist(aic2))
mean(z$unlist.aic2.)
min(z$unlist.aic2.)

plot <- ggplot(z, aes(x=unlist.aic2.)) + 
  geom_histogram(binwidth=2, fill="#00A08A", color="black") +
  geom_vline(xintercept = AIC(betareg_39), linetype="longdash", color="blue") +
  xlim(c(-380,-250)) +
  xlab("AIC") +
  ylab("Frequency") +
  theme(plot.margin=margin(t=0,r=0.5,b=0,l=0,"cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        axis.ticks.length=unit(0.2, "cm"),
        legend.position = "none")

jpeg('tempfig.jpeg', width = 400, height = 375)
plot
dev.off()

# LOOCV cannot be performed in normal way as too few betareg models can be built
train.control <- caret::trainControl(method = "LOOCV")
model <- caret::train(actual ~ predicted, data = car.predictedITH, method = "lm",
               trControl = train.control)
print(model)


