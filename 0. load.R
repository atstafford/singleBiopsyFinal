# Output from section 1 (copy number and data prep files) ======================
car.raw <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.raw.rds")
car.info <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.info.rds")
car.clonality <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.clonality.rds")
car.diversity <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.diversity.rds")
car.actualITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.actualITH.rds")

ad.raw <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.raw.rds")
ad.info <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.info.rds")
ad.clonality <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.clonality.rds")
ad.diversity <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.diversity.rds")
ad.actualITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.actualITH.rds")

validationCar.raw <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.raw.rds")
validationCar.info <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.info.rds")
validationCar.clonality <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.clonality.rds")
validationCar.diversity <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.diversity.rds")
validationCar.actualITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.actualITH.rds")

validationAd.raw <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.raw.rds")
validationAd.info <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.info.rds")
validationAd.clonality <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.clonality.rds")
validationAd.diversity <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.diversity.rds")
validationAd.actualITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationAd.actualITH.rds")

tracerx.raw <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.raw.rds")
tracerx.info <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.info.rds")
tracerx.clonality <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.clonality.rds")
tracerx.diversity <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.diversity.rds")
tracerx.actualITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.actualITH.rds")

coad.raw <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/coad.raw.rds")
coad.info <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/coad.info.rds")

ad.matrices <- ("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.matrices.rds")
car.matrices <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.matrices.rds")
tracerx.matrices <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.matrices.rds")

# Output from section 6 (heatmep correlation matrix) ===========================
cna.corr.car <- readRDS('~/Documents/CNA/Github/singleBiopsyITH/Data/cna.corr.car.rds')
cna.corr.tracerx <- readRDS('~/Documents/CNA/Github/singleBiopsyITH/Data/cna.corr.tracerx.rds')

# Output from section 7 (CNA clustering) =======================================
sig.hclust <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/sig.hclust.rds")

# Output from section 8 (single bin models) ====================================
uniReg.out.list <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/uniReg.out.list.rds")

# Output from section 9 (candidate and representative bins) ====================
candidate.bins <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/candidate.bins.rds")
representitive.bins <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins.rds")

# Output from section 10 (multibin models) =====================================
betareg_52 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betareg_51.rds")
betareg_40 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betareg_40.rds")
betareg_39 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betareg_39.rds")

# Output from section 11 (hg19predictors and predicted ITH) ====================
hg19predictors <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/hg19predictors.rds")
car.predictedITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/car.predictedITH.rds")
validationCar.predictedITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/validationCar.predictedITH.rds")
tracerx.predictedITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/tracerx.predictedITH.rds")
ad.predictedITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/ad.predictedITH.rds")
coad.predictedITH <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/coad.predictedITH.rds")

# Output from section 12 (alternative representative bins) =====================
mods_secondrep <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/mods_secondrep.rds")
bins_secondrep <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/bins_secondrep.rds")
clusters_secondrep <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/clusters_secondrep.rds")
representitive.bins_secondrep <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins_secondrep.rds")

mods_nonmax <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/mods_nonmax.rds")
bins_nonmax <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/bins_nonmax.rds")
clusters_nonmax <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/clusters_nonmax.rds")
representitive.bins_nonmax <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins_nonmax.rds")

mods_rand <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/mods_rand.rds")
bins_rand <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/bins_rand.rds")
clusters_rand <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/clusters_rand.rds")
representitive.bins_rand <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/representitive.bins_rand.rds")

# Output from section 16 (genes in 616 candidate bins and bootstrapping) =======
betaregCandGenes <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betaregCandGenes.rds")
candGOboot1 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/candGOboot1.rds")
candGOboot2 <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/candGOboot2.rds")

# Output from section 17 (genes represented in model) ==========================
betaregGenes <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betaregGenes.rds")
betaregRepGenes <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/betaregRepGenes.rds")
GOI.acrossclusters <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/GOI.acrossclusters.rds")

# Output from section 19 (differentially expressed genes) ======================
countData <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/countData.rds")
metaData <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/metaData.rds")
dds <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/dds.rds")
dde <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/dde.rds")
res <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/res.rds")
res.sig <- readRDS("~/Documents/CNA/Github/singleBiopsyITH/Data/res.sig.rds")




