#==============================================================================
# Supplementary R scripts for mater's thesis
# Title: Potential of FTIR spectroscopy for tropical dendrochronology
# Author: Lukas Erzfeld
# Insitution: University of Leipzig
# Date: August 26, 2025
#===============================================================================

# required packages may have to be installed primarily
# install.packages(c("prospectr", "mdatools", "readxl", "writexl"))
library(prospectr)
library(mdatools)
library(readxl) 
library(writexl)
source("Utils.R")
#-------------------------------------------------------------------------------
# NOTE: number of latent variables used in the model has to be chosen in line 240
#-------------------------------------------------------------------------------

# out path of 00_PLS-DA_data_prep_master.R
path = "data/FTIR_Spectra/PLS-DA/"
species = "Pinus_patula"
filename_meta <- paste0(species, "_file_meta_PLS-DA.csv")
filename_dat <- paste0(species, "_data_mat_PLS-DA.csv")

path_meta <- file.path(path, filename_meta)
path_spectra <- file.path(path, filename_dat)
  
# read saved data
file_meta <- read.csv(path_meta)
data <- read.csv(path_spectra)
colnames(data) <- gsub("^X", "", colnames(data))
colnames(data) <- as.numeric(colnames(data))
wn <- as.numeric(colnames(data))
#-------------------------------------------------------------------------------

# 1 spectrum of P. patula seems to be outlier --> rm
ref_band <- "1025.9445"
if (species == "Pinus_patula") {
  n_rm <- 1
  rm_ix_high <- tail(order(data[["1620.8765"]]), n_rm)
  data <- data[-rm_ix_high, ]
  file_meta <- file_meta[-rm_ix_high, ]
  rm_ix <- order(data[[ref_band]])[1:n_rm]
  data <- data[-rm_ix, ]
  file_meta <- file_meta[-rm_ix, ]
}
#-------------------------------------------------------------------------------

# slice data to fingerprint region 1800 - 800 cm-1
fingerprint_cols <- wn >= 800 & wn <= 1800
wn4000 <- wn
wn <- wn[fingerprint_cols]
spec_fingerprint <- data[, fingerprint_cols]

### preprocessing ###
mXr00 <- as.matrix(spec_fingerprint)
mXr01 <- sweep(mXr00, MARGIN=1, apply(mXr00, 1, function(x) (mean(x))), FUN="/") # mean scaling
mXr02 <- scale(mXr00, center = TRUE, scale = TRUE) # z-transform
mXr03 <- standardNormalVariate(mXr00)
mXr04 <- prospectr::detrend(mXr00, wn)

# SMOOTHING + DERIVATIVES + SNV
mXrSmoothed <- savitzkyGolay(mXr00, m=0, p=2, w=15)
mXrSG1 <- savitzkyGolay(mXrSmoothed, m=1, p=2, w=3)
mXrSG2 <- savitzkyGolay(mXrSmoothed, m=2, p=2, w=3)
mXr05 <- standardNormalVariate(mXrSG1)
mXr06 <- standardNormalVariate(mXrSG2)

# DERIVATIVES + SNV
mXrD1 <- savitzkyGolay(mXr00, m=1, p=1, w=3)
mXrD2 <- savitzkyGolay(mXr00, m=2, p=2, w=3)
if (species == "Pinus_patula") {
  mXrD1 <- savitzkyGolay(mXrD1, m=0, p=2, w=15) # smoothing
  mXrD2 <- savitzkyGolay(mXrD2, m=0, p=2, w=15) # smoothing
}
mXr07 <- standardNormalVariate(mXrD1)
mXr08 <- standardNormalVariate(mXrD2)

# DERIVATIVES + VECTOR NORMALIZATION (PAPER)
mXr09 <- t(apply(mXrD1, 1, vector_normalize))
mXr10 <- t(apply(mXrD2, 1, vector_normalize))

# BASELINE + SNV
mXrBase <- baseline(data)
colnames(mXrBase) <- wn4000
#mXrBaseBand <- manual_baseline_correction(data, wn=as.numeric(wn4000), anchors = c(3750, 1800, 1540, 785))
mXrBaseBand <- manual_baseline_correction(data, wn=as.numeric(wn4000), anchors = c(3750, 1780, 1485, 1185, 820), window = 50)
mXrBase <- mXrBase[, fingerprint_cols]
mXrBaseBand <- mXrBaseBand[, fingerprint_cols]

mXr11 <- standardNormalVariate(mXrBase)
mXr12 <- standardNormalVariate(mXrBaseBand)

# BASELINE + BAND NORMALIZATION 1424
mXr13 <- normalize_to_band(mXrBase, wn, band = 1424) # mal lieber mit 1026 bzw höchstem peak probieren?
mXr13 <- mXr13[, apply(mXr13, 2, sd) > 0] # remove this band, otherwise PLS-DA returns NaN
mXr14 <- normalize_to_band(mXrBaseBand, wn, band = 1424)
mXr14 <- mXr14[, apply(mXr14, 2, sd) > 0]

# BASELINE + BAND NORMALIZATION 1026
mXr17 <- normalize_to_band(mXrBase, wn, band = 1026)
mXr17 <- mXr17[, apply(mXr17, 2, sd) > 0] # remove this band, otherwise PLS-DA returns NaN
mXr18 <- normalize_to_band(mXrBaseBand, wn, band = 1026)
mXr18 <- mXr18[, apply(mXr18, 2, sd) > 0] # remove this band, otherwise PLS-DA returns NaN

# 2nd DERIVATIVE + SNV (poly 3)
mXrD2_p3 <- savitzkyGolay(mXr00, m=2, p=3, w=5)
if (species == "Pinus_patula") {
  mXrD2_p3 <- savitzkyGolay(mXrD2_p3, m=0, p=2, w=15) # smoothing
}
mXr19 <- standardNormalVariate(mXrD2)

# BASEBAND + DERIVATIVE
mXr20 <- savitzkyGolay(mXrBaseBand, m=1, p=1, w=3)
mXr21 <- savitzkyGolay(mXrBaseBand, m=2, p=2, w=3)
if (species == "Pinus_patula") {
  mXr20 <- savitzkyGolay(mXr20, m=0, p=2, w=15) # smoothing
  mXr21 <- savitzkyGolay(mXr21, m=0, p=2, w=15) # smoothing
}

# BASEBAND + DERIVATIVE + SNV
mXr22 <- standardNormalVariate(mXr20)
mXr23 <- standardNormalVariate(mXr21)

# MSC (+ DERIVATIVES)
mXr24 <- msc(mXr00)
mXr25 <- savitzkyGolay(mXr24, m=1, p=1, w=3) # paper
mXr26 <- savitzkyGolay(mXr24, m=2, p=2, w=3)
if (species == "Pinus_patula") {
  mXr25 <- savitzkyGolay(mXr25, m=0, p=2, w=15) # smoothing
  mXr26 <- savitzkyGolay(mXr26, m=0, p=2, w=15) # smoothing
}

#-------CHOOSE WHICH METHOD TO USE-----------
mXr_mod <- mXr05 # choose which one to use
main_title = "SG + 1st derivarive + SNV"
file_name <- gsub(" ", "_", main_title)
file_name <- paste0(file_name, ".png")

# plot spectra and mean
matplot(colnames(mXr_mod), t(mXr_mod), type="l", lty=1, xlim=c(1800, 800), 
        xlab="wavenumber [cm-1]", ylab="absorbance",
        col=rgb(0.5, 0.5, 0.5, 0.8))
matlines(colnames(mXr_mod), y = colMeans(mXr_mod), col = "red", lwd = 2)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## helpful PLS-DA tutorial: https://mdatools.com/docs/plsda--calibration.html

### TRAIN / TEST
# NAES sampling
# k = 70% training, remaining 30% validation
# for pc, you can also insert the var_explained, e.g. pc=0.95
#----------
# define class labels
class_labels <- file_meta$ring

# subset data for each class
mod_EW <- mXr_mod[class_labels == "EW", ]
mod_LW <- mXr_mod[class_labels == "LW", ]

# apply NAES separately to each class
n_LW <- table(file_meta$ring)["LW"][[1]] # calc number of target LW spectra
set.seed(127)
# 70% train, 30% test
split_EW <- naes(mod_EW, k = ceiling(n_LW * 0.7), pc = 0.95, iter.max = 10, method = 0, .center = TRUE, .scale = TRUE)
split_LW <- naes(mod_LW, k = ceiling(n_LW * 0.7), pc = 0.95, iter.max = 10, method = 0, .center = TRUE, .scale = TRUE)

# relative ixs (within each subset)
csc_EW <- split_EW$model
cst_EW <- split_EW$test
csc_LW <- split_LW$model
cst_LW <- split_LW$test

# map back to full dataset
ix_EW <- which(class_labels == "EW")
ix_LW <- which(class_labels == "LW")

csc_EW_abs <- ix_EW[csc_EW]
cst_EW_abs <- ix_EW[cst_EW]
csc_LW_abs <- ix_LW[csc_LW]
cst_LW_abs <- ix_LW[cst_LW]

# combine absolute ixs
csc <- c(csc_EW_abs, csc_LW_abs)
cst  <- c(cst_EW_abs, cst_LW_abs)

# (optionally) check class distribution to ensure balanced dataset
table(class_labels[csc])
table(class_labels[cst])
#-------------------------------------------------------------------------------

# PCA to be able to plot the distribution of the samples
pc <- prcomp(mXr_mod, center = T, scale. = T, rank. = 4) # rank = n of PCs

par(pty="s") # s = squared, m=maximized (plot space filled)
file_path <- file.path("data", "out", species, "PLS-DA_calibration", "cal_split", file_name)
png(file_path, width = 3600, height = 3600, res = 300)
matplot(pc$x[, 1], pc$x[, 2], pch = 16, col="black",
        xlab = "", ylab = "", main = main_title, cex.main = 3, cex = 2, cex.axis = 2)
matpoints(pc$x[csc, 1], pc$x[csc, 2], pch = 16, col = "red", cex = 2) # cal samples in red
legend("topright", c("model", "test"), pch = 16, col=c("red", "black"), cex = 2.5)
# result looks fine: good spread of variation in train as well as test data set
dev.off()
#-------------------------------------------------------------------------------

Xc = mXr_mod[csc, ]
Xv = mXr_mod[cst, ]

cc_all = as.factor(file_meta$ring[csc])
cv_all = as.factor(file_meta$ring[cst])
cc_lw = cc_all == "LW" # only TRUE/FALSE needed
cv_lw = cv_all == "LW"

set.seed(127)
m.ring = mdatools::plsda(Xc, cc_lw, ncomp = 25, cv = 10, classname = "LW", scale = TRUE)
summary(m.ring)
getConfusionMatrix(m.ring$calres)

#plotPredictions(m.ring, ncomp=8)
#plotMisclassified(m.ring, ylim=c(0, 0.5), xlim=c(1, 25))
#plotSensitivity(m.ring)
#plotSpecificity(m.ring)
#plotRegcoeffs(m.ring, ncomp = 1, show.ci = FALSE, xlim=c(1100, 0))
#plotProbabilities(res, ylim=c(-0.1, 2))

# bias-variance trade-off
file_path <- file.path("data", "out", species, "PLS-DA_calibration", "n_lv_sel", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
par(pty = "m", cex = 1.8, mar = c(3, 3, 2, 2))
plotRMSE(m.ring, xlim = c(1, 25), cex = 1.2, show.grid = F, xlab = "", ylab = "", ylim = c(0, 1.2), 
         #cex.axis = 5, cex.names = 5, # plots without x and y-labs to make grid in latex
         main = main_title, col = c('darkorchid1', 'blue2'))
# save with 1600 minimum
dev.off()


#-------CHOOSE HOW MANY LVs SHOULD BE USED-----------
nlv = 5 # choose best nlv

#-------------------------------------------------------------------------------

# prediction
res = predict(m.ring, Xv, cv_lw)
summary(res)

getConfusionMatrix(res, ncomp=nlv)

# plot predictions
file_path <- file.path("data", "out", species, "PLS-DA_calibration", "pred_test", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
par(pty = "m", cex = 1.8, mar = c(3, 2, 2, 1))
plotPredictions(res, ncomp = nlv, ylim = c(0.8, 2.4), cex = 1.8, pch = 124, main = "", lab.col = "black", lab.cex = 1.2)
# ylab="Ring class", xlab="Test set spectra"
title(main = main_title, line = 0.7, cex.main = 1.6)

dev.off()
#-------------------------------------------------------------------------------

# calculate variable importance using VIP function
# requires PLS object and optimal number of latent variables
vip <- vipscores(m.ring, ncomp = nlv)

file_path <- file.path("data", "out", species, "PLS-DA_calibration", "vip", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
par(mgp = c(3, 2, 0))
matplot(colnames(mXr_mod), vip, type="l", lty=1, xlim=c(1800, 800), xlab = "", ylab = "",
        cex.axis = 2.5, tck = -0.02,
        main = main_title, cex.main = 3)
abline(h=1, col="red")

dev.off()
#-------------------------------------------------------------------------------

# performance metrics
TP = getConfusionMatrix(res, ncomp=nlv)["LW", "LW"]
TN = getConfusionMatrix(res, ncomp=nlv)["None", "None"]
FP = getConfusionMatrix(res, ncomp=nlv)["LW", "None"]
FN = getConfusionMatrix(res, ncomp=nlv)["None", "LW"]
OA = round((TP + TN) / sum(TP, TN, FP, FN), 2)  # overall accuracy
Recall_LW = round((TP / (TP + FN)), 2)  # recall for class "LW"
Precision_LW = round((TP / (TP + FP)), 2)  # precision for class "LW"
Specificity_LW = round((TN / (TN + FP)), 2)  # specificity for class "LW"

Recall_EW <- round((TN / (TN + FP)), 2)  # sensitivity for EW
Precision_EW <- round((TN / (TN + FN)), 2)  # precision for EW
Specificity_EW <- round((TP / (TP + FN)), 2) # specificity for EW


metrics_row <- data.frame(
  method = main_title,
  nlv = nlv,
  
  Recall_LW = Recall_LW,
  Recall_EW = Recall_EW,

  Precision_LW = Precision_LW,
  Precision_EW = Precision_EW,
  
  Specificity_LW = Specificity_LW,
  Specificity_EW = Specificity_EW,
  OA = OA
)

# creation of results table for different preprocessing methods (the following adds 1 row for the current method used)
# read existing data (if file exists)
if (file.exists(file.path("data", "out", species, "classification_metrics.xlsx"))) {
  existing_data <- read_xlsx(file.path("data", "out", species, "classification_metrics.xlsx"))
} else {
  existing_data <- NULL
}

# combine with new row
metrics_row <- as.data.frame(metrics_row)
updated_data <- rbind(existing_data, metrics_row)

write_xlsx(updated_data, path = file.path("data", "out", species, paste0(species, "_classification_metrics.xlsx")))
#-------------------------------------------------------------------------------

### print accuracy metrics
# recall
cat("Recall (LW):", Recall_LW, "\n")
cat("Recall (EW):", Recall_EW, "\n")
# precision
cat("Precision (LW):", Precision_LW, "\n")
cat("Precision (EW):", Precision_EW, "\n")
# specifity
cat("Specificity (LW):", Specificity_LW, "\n")
cat("Specificity (EW):", Specificity_EW, "\n")
# OA
cat("Overall Accuracy (OA):", OA, "\n")
