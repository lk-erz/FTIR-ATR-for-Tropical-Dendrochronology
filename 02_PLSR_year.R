#==============================================================================
# Supplementary R scripts for mater's thesis
# Title: Potential of FTIR spectroscopy for tropical dendrochronology
# Author: Lukas Erzfeld
# Insitution: University of Leipzig
# Date: August 26, 2025
#===============================================================================

# required packages may have to be installed primarily
# install.packages(c("prospectr", "pls", "plsVarSel", "readxl", "writexl"))
library(prospectr)
library(pls)
library(plsVarSel)
library(readxl) 
library(writexl)
source("Utils.R")
#-------------------------------------------------------------------------------
# NOTE: number of latent variables used in the models have to be chosen in lines 248 and 335 
#-------------------------------------------------------------------------------

# out path of 00_PLS-DA_data_prep_master.R
path = "data/FTIR_Spectra/PLSR/"
species = "Pinus_patula"
filename_meta <- paste0(species, "_file_meta_PLSR.csv")
filename_dat <- paste0(species, "_data_mat_PLSR.csv")

path_meta <- file.path(path, filename_meta)
path_spectra <- file.path(path, filename_dat)

# read saved data
file_meta <- read.csv(path_meta)
data <- read.csv(path_spectra)
colnames(data) <- gsub("^X", "", colnames(data))
colnames(data) <- as.numeric(colnames(data))
wn <- as.numeric(colnames(data))
#-------------------------------------------------------------------------------

# a few spectra show very low absorption values at band 1025.9445,
# which usually exhibits the highest absorption values --> remove the lowest few spec
# (10 for B. simaruba, 4 for M. tinctoria, 2 for P. patula)
ref_band <- "1025.9445"
if (species == "Bursera_simaruba") {
  n_rm <- 10
} else if (species == "Maclura_tinctoria") {
  rm_outlier <- c(100, 111, 113, 114, 180, 189) # PCA outliers
  data <- data[-rm_outlier, ]
  file_meta <- file_meta[-rm_outlier, ]
  n_rm <- 4
  
  
  rm_ix_high <- order(data[["1178.2933"]], decreasing = TRUE)[1:n_rm]
  data <- data[-rm_ix_high, ]
  file_meta <- file_meta[-rm_ix_high, ]
} else if (species == "Pinus_patula") {
  n_rm <- 2
  rm_ix_high <- tail(order(data[["1620.8765"]]), n_rm)
  data <- data[-rm_ix_high, ]
  file_meta <- file_meta[-rm_ix_high, ]
}
rm_ix <- order(data[[ref_band]])[1:n_rm]
data <- data[-rm_ix, ]
file_meta <- file_meta[-rm_ix, ]
#-------------------------------------------------------------------------------

# check distribution of response variable
growth_years <- as.numeric(file_meta$growth_years)
hist(growth_years, freq=F, 
     breaks=c(seq(0, 100, by=5)),
     xlim=c(0, 100),# ylim=c(0, 0.02),
     col="blue", xlab="Years of growth", ylab="Probab. density")
lines(density(growth_years))

boxplot(growth_years, main="Boxplot of growth years")#, ylim=c(1945, 2025))

qqnorm(growth_years) # looks okay
qqline(growth_years, col="red")
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
mXrSG1 <- savitzkyGolay(mXrSmoothed, m=1, p=1, w=3)
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
mXrBaseBand <- manual_baseline_correction(data, wn=as.numeric(wn4000), anchors = c(3750, 1780, 1485, 1185, 820), window = 50) # for Pinus p. 1530, not 1485
mXrBase <- mXrBase[, fingerprint_cols]
mXrBaseBand <- mXrBaseBand[, fingerprint_cols]

mXr11 <- standardNormalVariate(mXrBase)
mXr12 <- standardNormalVariate(mXrBaseBand)

# BASELINE + BAND NORMALIZATION 1424
mXr13 <- normalize_to_band(mXrBase, wn, band = 1424)
mXr13 <- mXr13[, apply(mXr13, 2, sd) > 0] # remove this band, otherwise PLSR returns NaN
mXr14 <- normalize_to_band(mXrBaseBand, wn, band = 1424)
mXr14 <- mXr14[, apply(mXr14, 2, sd) > 0]

# BASELINE + MEAN SCALING
mXr15 <- sweep(mXrBase, MARGIN=1, apply(mXrBase, 1, function(x) (mean(x))), FUN="/")
mXr16 <- sweep(mXrBaseBand, MARGIN=1, apply(mXrBaseBand, 1, function(x) (mean(x))), FUN="/")

# BASELINE + BAND NORMALIZATION 1026
mXr17 <- normalize_to_band(mXrBase, wn, band = 1026)
mXr17 <- mXr17[, apply(mXr17, 2, sd) > 0] # remove this band, otherwise PLSR returns NaN
mXr18 <- normalize_to_band(mXrBaseBand, wn, band = 1026)
mXr18 <- mXr18[, apply(mXr18, 2, sd) > 0] # remove this band, otherwise PLSR returns NaN

# DERIVATIVES + SNV (poly 3)
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
mXr_mod <- mXr00 # choose which one to use
main_title = "band auschecken"
coef_leg <- "topleft"
coef_leg <- "bottomleft"
#coef_leg <- "bottomright"
file_name <- gsub(" ", "_", main_title)
file_name <- paste0(file_name, ".png")

# plot spectra and mean
matplot(colnames(mXr_mod), t(mXr_mod), type="l", lty=1, xlim=c(1800, 800), 
        xlab="wavenumber [cm-1]", ylab="absorbance",
        col=rgb(0.5, 0.5, 0.5, 0.8))
matlines(colnames(mXr_mod), y = colMeans(mXr_mod), col = "red", lwd = 2)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### TRAIN / TEST
# NAES sampling
# k = 70% training, remaining 30% validation
# for pc, you can also insert the var_explained, e.g. pc=0.95
n_spec <- nrow(data)
set.seed(127)
split <- naes(mXr_mod, k = n_spec * 0.7, pc = 0.95, iter.max = 10, method = 0, .center = TRUE, .scale = TRUE)

csc <- split$model # row numbers calibration set
cst <- split$test # row numbers test set
#-------------------------------------------------------------------------------

# PCA to be able to plot the distribution of the samples
pc <- prcomp(mXr_mod, center = T, scale. = T, rank. = 4) # rank = n of PCs

par(pty="s")

test_idx <- setdiff(seq_len(nrow(pc$x)), csc) # only for plotting
file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_all_bands", "cal_split", file_name)
png(file_path, width = 3600, height = 3600, res = 300)
matplot(pc$x[test_idx, 1], pc$x[test_idx, 2], pch = 16, col = "black", 
        xlab = "", ylab = "", main = main_title, cex.main = 3, cex = 2, cex.axis = 2)
matpoints(pc$x[csc, 1], pc$x[csc, 2], pch = 16, col = "red", cex = 2)
legend("topright", c("model", "test"), pch = 16, col = c("red", "black"), cex = 2.5)
# result looks fine: good spread of variation in train as well as test data set
dev.off()
#-------------------------------------------------------------------------------

### define response and predictors
Y <- as.matrix(as.numeric(file_meta$growth_years[csc])) # plsr function needs matrix
X <- as.matrix(mXr_mod[csc, ]) # slice only for training

set.seed(127)
mod <- plsr(Y~X, ncomp = 25, method = "oscorespls", 
            validation = "CV", # method of validation, see documentation 
            segments = 10, # for B. simaruba: 691 spectra, so ~69 per segment
            jackknife = T, # for significance assessment of bands
            scale = T 
            )

#par(pty="m")
#matplot(colnames(mXr_mod), mod$coefficients[,,], type="l", lty=1, xlim=c(1800, 800)) # plot all latent variables ([,,])
#matplot(colnames(mXr_mod), mod$coefficients[,,2], type="l", lty=1, xlim=c(1800, 800))

# RMSE
rmse <- RMSEP(mod, "all")
rms <- rmse$val[1:2, , -1]
colnames(rms) <- c(1:ncol(rms))
rownames(rms) <- c("RMSEcal", "RMSEval")

# bias-variance trade-off
file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_all_bands", "n_lv_sel", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
#barplot(rms[2, 1:25], space=0, col='blue2', ylim=c(0, max(rms[2, ])), las=2, cex.names=0.6, xlab="# latent vectors", ylab="RMSE")
barplot(rms[2, 1:25], space = 0, col = 'blue2', ylim = c(0, ceiling(max(rms[2, ]) / 5) * 5), las = 2, xlab = "", ylab = "", 
        cex.axis = 2.5, cex.names = 2,
        main = main_title, cex.main = 3)
barplot(rms[1, 1:25], space = 0, col = 'darkorchid1', add = T, axes = F, axisnames = F)
legend("topright", legend=c("cal", "val"), pch = 15, col = c('darkorchid1', 'blue2'), bty = "n", horiz = T, cex = 3, pt.cex = 3.5)
      # save with 1600 minimum
dev.off()


#-------CHOOSE nLV: LOWEST RMSE WHILE CONSIDERING MODEL COMPLEXITY-----------
nlv <- 5

#-------------------------------------------------------------------------------

fit.val <- mod$validation$pred[, , nlv]
fit.cal <- mod$fitted.values[, , nlv]

file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_all_bands", "pred_cal_CV", file_name)
png(file_path, width = 3600, height = 3600, res = 300)
par(mgp = c(3, 2, 0))
#plot(growth_years[csc], fit.cal, xlab="Counted distance to pith (years)", ylab="Predicted distance to pith (years)", pch=19, col='blue2', xlim=c(0, 100), ylim=c(0, 100))
plot(growth_years[csc], fit.cal, xlab = "", ylab = "", pch = 19, col = 'blue2', xlim = c(-5, 80), ylim = c(-5, 80), 
     cex.axis = 3, tck = -0.02, cex = 1.5) # plots without x and y-labs to make grid in latex
title(main = main_title, cex.main = 3)
abline(0, 1, col = "black", lty = 2, lwd = 2)
points(growth_years[csc], fit.val, pch = 19, col = 'darkorchid1', cex = 1.5)
pred.obs.cal <- lm(fit.cal ~ growth_years[csc])
pred.obs.val <- lm(fit.val ~ growth_years[csc])
abline(pred.obs.cal, col='blue2', lwd = 2)
abline(pred.obs.val, col='darkorchid1', lwd = 2)
leg <- c(
  paste0("RMSEcal = ", round(rms[1, nlv], 2)),
  paste0("RMSECV = ", round(rms[2, nlv], 2)),
  bquote(R^2*"cal ="~.(round(summary(pred.obs.cal)$r.squared, 2))),
  bquote(R^2*"CV ="~.(round(summary(pred.obs.val)$r.squared, 2)))
)
legend("topleft", legend=leg, text.col=c('blue2', 'darkorchid1', 'blue2', 'darkorchid1'), bty="n", cex = 3)
# then export with 1600 and then crop
dev.off()
#-------------------------------------------------------------------------------

# regression coefficients
reg.coef <- coef(mod, ncomp=nlv)
p <- 0.05 # define the threshold for the level of significance
jt <- suppressWarnings(jack.test(mod, ncomp=nlv))$pvalues
sig <- jt < p


wn <- as.numeric(colnames(mXr_mod))
file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_all_bands", "coefficients", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
par(mgp = c(3, 2, 0))
#plot(wn, reg.coef, type="h", xlab="Wavenumber / cm-1", ylab="Coefficient", xlim=c(1800, 800))#, ylim=c(-0.45, 0.45), lab=c(10, 5, 5))
# plots without x and y-labs to make grid in latex
plot(wn, reg.coef, type = "h", xlab = "", ylab = "", xlim = c(1800, 800), cex.axis = 2.5, tck = -0.02) 
title(main = main_title, cex.main = 3)
points(wn, reg.coef, pch=c(1, 19)[sig+1])
grid()
legend(coef_leg, legend = paste0(c('p<', 'p>='), p), pch = c(19, 1), bty = "n", cex = 2.5)
#abline(v=c(1593, 1505, 1460, 1424, 1267, 1225, 960, 859, 818), col="green") # important bands for lignin content
  # save!
dev.off()
#-------------------------------------------------------------------------------

# optimizing the model through backward selection of spectral bands
wn2 <- wn[sig] # only significant bands
wn2_char <- as.character(wn2)
X2 <- X[, wn2_char, drop=FALSE]

set.seed(127)
mod2 <- plsr(Y~X2, ncomp = 25, method="oscorespls", 
            validation="CV", # method of validation, see documentation 
            segments = 10,
            jackknife = T,
            scale = T 
)

rmse2 <- RMSEP(mod2, "all")
rms2 <- rmse2$val[1:2, , -1]
colnames(rms2) <- c(1:ncol(rms2))
rownames(rms2) <- c("RMSEcal", "RMSEval")


file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_backwards_selection", "n_lv_sel", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
# uncomment for axis labels
#barplot(rms2[2, 1:25], space=0, col='blue2', ylim=c(0, max(rms2[2, ])), las=2, xlab="# latent vectors", ylab="RMSE (years)")
barplot(rms2[2, 1:25], space = 0, col = 'blue2', ylim = c(0, ceiling(max(rms2[2, ]) / 5) * 5), las = 2, xlab = "", ylab = "", 
        cex.axis = 2.5, cex.names = 2, 
        main = main_title, cex.main = 3)
barplot(rms2[1, 1:25], space = 0, col = 'darkorchid1', add = T, axes = F, axisnames = F)
legend("topright", legend = c("cal", "val"), pch = 15, col = c('darkorchid1', 'blue2'), bty = "n", horiz = T, cex = 3, pt.cex = 3.5) 
      # save with 1600
dev.off()


#-------CHOOSE nLV: LOWEST RMSE WHILE CONSIDERING MODEL COMPLEXITY-----------
nlv2 <- 7

#-------------------------------------------------------------------------------
fit.val2 <- mod2$validation$pred[, , nlv2]
fit.cal2 <- mod2$fitted.values[, , nlv2]

file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_backwards_selection", "pred_cal_CV", file_name)
png(file_path, width = 3600, height = 3600, res = 300)
par(mgp = c(3, 2, 0))
# uncomment for axis labels
#plot(growth_years[csc], fit.cal2, xlab="Counted distance to pith (years)", ylab="Predicted distance to pith (years)", pch=19, col='blue2', xlim=c(-5, 100), ylim=c(-5, 100), cex.axis = 1.8)
plot(growth_years[csc], fit.cal2, xlab = "", ylab = "", pch = 19, col = 'blue2', xlim = c(-5, 80), ylim = c(-5, 80),
     cex.axis = 3, tck = -0.02, cex = 1.5)
title(main = main_title, cex.main = 3)

abline(0, 1, col = "black", lty = 2, lwd = 2)
points(growth_years[csc], fit.val2, pch = 19, col = 'darkorchid1', cex = 1.5)
pred.obs.cal2 <- lm(fit.cal2 ~ growth_years[csc])
pred.obs.val2 <- lm(fit.val2 ~ growth_years[csc])
abline(pred.obs.cal2, col = 'blue2', lwd = 2)
abline(pred.obs.val2, col = 'darkorchid1', lwd = 2)
leg <- list(
  paste0("RMSEcal = ", round(rms2[1, nlv], 2)),
  paste0("RMSECV = ", round(rms2[2, nlv], 2)),
  bquote(R^2*"cal ="~.(round(summary(pred.obs.cal2)$r.squared, 2))),
  bquote(R^2*"CV ="~.(round(summary(pred.obs.val2)$r.squared, 2)))
)
legend("topleft", legend = leg, text.col = c('blue2', 'darkorchid1', 'blue2', 'darkorchid1'), bty = "n", cex = 3)
dev.off()

#-------------------------------------------------------------------------------
# VIP
par(pty="m")
# calculate variable importance using VIP function
# requires PLS object and optimal number of latent variables
vi <- VIP(mod, opt.comp=nlv)
# result is vector with 1025 values (so importance assigned to each band)
file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_all_bands", "vip", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
par(mgp = c(3, 2, 0)) 
#matplot(colnames(mXr_mod), vi, type="l", lty=1, xlim=c(1800, 800), xlab="wavenumber [cm-1]", ylab="VIP score")
matplot(colnames(mXr_mod), vi, type="l", lty=1, xlim=c(1800, 800), xlab = "", ylab = "", 
  cex.axis = 2.5, tck = -0.02, # plots without x and y-labs to make grid in latex
  main = main_title, cex.main = 3)
abline(h = 1, col="red") # mostly between 0.8-1.2 (mostly used is 1), in the documentation is some info
dev.off()

file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_backwards_selection", "vip", file_name)
png(file_path, width = 3600, height = 2276, res = 300)
par(mgp = c(3, 2, 0))
# of the backwards selected model
vi2 <- VIP(mod2, opt.comp=nlv2)
# result is vector with 1025 values (so importance assigned to each band)
matplot(wn2, vi2, type = "l", lty = 1, xlim = c(1800, 800), xlab = "", ylab = "", 
        cex.axis = 2.5, tck = -0.02, # plots without x and y-labs to make grid in latex
        main = main_title, cex.main = 3)
abline(h = 1, col = "red")
dev.off()
#-------------------------------------------------------------------------------

# prediction
Ynew <- as.matrix(as.numeric(file_meta$growth_years[cst]))
Xnew <- as.matrix(mXr_mod[cst,])

Ypred <- predict(mod, ncomp=nlv, newdata=Xnew)

# quantify outcome
np <- nrow(Ynew)
SST <- sum((Ypred[, 1, 1] - mean(Ynew)) ^ 2)
SSE <- sum((Ynew - Ypred[, 1, 1]) ^ 2)
R2 <- round(1 - (SSE / SST), digits=4)
R2

RMSE <- round(sqrt(sum((Ypred[, 1, 1] - Ynew) ^ 2) / np), 
              digits=4)
RMSE

RPD <- sd(Ynew) / RMSE
RPD

# plot prediction
par(pty="s")

file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_all_bands", "pred_test", file_name)
png(file_path, width = 3600, height = 3600, res = 300)
par(mgp = c(3, 2, 0))
#matplot(Ynew, Ypred[, 1, 1], pch = 16, 
#        xlab = "Counted distance to pith (years)", ylab="Predicted distance to pith (years)", 
#        xlim = c(-5, 100), ylim = c(-5, 100),
#        cex.axis = 1.8)
matplot(Ynew, Ypred[, 1, 1], pch = 16, 
        xlab = "", ylab = "", # plots without x and y-labs to make grid in latex
        xlim = c(-5, 80), ylim = c(-5, 80),
        cex.axis = 3, tck = -0.02, cex = 2)
title(main = main_title, cex.main = 3.5)
abline(0, 1, lty = 2, lwd = 3)
vm <- lm(Ypred[, 1, 1] ~ Ynew[, 1])
abline(vm, lty = 1, lwd = 3)
coef_vm <- coef(vm)
intercept <- round(coef_vm[1], 2)
slope <- round(coef_vm[2], 2)
leg <- list(
  bquote(R^2*" ="~.(round(R2, 2))),
  paste0("RMSE = ", (round(RMSE, 2))),
  paste0("RPD = ", (round(RPD, 2))),
  bquote(y == .(slope) * x + .(intercept))
)
legend("topleft", legend=leg, bty="n", cex = 3.5)
dev.off()


metrics_row <- data.frame(
  method = main_title,
  nlv = nlv,
  
  RMSEC = round(rms[1, nlv], 2),
  R2C = round(summary(pred.obs.cal)$r.squared, 2),
  
  RMSECV = round(rms[2, nlv], 2),
  R2CV = round(summary(pred.obs.val)$r.squared, 2),
  
  RMSE = round(RMSE, 2),
  R2 = round(R2, 2),
  RPD = round(RPD, 2)
)

# creation of results table for different preprocessing methods (the following adds 1 row for the current method used)
# read existing data (if file exists)
if (file.exists(file.path("data", "out", species, "regression_metrics.xlsx"))) {
  existing_data <- read_xlsx(file.path("data", "out", species, "regression_metrics.xlsx"))
} else {
  existing_data <- NULL
}

# combine with new row
metrics_row <- as.data.frame(metrics_row)
updated_data <- rbind(existing_data, metrics_row)

write_xlsx(updated_data, path = file.path("data", "out", species, "regression_metrics.xlsx"))
#-------------------------------------------------------------------------------

### same for model with backward selected bands
# prediction
Ynew2 <- as.matrix(as.numeric(file_meta$growth_years[cst]))
Xnew2 <- as.matrix(mXr_mod[cst,])
Xnew2 <- Xnew2[, wn2_char, drop=FALSE]

Ypred2 <- predict(mod2, ncomp=nlv2, newdata=Xnew2)

# quantify outcome
np2 <- nrow(Ynew2)
SST2 <- sum((Ypred2[, 1, 1] - mean(Ynew2)) ^ 2)
SSE2 <- sum((Ynew2 - Ypred2[, 1, 1]) ^ 2)
R22 <- round(1 - (SSE2 / SST2), digits=4)
R22

RMSE2 <- round(sqrt(sum((Ypred2[, 1, 1] - Ynew2) ^ 2) / np2), 
               digits=4)
RMSE2

RPD2 <- sd(Ynew2) / RMSE2
RPD2

# plot prediction
par(pty="s")

file_path <- file.path("data", "out", species, "PLSR_calibration", "mod_backwards_selection", "pred_test", file_name)
png(file_path, width = 3600, height = 3600, res = 300)
par(mgp = c(3, 2, 0))
#matplot(Ynew2, Ypred2[, 1, 1], pch = 16, 
#        xlab = "Counted distance to pith (years)", ylab="Predicted distance to pith (years)", 
#        xlim = c(-5, 100), ylim = c(-5, 100),
#        cex.axis = 1.8)
matplot(Ynew2, Ypred2[, 1, 1], pch = 16, 
        xlab = "", ylab = "", # plots without x and y-labs to make grid in latex
        xlim = c(-5, 80), ylim = c(-5, 80),
        cex.axis = 3, tck = -0.02, cex = 2)
title(main = main_title, cex.main = 3.5)
abline(0, 1, lty=2, lwd = 3)
vm <- lm(Ypred2[, 1, 1] ~ Ynew2[, 1])
abline(vm, lty=1, lwd = 3)
coef_vm <- coef(vm)
intercept <- round(coef_vm[1], 2)
slope <- round(coef_vm[2], 2)
leg <- list(
  bquote(R^2*" ="~.(round(R22, 2))),
  paste0("RMSE = ", (round(RMSE2, 2))),
  paste0("RPD = ", (round(RPD2, 2))),
  bquote(y == .(slope) * x + .(intercept))
)
legend("topleft", legend=leg, bty="n", cex = 3.5)
dev.off()


metrics_row2 <- data.frame(
  method = main_title,
  n_bands = length(wn2),
  nlv2 = nlv2,
  
  RMSEC2 = round(rms2[1, nlv2], 2),
  R2C2 = round(summary(pred.obs.cal2)$r.squared, 2),
  
  RMSECV2 = round(rms2[2, nlv2], 2),
  R2CV2 = round(summary(pred.obs.val2)$r.squared, 2),
  
  RMSE2 = round(RMSE2, 2),
  R22 = round(R22, 2),
  RPD2 = round(RPD2, 2)
)

# creation of results table for different preprocessing methods (the following adds 1 row for the current method used)
# read existing data (if file exists)
if (file.exists(file.path("data", "out", species, "regression_metrics_backward_selection.xlsx"))) {
  existing_data <- read_xlsx(file.path("data", "out", species, "regression_metrics_backward_selection.xlsx"))
} else {
  existing_data <- NULL
}

# combine with new row
metrics_row2 <- as.data.frame(metrics_row2)
updated_data2 <- rbind(existing_data, metrics_row2)

write_xlsx(updated_data2, path = file.path("data", "out", species, "regression_metrics_backward_selection.xlsx"))
