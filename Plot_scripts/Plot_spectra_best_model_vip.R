# set wd to root dir
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # assuming use of RStudio
parent_dir <- dirname(current_dir)
setwd(parent_dir) # or set wd manually here

library(prospectr)
library(patchwork)
library(tidyverse)
source("Utils.R")
#-------------------------------------------------------------------------------
# out path of 00_PLS-DA_data_prep_master.R
path = "data/FTIR_Spectra/PLSR/"
species = "Bursera_simaruba"
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
     xlim=c(0, 100), ylim=c(0, 0.02),
     col="blue", xlab="Years of growth", ylab="Probab. density")
lines(density(growth_years))

boxplot(growth_years, main="Boxplot of growth years")#, ylim=c(1945, 2025))

qqnorm(growth_years) # looks okay for all
qqline(growth_years, col="red")
#-------------------------------------------------------------------------------
# slice data to fingerprint region 800 - 1800 cm-1
fingerprint_cols <- wn >= 800 & wn <= 1800
wn4000 <- wn
wn <- wn[fingerprint_cols]
spec_fingerprint <- data[, fingerprint_cols]

### preprocessing ###
mXr00 <- as.matrix(spec_fingerprint)
mXr01 <- sweep(mXr00, MARGIN=1, apply(mXr00, 1, function(x) (mean(x))), FUN="/") # mean scaling
mXr001 <- scale(mXr00, center = TRUE, scale = TRUE) # z-transform
mXr02 <- standardNormalVariate(mXr00)
mXr03 <- prospectr::detrend(mXr00, wn)

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
mXrBaseBand <- manual_baseline_correction(data, wn=as.numeric(wn4000), anchors = c(3750, 1780, 1485, 1185, 820), window = 50) # for P. patula 1530, not 1485
mXrBase <- mXrBase[, fingerprint_cols]
mXrBaseBand <- mXrBaseBand[, fingerprint_cols]

mXr11 <- standardNormalVariate(mXrBase)
mXr12 <- standardNormalVariate(mXrBaseBand)

# BASELINE + BAND NORMALIZATION 1424
mXr13 <- normalize_to_band(mXrBase, wn, band = 1424)
mXr13 <- mXr13[, apply(mXr13, 2, sd) > 0]
mXr14 <- normalize_to_band(mXrBaseBand, wn, band = 1424)
mXr14 <- mXr14[, apply(mXr14, 2, sd) > 0]

# BASELINE + MEAN SCALING
mXr15 <- sweep(mXrBase, MARGIN=1, apply(mXrBase, 1, function(x) (mean(x))), FUN="/")
mXr16 <- sweep(mXrBaseBand, MARGIN=1, apply(mXrBaseBand, 1, function(x) (mean(x))), FUN="/")

# BASELINE + BAND NORMALIZATION 1026
mXr17 <- normalize_to_band(mXrBase, wn, band = 1026)
mXr17 <- mXr17[, apply(mXr17, 2, sd) > 0]
mXr18 <- normalize_to_band(mXrBaseBand, wn, band = 1026)
mXr18 <- mXr18[, apply(mXr18, 2, sd) > 0]

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

# plotting function
plot_spectra_with_bands <- function(data_raw, data_best_mod, bands_df) {
  
  # --- raw spectra ---
  df_raw <- as.data.frame(data_raw) %>%
    mutate(Sample = row_number()) %>%
    pivot_longer(-Sample, names_to = "Wavenumber", values_to = "Absorbance") %>%
    mutate(Wavenumber = as.numeric(Wavenumber))
  
  df_raw_mean <- df_raw %>%
    group_by(Wavenumber) %>%
    summarise(Absorbance = mean(Absorbance), .groups = "drop")
  
  p_raw <- ggplot() +
    geom_line(data = df_raw, aes(x = Wavenumber, y = Absorbance, group = Sample),
              color = rgb(0.5, 0.5, 0.5, 0.6), linewidth = 0.8) +
    geom_line(data = df_raw_mean, aes(x = Wavenumber, y = Absorbance),
              color = "red", linewidth = 0.6) +
    geom_vline(data = bands_df, aes(xintercept = Wavenumber, color = Type),
               linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_text(data = bands_df,
              aes(x = Wavenumber, y = Inf, label = Wavenumber, color = Type, fontface = "bold", vjust = Vjust),
              angle = 90, hjust = 1.5, size = 4.5, show.legend = FALSE) +
    scale_color_manual(values = c("VIP" = "black", "Positive" = "black", "Negative" = "black")) +
    scale_x_reverse(limits = c(1800, 800),
                    breaks = seq(1800, 800, by = -200),
                    expand = c(0, 20)) +
    labs(x = NULL, y = "Absorbance (a.u.)") +
    annotate("text", x = max(df_raw$Wavenumber), y = max(df_raw$Absorbance),
             label = "a", hjust = -2.0, vjust = 0.5, size = 6, fontface = "bold") +
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_line(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1.3),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # --- best model spectra ---
  df_best_mod <- as.data.frame(data_best_mod) %>%
    mutate(Sample = row_number()) %>%
    pivot_longer(-Sample, names_to = "Wavenumber", values_to = "Absorbance") %>%
    mutate(Wavenumber = as.numeric(Wavenumber))
  
  df_best_mod_mean <- df_best_mod %>%
    group_by(Wavenumber) %>%
    summarise(Absorbance = mean(Absorbance), .groups = "drop")
  
  p_best_mod <- ggplot() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3, alpha = 0.8) +
    geom_line(data = df_best_mod, aes(x = Wavenumber, y = Absorbance, group = Sample),
              color = rgb(0.5, 0.5, 0.5, 0.6), linewidth = 0.8) +
    geom_line(data = df_best_mod_mean, aes(x = Wavenumber, y = Absorbance),
              color = "red", linewidth = 0.6) +
    geom_vline(data = bands_df, aes(xintercept = Wavenumber, color = Type),
               linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_text(data = bands_df,
              aes(x = Wavenumber, y = Inf, label = Wavenumber, color = Type, fontface = "bold", vjust = Vjust),
              angle = 90, hjust = 1.5, size = 4.5, show.legend = FALSE) +
    scale_color_manual(values = c("VIP" = "black", "Positive" = "black", "Negative" = "black")) +
    scale_x_reverse(limits = c(1800, 800),
                    breaks = seq(1800, 800, by = -200),
                    expand = c(0, 20)) +
    scale_y_continuous(breaks = c(0)) +
    labs(x = "Wavenumber (cm⁻¹)", y = "Absorbance (a.u.)") +
    annotate("text", x = max(df_raw$Wavenumber), y = max(df_best_mod$Absorbance),
             label = "b", hjust = -2.0, vjust = 0.5, size = 6, fontface = "bold") +
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_text(size = 18, color = "black"),
      axis.ticks.x = element_line(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text.y = element_text(color = "black"),
      axis.ticks.y = element_line(),
      axis.title.y = element_text(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1.3),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # --- combine plots ---
  combined_plot <- p_raw / p_best_mod +
    plot_layout(heights = c(1, 1)) & 
    theme(plot.margin = margin(1, 0, 1, 0))  # remove space
  
  return(combined_plot)
}

#-------------------------------------------------------------------------------------------
### Bursera simaruba
# define important bands
important_wavenumbers_vip <- c(1022, 1069, 1228, 1481)
important_wavenumbers_coeff_positive <- c(1327, 1786, 1420, 956)
important_wavenumbers_coeff_negative <- c(1456, 1057, 893, 828)
label_right_bands <- c(1057, 1022)

# df creation
vip_df <- tibble(Wavenumber = important_wavenumbers_vip, Type = "VIP")
pos_df <- tibble(Wavenumber = important_wavenumbers_coeff_positive, Type = "Positive")
neg_df <- tibble(Wavenumber = important_wavenumbers_coeff_negative, Type = "Negative")
bands_df <- bind_rows(vip_df, pos_df, neg_df)
bands_df <- bands_df %>%
  mutate(Vjust = ifelse(Wavenumber %in% label_right_bands, 1.1, -0.2))

plot_spectra_with_bands(mXr00, mXr05, bands_df)
ggsave("data/out/Bursera_simaruba/Results_BS_model_important_bands.png") # or save manually via >Export

#-------------------------------------------------------------------------------------------
### Maclura tinctoria
# define important bands
important_wavenumbers_vip <- c(1284, 831, 1145, 1767, 1124, 1313)
important_wavenumbers_coeff_positive <- c(948, 831, 972, 1460)
important_wavenumbers_coeff_negative <- c(807, 1767, 1231, 1275, 1321, 1140)
label_right_bands <- c(1275, 1124, 1313, 1140)

# df creation
vip_df <- tibble(Wavenumber = important_wavenumbers_vip, Type = "VIP")
pos_df <- tibble(Wavenumber = important_wavenumbers_coeff_positive, Type = "Positive")
neg_df <- tibble(Wavenumber = important_wavenumbers_coeff_negative, Type = "Negative")
bands_df <- bind_rows(vip_df, pos_df, neg_df)
bands_df <- bands_df %>%
  mutate(Vjust = ifelse(Wavenumber %in% label_right_bands, 1.1, -0.2))

plot_spectra_with_bands(mXr00, mXr09, bands_df)
ggsave("data/out/Maclura_tinctoria/Results_MT_model_important_bands.png") # or save manually via >Export

#-------------------------------------------------------------------------------------------
### Pinus patula
# define important bands
important_wavenumbers_vip <- c(1122, 1663, 1788, 938, 1254)
important_wavenumbers_coeff_positive <- c(1167, 865, 1046, 1551)
important_wavenumbers_coeff_negative <- c(938, 924, 830, 1114, 1254, 1145)
label_right_bands <- c(1114, 924)

# df creation
vip_df <- tibble(Wavenumber = important_wavenumbers_vip, Type = "VIP")
pos_df <- tibble(Wavenumber = important_wavenumbers_coeff_positive, Type = "Positive")
neg_df <- tibble(Wavenumber = important_wavenumbers_coeff_negative, Type = "Negative")
bands_df <- bind_rows(vip_df, pos_df, neg_df)
bands_df <- bands_df %>%
  mutate(Vjust = ifelse(Wavenumber %in% label_right_bands, 1.1, -0.2))

plot_spectra_with_bands(mXr00, mXr21, bands_df)
ggsave("data/out/Pinus_patula/Results_PP_model_important_bands.png") # or save manually via >Export