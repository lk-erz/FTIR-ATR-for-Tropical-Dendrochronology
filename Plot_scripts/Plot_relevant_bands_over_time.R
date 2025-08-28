# set wd to root dir
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # assuming use of RStudio
parent_dir <- dirname(current_dir)
setwd(parent_dir) # or set wd manually here

library(ggplot2)
library(prospectr)
library(dplyr)
library(patchwork)
library(tidyr)
source("Utils.R")
#-------------------------------------------------------------------------------

path = "data/FTIR_Spectra/Bands_over_time/"
species = "Bursera_simaruba"
filename_meta <- paste0(species, "_file_meta_valid_years.csv")
filename_dat <- paste0(species, "_data_mat_valid_years.csv")

path_meta <- file.path(path, filename_meta)
path_spectra <- file.path(path, filename_dat)

# read saved data
file_meta <- read.csv(path_meta)
data <- read.csv(path_spectra)
colnames(data) <- gsub("^X", "", colnames(data))
colnames(data) <- as.numeric(colnames(data))
wn <- as.numeric(colnames(data))

#-------------------------------------------------------------------------------
# slice data to fingerprint region 1800 - 800 cm-1
fingerprint_cols <- wn >= 800 & wn <= 1800
wn4000 <- wn
wn <- wn[fingerprint_cols]
spec_fingerprint <- data[, fingerprint_cols]

### preprocessing ###
mXr00 <- as.matrix(spec_fingerprint)
mXr02 <- standardNormalVariate(mXr00)

# BASELINE + SNV
mXrBase <- baseline(data)
colnames(mXrBase) <- wn4000
mXrBaseBand <- manual_baseline_correction(data, wn=as.numeric(wn4000), anchors = c(3750, 1780, 1530, 1185, 820), window = 50) # for Pinus p. 1530, not 1485
mXrBase <- mXrBase[, fingerprint_cols]
mXrBaseBand <- mXrBaseBand[, fingerprint_cols]

#-------------------------------------------------------------------------------
### absorbance mean by growth_years
combined_data <- data.frame(growth_years = file_meta$growth_years, mXr02, check.names = FALSE)

# calc mean + SD for each group
stats_spectra <- combined_data %>%
  group_by(growth_years) %>%
  summarise(across(everything(), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))  # Replace NA with 0

avg_mXr02 <- as.matrix(stats_spectra[, grep("_mean$", names(stats_spectra))])
sd_mXr02 <- as.matrix(stats_spectra[, grep("_sd$", names(stats_spectra))])

# assign column names based on wavenumbers
colnames(avg_mXr02) <- as.character(colnames(mXr02))
colnames(sd_mXr02) <- as.character(colnames(mXr02))

avg_growth_years <- stats_spectra$growth_years
#-----------------------------------------------------------------------------------------------
### for lignin-to-polysaccharide ratios, use baseline-corrected spectra
bands_ratio <- c("1508.0612", "1374.9969", "1159.0087") # bands used for ratios
mXrBaseBand_ratio <- mXrBaseBand[, bands_ratio, drop = FALSE] # use baseline-corrected spec
combined_data_ratio <- data.frame(growth_years = file_meta$growth_years, mXrBaseBand_ratio, check.names = FALSE)

stats_ratio <- combined_data_ratio %>%
  group_by(growth_years) %>%
  summarise(across(everything(), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

avg_ratio <- as.matrix(stats_ratio[, grep("_mean$", names(stats_ratio))])
sd_ratio <- as.matrix(stats_ratio[, grep("_sd$", names(stats_ratio))])
colnames(avg_ratio) <- bands_ratio
colnames(sd_ratio) <- bands_ratio

#   Lignin:      1593, 1508
#   cellulose:   1375, 893; alternatives: 1057, 1280, 1157, 1315 (Popescu et al. 2007), 1180)
#   Ratio:       1375, 1159
# Popescu et al. (2007): important cellulose bands: 1315 (frequently used by models), 1280 (used by Maclura models), 1180, 1060 (used by models)

#####################
# FINAL BANDS
#####################

if (species == "Bursera_simaruba") {
  bands_model <- c("1420.3159", "1327.7495") # Bursera (insert two important bands of the respective model here)
  band_assignment_1 = "1420 cm-1: lignin/cellulose"
  band_assignment_2 = "1327 cm-1: lignin"
} else if (species == "Maclura_tinctoria") {
  bands_model <- c("1767.4399", "1280.5021") # Maclura (insert two important bands of the respective model here)
  band_assignment_1 = "1767 cm-1: ester/hemicellulose"
  band_assignment_2 = "1280 cm-1: cellulose/lignin"
} else if (species == "Pinus_patula") {
  bands_model <- c("1121.4036", "938.1992") # Pinus (insert two important bands of the respective model here)
  band_assignment_1 = "1121 cm-1: lignin/cellulose"
  band_assignment_2 = "938 cm-1: lignin"
}

# df creation
df <- data.frame(
  avg_growth_years,
  `1593 cm-1: lignin` = avg_mXr02[, "1592.9137"],
  `1508 cm-1: lignin` = avg_mXr02[, "1508.0612"],
  
  `1375 cm-1: (hemi)cellulose` = avg_mXr02[, "1374.9969"],
  `893 cm-1: (hemi)cellulose` = avg_mXr02[, "892.8802"],

  setNames(list(avg_mXr02[, bands_model[1]]), band_assignment_1),
  setNames(list(avg_mXr02[, bands_model[2]]), band_assignment_2),

  `ratio 1508 / 1375 cm-1` = (avg_ratio[, "1508.0612"] / avg_ratio[, "1374.9969"]),
  `ratio 1508 / 1159 cm-1` = (avg_ratio[, "1508.0612"] / avg_ratio[, "1159.0087"]),
  
  check.names=F
) %>%
  tidyr::pivot_longer(-avg_growth_years, names_to = "Wavenumber", values_to = "Absorbance")

# df of SDs
df_sd <- data.frame(
  avg_growth_years,
  `1593 cm-1: lignin` = sd_mXr02[, "1592.9137"],
  `1508 cm-1: lignin` = sd_mXr02[, "1508.0612"],
  
  `1375 cm-1: (hemi)cellulose` = sd_mXr02[, "1374.9969"], # 1738.5129, 1374.9969, 1158.0444
  `893 cm-1: (hemi)cellulose` = sd_mXr02[, "892.8802"],
  
  # alternatively 1420.3159 or 1455.9926 1785.7604 1228.4335, 1048.1218, 1083.7985
  setNames(list(sd_mXr02[, bands_model[1]]), band_assignment_1),
  setNames(list(sd_mXr02[, bands_model[2]]), band_assignment_2),
  
  `ratio 1508 / 1375 cm-1` = (sd_ratio[, "1508.0612"] / sd_ratio[, "1374.9969"]),
  `ratio 1508 / 1159 cm-1` = (sd_ratio[, "1508.0612"] / sd_ratio[, "1159.0087"]),
  
  check.names=F
) %>%
  tidyr::pivot_longer(-avg_growth_years, names_to = "Wavenumber", values_to = "SD")

# preserve order
wavenumber_order <- c(
  "1593 cm-1: lignin",
  "1508 cm-1: lignin",
  "1375 cm-1: (hemi)cellulose",
  "893 cm-1: (hemi)cellulose",
  #"1767 cm-1: ester/hemicellulose",
  #"1280 cm-1: cellulose/lignin",
  #"1420 cm-1: lignin/cellulose",
  #"1327 cm-1: lignin",
  band_assignment_1,
  band_assignment_2,
  "ratio 1508 / 1375 cm-1",
  "ratio 1508 / 1159 cm-1"
)
df <- df %>%
  mutate(Wavenumber = factor(Wavenumber, levels = wavenumber_order))

df_sd <- df_sd %>%
  mutate(Wavenumber = factor(Wavenumber, levels = wavenumber_order))

# subplot annotations
tag_labels <- letters[1:8]  # "a" to "h"
names(tag_labels) <- wavenumber_order
df$tag <- tag_labels[as.character(df$Wavenumber)]

#--------------------------
# merge mean and SD data
df <- left_join(df, df_sd, by = c("avg_growth_years", "Wavenumber"))

ylabel <- function(label1, label2, label3, label4) {
  L1 <- nchar(label1)
  L2 <- nchar(label2)
  L3 <- nchar(label3)
  L4 <- nchar(label4)
  
  total_length <- L1 + L2 + L3 + L4
  scaler <- ifelse(total_length > 20, 4, 0)
  
  space1 <- paste0(rep(" ", 18 - floor(L1 / 2)), collapse = "")
  space2 <- paste0(rep(" ", 26 - floor((L1 + L2) / 2) - scaler), collapse = "")
  space3 <- paste0(rep(" ", 26 - floor((L2 + L3) / 2) - scaler), collapse = "")
  space4 <- paste0(rep(" ", 26 - floor((L3 + L4) / 2) - scaler), collapse = "")
  space5 <- paste0(rep(" ", 18 - floor(L4 / 2)), collapse = "")
  
  paste0(space1, label1, space2, label2, space3, label3, space4, label4, space5)
}
ylabs <- ylabel("Absorbance (a.u.)", "Absorbance (a.u.)", "Absorbance (a.u.)", "Absorbance (a.u.)")

# ticks for different species
if (species == "Bursera_simaruba") {
  x_breaks <- seq(0, 100, by = 5)
  x_labels <- ifelse(x_breaks %% 10 == 0, x_breaks, "")
} else if (species == "Maclura_tinctoria") {
  x_breaks <- seq(0, 65, by = 5)
  x_labels <- ifelse(x_breaks %% 5 == 0, x_breaks, "")
} else if (species == "Pinus_patula") {
  x_breaks <- seq(0, 20, by = 2)
  x_labels <- ifelse(x_breaks %% 2 == 0, x_breaks, "")
}

# plot
p_main <- ggplot(df, aes(x = avg_growth_years, y = Absorbance)) +
  geom_smooth(method = "loess", color = "red", linetype = "dashed", se = T, size = 0.7) +  # trend line
  geom_point(color = "black", size = 1.5) +
  #geom_errorbar(aes(ymin = Absorbance - SD, ymax = Absorbance + SD), width = 1, color = "black") +  ### rm comment for error bars with SDs
  geom_text(data = df %>% group_by(Wavenumber) %>% slice(1), 
            aes(x = 0, y = Inf, label = tag), 
            hjust = 1.8, vjust = 1.0, size = 5, inherit.aes = FALSE) +
  facet_wrap(~Wavenumber, ncol = 2, scales = "free") +
  scale_x_continuous(breaks = x_breaks, labels = rep("", length(x_breaks))) +
  ylab(ylabs) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_blank(),        # remove x-axis labels
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length = unit(8, "pt"),  # increase length, default is ~5
    axis.ticks.x = element_line(color = "black"),
    axis.line = element_line()
  )

# dummy plots to add x-axis labels below each column
p_x1 <- ggplot(data.frame(x = x_breaks, y = 0), aes(x, y)) +
  geom_blank() +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0, max(x_breaks))) +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 16, hjust = 0.5, margin = margin(t = -20)),
    axis.text.x = element_text(size = 16, color = "black", margin = margin(t = -48, r = 25)),
  ) +
  labs(x = "Distance to pith (years)")

p_x2 <- p_x1 # identical layout for the second column

# combine using patchwork
final_plot <- (p_main) / (p_x1 + p_x2) + plot_layout(heights = c(20, 2))
print(final_plot)

ggsave(file.path("data/out", paste0(species, "/Results_all_bands_over_time.png")), width=12, height=8)
