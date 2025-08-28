# set wd to root dir
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # assuming use of RStudio
parent_dir <- dirname(current_dir)
setwd(parent_dir) # or set wd manually here

library(ggplot2)
library(reshape2)
library(dplyr)
library(prospectr)
source("Utils.R")
#-------------------------------------------------------------------------------

# paths to spectra per species
species_paths <- c(
  "data/FTIR_Spectra/Bursera_Simaruba/",
  "data/FTIR_Spectra/Maclura_Tinctoria/",
  "data/FTIR_Spectra/Pinus_Patula/"
)

colors <- c("blue", "green", "#d27570")

# list to store mean spectra
mean_spectra_list <- list()
wn_reference <- NULL
ref_band <- 1026 # to normalize

# for each species
for (i in seq_along(species_paths)) {
  path <- species_paths[i]
  
  spec_wn <- read_spec(path)
  data_mat <- spec_wn$spec
  wn <- spec_wn$wn
  
  # Slice to fingerprint region: 1800–800 cm⁻¹
  fingerprint_cols <- wn >= 800 & wn <= 1800
  wn_fingerprint <- wn[fingerprint_cols]
  spec_fingerprint <- data_mat[, fingerprint_cols, drop = FALSE]
  
  if (is.null(wn_reference)) {
    wn_reference <- wn_fingerprint
  }
  
  # baseline correction
  spec_baseline_corrected <- baseline(spec_fingerprint)
  
  # ix of ref band (e.g. 1026 cm⁻¹)
  band_idx <- which.min(abs(wn_fingerprint - ref_band))
  spec_normalized <- t(apply(spec_baseline_corrected, 1, function(sp) sp / sp[band_idx]))

  # calc mean spec
  mean_spectrum <- colMeans(spec_normalized)
  mean_spectra_list[[i]] <- mean_spectrum
}

#-------------------------------------------------------------------------------
# df creation
df_plot <- data.frame(
  wavenumber = wn_reference,
  Bursera = mean_spectra_list[[1]],
  Maclura = mean_spectra_list[[2]],
  Pinus = mean_spectra_list[[3]]
)

df_long <- melt(df_plot, id.vars = "wavenumber", variable.name = "Species", value.name = "Absorbance")

species_colors <- c("Bursera" = "blue", "Maclura" = "forestgreen", "Pinus" = "#d27570") # darkred

p <- ggplot(df_long, aes(x = wavenumber, y = Absorbance, color = Species)) +
  geom_line(linewidth = 1.2) +
  scale_x_reverse(breaks = seq(1800, 800, by = -200)) +  
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = species_colors) +
  labs(x = "Wavenumber (cm⁻¹)", y = "Absorbance (a.u.)") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),  # <-- remove black bottom/left lines
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1.3),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#-------------------------------------------------------------------------------
# add vertical dashed lines based on literature search
# this is done manually to be able to make spatial adjustments when text and spectra are overlapping
p <- p + 
  geom_text(data = NULL, aes(x = 1770, y = 0.06, label = "1770"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1770, xend = 1770, y = 0, yend = 0.04), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1731, y = 0.145, label = "1731"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1731, xend = 1731, y = 0, yend = 0.125), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1690, y = 0.14, label = "1690"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1690, xend = 1690, y = 0, yend = 0.12), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1660, y = 0.255, label = "1660"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1660, xend = 1660, y = 0, yend = 0.235), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1639, y = 0.4, label = "1639"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1639, xend = 1639, y = 0, yend = 0.38), linetype = "dashed", color = "gray33", linewidth=0.2) +

  geom_text(data = NULL, aes(x = 1618, y = 0.365, label = "1618"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1618, xend = 1618, y = 0, yend = 0.355), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1593, y = 0.4, label = "1593"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1593, xend = 1593, y = 0, yend = 0.38), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1510, y = 0.18, label = "1510"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1510, xend = 1510, y = 0, yend = 0.16), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1485, y = 0.235, label = "1465"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1465, xend = 1465, y = 0, yend = 0.235), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1457, y = 0.295, label = "1457"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1457, xend = 1457, y = 0, yend = 0.275), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1448, y = 0.265, label = "1452"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1452, xend = 1452, y = 0, yend = 0.245), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1425, y = 0.235, label = "1425"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1425, xend = 1425, y = 0, yend = 0.22), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1375, y = 0.22, label = "1375"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1375, xend = 1375, y = 0, yend = 0.21), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1334, y = 0.31, label = "1334"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1334, xend = 1334, y = 0, yend = 0.29), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1323, y = 0.275, label = "1328"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1328, xend = 1328, y = 0, yend = 0.255), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1310, y = 0.245, label = "1315"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1320, xend = 1320, y = 0, yend = 0.235), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1264, y = 0.29, label = "1264"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1264, xend = 1264, y = 0, yend = 0.27), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1244, y = 0.33, label = "1244"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1244, xend = 1244, y = 0, yend = 0.31), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1224, y = 0.29, label = "1224"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1224, xend = 1224, y = 0, yend = 0.27), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1157, y = 0.305, label = "1157"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1157, xend = 1157, y = 0, yend = 0.295), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1145, y = 0.35, label = "1140"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1140, xend = 1140, y = 0, yend = 0.33), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1106, y = 0.5, label = "1106"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1106, xend = 1106, y = 0, yend = 0.5), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1084, y = 0.65, label = "1084"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1084, xend = 1084, y = 0, yend = 0.65), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1060, y = 0.94, label = "1048"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1048, xend = 1048, y = 0, yend = 0.92), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1055, y = 1, label = "1041"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1041, xend = 1041, y = 0, yend = 1.0), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 1002, y = 1.0, label = "1007"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 1007, xend = 1007, y = 0, yend = 1.0), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 984, y = 0.8, label = "984"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 984, xend = 984, y = 0, yend = 0.78), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 963, y = 0.62, label = "968"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 968, xend = 968, y = 0, yend = 0.6), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 925, y = 0.15, label = "921"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 921, xend = 921, y = 0, yend = 0.13), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 907, y = 0.125, label = "907"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 907, xend = 907, y = 0, yend = 0.115), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 898, y = 0.18, label = "898"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 898, xend = 898, y = 0, yend = 0.16), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 860, y = 0.06, label = "860"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 860, xend = 860, y = 0, yend = 0.05), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 832, y = 0.08, label = "832"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 832, xend = 832, y = 0, yend = 0.06), linetype = "dashed", color = "gray33", linewidth=0.2) +
  
  geom_text(data = NULL, aes(x = 805, y = 0.06, label = "805"), size = 3.5, color = "black") +
  geom_segment(data = NULL, aes(x = 805, xend = 805, y = 0, yend = 0.05), linetype = "dashed", color = "gray33", linewidth=0.2)

print(p)

ggsave("Results_average_spectra.png", width = 2800 / 300, height = 1681 / 300, dpi = 300)