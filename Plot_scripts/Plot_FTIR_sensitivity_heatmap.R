# set wd to root dir
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # assuming use of RStudio
parent_dir <- dirname(current_dir)
setwd(parent_dir) # or set wd manually here

library(ggplot2)
library(reshape2)
#-------------------------------------------------------------------------------
# sensitivity analysis
n_scan <- c(10, 20, 40, 60, 80, 100, 120)
spec_res <- c(1, 2, 4, 8, 16)

time_mat <- matrix(
  nrow = length(spec_res),
  ncol = length(n_scan),
  dimnames = list(
      spec_res,
      n_scan)
  )

time_mat[,] <- c( # column-wise filling!
  28, 17, 12, 9, 8,
  56, 35, 24, 19, 16,
  113, 70, 48, 37, 32,
  169, 105, 72, 56, 48,
  226, 139, 96, 75, 64,
  282, 174, 120, 93, 80,
  339, 209, 144, 112, 96
)

time_long <- melt(time_mat)
colnames(time_long) <- c("spec_res", "n_scan", "Time")
time_long$spec_res <- factor(time_long$spec_res, levels = rev(rownames(time_mat)))
time_long$n_scan <- factor(time_long$n_scan, levels = colnames(time_mat))

ggplot(time_long, aes(x = n_scan, y = spec_res, fill = Time)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Time), color = "black", size = 5) +
  scale_fill_gradient(low = "yellow", high = "red", name = "Time [s]") +
  labs(
    x = "n of averaged scans",
    y = "Spectral resolution (cm-1)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1.5, "cm"),
    panel.grid = element_blank()
  )
  