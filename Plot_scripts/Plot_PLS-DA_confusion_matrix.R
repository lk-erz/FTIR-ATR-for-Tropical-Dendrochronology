# set wd to root dir
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # assuming use of RStudio
parent_dir <- dirname(current_dir)
setwd(parent_dir) # or set wd manually here

library(caret)
library(ggplot2)
library(patchwork)
#-------------------------------------------------------------------------------
# relabel function
relab <- function(x) factor(ifelse(x == 1, "LW", "EW"), levels = c("LW", "EW"))

# Create confusion matrices
cm1 <- confusionMatrix(relab(c(1,0,1,1,0)), relab(c(1,0,1,0,0)))
cm2 <- confusionMatrix(relab(c(1,0,0,1,1)), relab(c(1,1,0,0,1)))
cm3 <- confusionMatrix(relab(c(0,1,1,1,0)), relab(c(0,0,1,1,0)))
cm4 <- confusionMatrix(relab(c(1,1,0,0,1)), relab(c(1,0,0,1,1)))

# table has to be modified manually according to classification output!
cm1$table["LW", "LW"] <- 57; cm1$table["LW", "EW"] <- 22
cm1$table["EW", "LW"] <- 10; cm1$table["EW", "EW"] <- 69

cm2$table["LW", "LW"] <- 70; cm2$table["LW", "EW"] <- 9
cm2$table["EW", "LW"] <- 6; cm2$table["EW", "EW"] <- 73

cm3$table["LW", "LW"] <- 74; cm3$table["LW", "EW"] <- 5
cm3$table["EW", "LW"] <- 12; cm3$table["EW", "EW"] <- 67

cm4$table["LW", "LW"] <- 68; cm4$table["LW", "EW"] <- 11
cm4$table["EW", "LW"] <- 14; cm4$table["EW", "EW"] <- 65

max_freq <- max(cm1$table, cm2$table, cm3$table, cm4$table) # for color scale
# plot
plot_cm <- function(cm, title, max_freq, show_legend = TRUE, show_y_axis = TRUE) {
  df <- as.data.frame(as.table(cm$table))
  
  border_rect <- data.frame(xmin = 0.5, xmax = 2.5, ymin = 0.5, ymax = 2.5)
  
  p <- ggplot(df, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq)) +
    geom_text(aes(label = Freq), size = 5) +
    geom_rect(data = border_rect,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              color = "black", fill = NA, linewidth = 0.6) +
    scale_fill_gradient(
      low = "white", 
      high = "steelblue", 
      limits = c(0, 75),
      breaks = seq(0, 75, by = 15)
    ) +
    scale_y_discrete(limits = rev(levels(df$Prediction))) +
    coord_fixed() +
    labs(x = "True class", y = "Predicted class") +
    ggtitle(title) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.ticks.length = unit(0, "cm"),
      axis.ticks.x.top = element_blank(),
      axis.text.x.top = element_blank(),
      axis.title.x.top = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(color = "black", size = 12)
    ) +
    guides(fill = guide_colorbar(
      title = "Count", 
      title.position = "top", 
      label.theme = element_text(size = 14, face = "plain"),
      barwidth = 2,  # Make the color bar wider
      barheight = 10
    ))
  
  if (!show_y_axis) {
    p <- p + theme(
      axis.title.y.left = element_blank(),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank()
    )
  }
  
  return(p)
}

# subplots
p1 <- plot_cm(cm1, "Raw data", max_freq, FALSE, TRUE)
p2 <- plot_cm(cm2, "2nd derivative + \nSNV", max_freq, FALSE, FALSE)
p3 <- plot_cm(cm3, "Baseband + \nnorm 1424", max_freq, FALSE, FALSE)
p4 <- plot_cm(cm4, "1st derivative + \nvector norm", max_freq, TRUE, FALSE)

# combine
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
print(combined_plot)
