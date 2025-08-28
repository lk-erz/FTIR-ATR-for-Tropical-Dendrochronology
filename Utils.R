#==============================================================================
# Supplementary R scripts for mater's thesis
# Title: Potential of FTIR spectroscopy for tropical dendrochronology
# Author: Lukas Erzfeld
# Insitution: University of Leipzig
# Date: August 26, 2025
#===============================================================================

read_spec <- function(path) {

  # set current directory to working directory
  wd = dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(wd)
  
  files <- list.files(path, full.names = TRUE, recursive = TRUE) # recursive = also subdirs
  files <- files[!file.info(files)$isdir]
  
  data_list <- lapply(files, function(file) {
    df <- read.csv2(file, 
                    header = FALSE,
                    nrows = 3736,
                    skip = 19,
                    dec = ",",
                    col.names = c("wn", "absorbance"))
    df$absorbance
  })
  
  # use wn values as colnames
  wn <- read.csv2(files[1], 
                  header = FALSE,
                  nrows = 3736,
                  skip = 19,
                  dec = ",",
                  col.names = c("wn", "absorbance"))$wn
  
  data_mat <- do.call(rbind, data_list)
  colnames(data_mat) <- wn
  
  return(list(spec = data_mat,wn = wn))

}


read_file_metadat <- function(path) {
  # list all files including in subdirs
  files <- list.files(path, full.names = TRUE, recursive = TRUE)
  files <- files[!file.info(files)$isdir]
  
  # extract information from file names
  file_names <- basename(files)
  
  # split file names by "_"
  file_info_list <- strsplit(file_names, "_")
  
  # create a data frame with extracted information
  file_info_df <- data.frame(
    species = sapply(file_info_list, function(x) x[1]),
    sample = sapply(file_info_list, function(x) x[2]),
    ind = sapply(file_info_list, function(x) x[3]),
    loc = sapply(file_info_list, function(x) x[4]),
    ring = sapply(file_info_list, function(x) x[5]),
    year = sapply(file_info_list, function(x) gsub("\\..*$", "", x[6])),
    comment = sapply(file_info_list, function(x) {
      if (length(x) >= 7) {
        gsub("\\..*$", "", x[7]) # rm file extension if still present
      } else {
        NA
      }
    }),
    stringsAsFactors = FALSE
  )
  return(file_info_df)
}

manual_baseline_correction <- function(spectra, wn, anchors, window = 30) {
  # spectra: matrix (rows = samples, columns = intensities at wn)
  # wn: numeric vector of wavenumbers (length must match ncol(spectra))
  # anchors: vector of wavenumbers (cm⁻¹) to define baseline region centers
  # window: half-width of range (in cm⁻¹) to search for local minima
  
  corrected <- t(apply(spectra, 1, function(spectrum) {
    # find local minima near each anchor point
    minima_indices <- sapply(anchors, function(anchor) {
      range_idx <- which(wn >= (anchor - window) & wn <= (anchor + window))
      local_min_idx <- range_idx[which.min(spectrum[range_idx])]
      return(local_min_idx)
    })
    
    baseline <- approx(x = wn[minima_indices], 
                       y = spectrum[minima_indices], 
                       xout = wn, method = "linear", rule = 2)$y
    spectrum - baseline
  }))
  
  return(corrected)
}


normalize_to_band <- function(spectra, wn, band) {
  band_idx <- which.min(abs(wn - band))
  t(apply(spectra, 1, function(sp) sp / sp[band_idx]))
}

vector_normalize <- function(x) {
  x / sqrt(sum(x^2))
}
