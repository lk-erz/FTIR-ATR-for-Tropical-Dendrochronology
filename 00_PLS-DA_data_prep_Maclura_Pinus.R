#==============================================================================
# Supplementary R scripts for mater's thesis
# Title: Potential of FTIR spectroscopy for tropical dendrochronology
# Author: Lukas Erzfeld
# Insitution: University of Leipzig
# Date: August 26, 2025
#===============================================================================

# required packages may have to be installed primarily
# install.packages(c("dplyr", "dplR", "writexl"))
library(dplyr)
library(dplR)
library(writexl)
source("Utils.R")
source("Utils_PLS-DA_data_prep.R")
#-------------------------------------------------------------------------------
preprocess_maclura_pinus <- function(path_data, path_out, first_year_dict, model, species_name, heartwood_years=NA) {
  cat("Preparing data...\n")
  
  if (model == "PLS-DA") {
    
    # read data and meta
    spec_wn <- read_spec(path_data)
    data_mat <- spec_wn$spec
    wn <- spec_wn$wn
    file_meta <- read_file_metadat(path_data)
    #------------
    # remove (slightly) mixed spectra
    result <- remove_mixed_spectra(file_meta, data_mat)
    file_meta <- result$metadata
    data_mat  <- result$spectra_matrix
    #------------
    ### calc age
    file_meta <- calc_age_from_dict(file_meta, first_year_dict)
    #------------
    # filter heartwood
    heartwood <- heartwood_years # choose how many years are removed as heartwood
    result <- filter_heartwood(file_meta, data_mat, heartwood) 
    file_meta <- result$metadata
    data_mat  <- result$spectra_matrix
    #------------
    # until here reasonable removal of spectra -> still remaining more EW than LW -> random sampling
    result <- balance_EW_LW(file_meta, data_mat)
    file_meta <- result$metafile
    data_mat <- result$spectra_matrix
    #------------
    ### reorder data according to earlywood/latewood for better visualization of classification
    ix_sort <- with(file_meta, order(ring, loc, sample, ind, year))
    file_meta <- file_meta[ix_sort, ]
    data_mat  <- data_mat[ix_sort, ]
  #-----------
    # write to file
    cat("Writing data...\n")
    write_sampled_data(file_meta, data_mat, path_out, species_name, method = model)
    cat("Data preparation for", species_name, "data successful. \nReady for PLS-DA ring classification.")
    
  } else if (model == "PLSR") {
    
    # read data and meta
    spec_wn <- read_spec(path_data)
    data_mat <- spec_wn$spec
    wn <- spec_wn$wn
    file_meta <- read_file_metadat(path_data)
    #------------
    # remove (slightly) mixed spectra
    #result <- remove_mixed_spectra(file_meta, data_mat)
    #file_meta <- result$metadata
    #data_mat  <- result$spectra_matrix
    #------------
    ### calc age
    file_meta <- calc_age_from_dict(file_meta, first_year_dict)
    #------------
    # filter heartwood
    heartwood <- heartwood_years # choose how many years are removed as heartwood
    result <- filter_heartwood(file_meta, data_mat, heartwood) 
    file_meta <- result$metadata
    data_mat  <- result$spectra_matrix
    #------------
    ### reorder data according to earlywood/latewood for better visualization of classification
    #ix_sort <- with(file_meta, order(ring, loc, sample, ind, year))
    #file_meta <- file_meta[ix_sort, ]
    #data_mat  <- data_mat[ix_sort, ]
    #-----------
    
    write_sampled_data(file_meta, data_mat, path_out, species_name, method = model)
    cat("Data preparation for", species_name, "data successful. \nReady for PLSR age prediction.")
  } else {
    stop("Must choose PLS-DA or PLSR as model parameter.")
  }
}

preprocess_maclura_pinus_valid_years <- function(path_data, path_out, first_year_dict, species_name) {
  cat("Preparing data...\n")
  # species_name, heartwood_years=NA
  # read data and meta
  spec_wn <- read_spec(path_data)
  data_mat <- spec_wn$spec
  wn <- spec_wn$wn
  file_meta <- read_file_metadat(path_data)
  #------------
  #metadata$comment[is.na(metadata$comment)] <- "valid"
  
  #------------
  ### calc age
  file_meta <- calc_age_from_dict(file_meta, first_year_dict)
  
  #------------
  ### reorder data according to earlywood/latewood for better visualization of classification
  #ix_sort <- with(file_meta, order(ring, loc, sample, ind, year))
  #file_meta <- file_meta[ix_sort, ]
  #data_mat  <- data_mat[ix_sort, ]
  #-----------
  # write to file
  cat("Writing data...\n")
  write_sampled_data(file_meta, data_mat, path_out, species_name, method = "valid_years")
  cat("Data preparation for", species_name, "data successful. \nReady for plotting band values over time.")
}
