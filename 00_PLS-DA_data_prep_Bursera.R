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
preprocess_bursera_PLSDA <- function(path_data, path_out, model, heartwood_years=NA) {
  #------------
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
  # divide file_meta by locations Umpalá and La Fuente and
  # calc age of tree of respective spectrum if a year is provided
  # (division by location is necessary because samples have same labels)
  
  # keep spectra without year information aside
  noYear_idx <- which(file_meta$year == "noYear")
  file_meta_noYear <- file_meta[noYear_idx, ]
  data_mat_noYear   <- data_mat[noYear_idx, ]
  
  # rm from data set to calc age
  file_meta <- file_meta[-noYear_idx, ]
  data_mat  <- data_mat[-noYear_idx, ]
  
  ### calc age
  # Umpalá
  ump_ix <- which(file_meta$loc == "Ump")
  file_meta_ump <- file_meta[ump_ix, ]
  data_mat_ump <- data_mat[ump_ix, ]
  
  rw <- read.rwl("data/Bursera_simaruba_ring_widths/Bursera_Umpala.rwl")
  # for each tree get year of growth start
  first_years <- apply(rw, 2, function(col) {
    years <- rownames(rw)
    first_non_na <- which(!is.na(col))[1]  # index of first non-NA
    if (!is.na(first_non_na)) {
      return(years[first_non_na])
    } else {
      return(NA)  # if the whole column is NA
    }
  })
  
  names(first_years) <- substring(names(first_years), 4)
  first_years
  # create col in meta data to match chronology data -> join
  file_meta_ump$sample_ind <- paste(file_meta_ump$sample, file_meta_ump$ind, sep = "_")
  file_meta_ump$first_year <- first_years[file_meta_ump$sample_ind]
  # create new col containing years after growth of tree
  file_meta_ump$year <- as.numeric(file_meta_ump$year)
  file_meta_ump$first_year <- as.numeric(file_meta_ump$first_year)
  file_meta_ump$growth_years <- file_meta_ump$year - file_meta_ump$first_year
  
  # filter heartwood
  heartwood <- heartwood_years # choose how many years are removed as heartwood
  result <- filter_heartwood(file_meta, data_mat, heartwood) 
  file_meta <- result$metadata
  data_mat  <- result$spectra_matrix
  
  # La Fuente
  laf_ix <- which(file_meta$loc == "LaF")
  file_meta_laf <- file_meta[laf_ix, ]
  data_mat_laf <- data_mat[laf_ix, ]
  
  rw <- read.rwl("data/Bursera_simaruba_ring_widths/Bursera_La_Fuente.rwl")
  # for each tree get year of growth start
  first_years <- apply(rw, 2, function(col) {
    years <- rownames(rw)
    first_non_na <- which(!is.na(col))[1]  # ix of first non-NA
    if (!is.na(first_non_na)) {
      return(years[first_non_na])
    } else {
      return(NA)  # if the whole column is NA
    }
  })
  
  names(first_years) <- substring(names(first_years), 4)
  first_years
  # create col in meta data to match chronology data -> join
  file_meta_laf$sample_ind <- paste0("A", file_meta_laf$sample, "_", file_meta_laf$ind)#, "_") # according to .rwl
  file_meta_laf$first_year <- first_years[file_meta_laf$sample_ind]
  # create new col containing years after growth of tree
  file_meta_laf$year <- as.numeric(file_meta_laf$year)
  file_meta_laf$first_year <- as.numeric(file_meta_laf$first_year)
  file_meta_laf$growth_years <- file_meta_laf$year - file_meta_laf$first_year
  
  # filter heartwood
  heartwood <- heartwood_years # choose how many years are removed as heartwood
  result <- filter_heartwood(file_meta, data_mat, heartwood) 
  file_meta <- result$metadata
  data_mat  <- result$spectra_matrix
  
  #------
  # recombine
  file_meta_ump$year <- as.character(file_meta_ump$year)
  file_meta_laf$year <- as.character(file_meta_laf$year)
  
  # use bind_rows to automatically fill missing columns with NA
  file_meta <- bind_rows(file_meta_ump, file_meta_laf, file_meta_noYear)
  data_mat <- bind_rows(as.data.frame(data_mat_ump), as.data.frame(data_mat_laf), as.data.frame(data_mat_noYear))
  
  # check dimensions
  print(nrow(file_meta_noYear) == nrow(data_mat_noYear))  # should be TRUE
  print(nrow(file_meta_ump)    == nrow(data_mat_ump))
  print(nrow(file_meta_laf)    == nrow(data_mat_laf))
  #------------
  # until here reasonable removal of spectra -> remaining 371 EW vs. 267 LW -> random sampling to balance dataset
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
  write_sampled_data(file_meta, data_mat, "data/FTIR_Spectra/PLS-DA", "Bursera_simaruba", method = "PLS-DA")
  
  cat("Data preparation for Bursera simaruba data successful. \nReady for PLS-DA ring classification.")
}


preprocess_bursera_PLSR <- function(path_data, path_out) {
  #------------
  # read data and meta
  spec_wn <- read_spec(path_data)
  data_mat <- spec_wn$spec
  wn <- spec_wn$wn
  file_meta <- read_file_metadat(path_data)
  #------------
  # remove not suitable spectra
  bad_comments <- c("weakSignal", "plusminus2jahre", "weirdShape 1500", "noAbsOver03")
  rows_to_keep <- !file_meta$comment %in% bad_comments
  file_meta <- file_meta[rows_to_keep, ]
  data_mat <- data_mat[rows_to_keep, ]
  
  #------------
  # divide file_meta by locations Umpalá and La Fuente and
  # calc age of tree of respective spectrum if a year is provided
  # (division by location is necessary because samples have same labels)
  
  # filter meta data for spectra that contain valid year information
  year_ix <- which(file_meta$year != "noYear")
  file_meta <- file_meta[year_ix, ]
  data_mat <- data_mat[year_ix, ]
  
  ### calc age
  # Umpalá
  ump_ix <- which(file_meta$loc == "Ump")
  file_meta_ump <- file_meta[ump_ix, ]
  data_mat_ump <- data_mat[ump_ix, ]
  
  rw <- read.rwl("Bursera_Umpala.rwl")
  # for each tree get year of growth start
  first_years <- apply(rw, 2, function(col) {
    years <- rownames(rw)
    first_non_na <- which(!is.na(col))[1]  # index of first non-NA
    if (!is.na(first_non_na)) {
      return(years[first_non_na])
    } else {
      return(NA)  # if the whole column is NA
    }
  })
  
  names(first_years) <- substring(names(first_years), 4)
  first_years
  # create col in meta data to match chronology data -> join
  file_meta_ump$sample_ind <- paste(file_meta_ump$sample, file_meta_ump$ind, sep = "_")
  file_meta_ump$first_year <- first_years[file_meta_ump$sample_ind]
  # create new col containing years after growth of tree
  file_meta_ump$year <- as.numeric(file_meta_ump$year)
  file_meta_ump$first_year <- as.numeric(file_meta_ump$first_year)
  file_meta_ump$growth_years <- file_meta_ump$year - file_meta_ump$first_year
  
  # filter for heartwood spectra (<= 10 years of growth)
  #keep_rows <- file_meta_ump$growth_years > 14
  #keep_rows <- !(file_meta_ump$ring == "EW" & file_meta_ump$growth_years < 15)
  #file_meta_ump <- file_meta_ump[keep_rows, ]
  #data_mat_ump  <- data_mat_ump[keep_rows, ]
  
  # La Fuente
  laf_ix <- which(file_meta$loc == "LaF")
  file_meta_laf <- file_meta[laf_ix, ]
  data_mat_laf <- data_mat[laf_ix, ]
  
  rw <- read.rwl("Bursera_La_Fuente.rwl")
  # for each tree get year of growth start
  first_years <- apply(rw, 2, function(col) {
    years <- rownames(rw)
    first_non_na <- which(!is.na(col))[1]  # index of first non-NA
    if (!is.na(first_non_na)) {
      return(years[first_non_na])
    } else {
      return(NA)  # if the whole column is NA
    }
  })
  
  names(first_years) <- substring(names(first_years), 4)
  first_years
  # create col in meta data to match chronology data -> join
  file_meta_laf$sample_ind <- paste0("A", file_meta_laf$sample, "_", file_meta_laf$ind)#, "_") # according to .rwl
  file_meta_laf$first_year <- first_years[file_meta_laf$sample_ind]
  # create new col containing years after growth of tree
  file_meta_laf$year <- as.numeric(file_meta_laf$year)
  file_meta_laf$first_year <- as.numeric(file_meta_laf$first_year)
  file_meta_laf$growth_years <- file_meta_laf$year - file_meta_laf$first_year
  
  # filter for heartwood spectra (<= 10 years of growth)
  #keep_rows <- file_meta_laf$growth_years > 14
  #keep_rows <- !(file_meta_laf$ring == "EW" & file_meta_laf$growth_years < 15)
  
  #file_meta_laf <- file_meta_laf[keep_rows, ]
  #data_mat_laf  <- data_mat_laf[keep_rows, ]
  
  #------
  # recombine
  file_meta_ump$year <- as.character(file_meta_ump$year)
  file_meta_laf$year <- as.character(file_meta_laf$year)
  
  # use bind_rows to automatically fill missing columns with NA
  file_meta <- bind_rows(file_meta_ump, file_meta_laf)
  data_mat <- bind_rows(as.data.frame(data_mat_ump), as.data.frame(data_mat_laf))
  
  # check dimensions
  print(nrow(file_meta_ump)    == nrow(data_mat_ump))
  print(nrow(file_meta_laf)    == nrow(data_mat_laf))
  print(file_meta)
  #-----------
  # write to file
  cat("Writing data...\n")
  write_sampled_data(file_meta, data_mat, "data/FTIR_Spectra/PLSR", "Bursera_simaruba", method = "PLSR")
  
  cat("Data preparation for Bursera simaruba data successful. \nReady for PLSR age prediction.")
}



preprocess_bursera_valid_years <- function(path_data, path_out) {
  #------------
  # read data and meta
  spec_wn <- read_spec(path_data)
  data_mat <- spec_wn$spec
  wn <- spec_wn$wn
  file_meta <- read_file_metadat(path_data)
  #------------
  # remove not suitable spectra
  #bad_comments <- c("weakSignal", "plusminus2jahre", "weirdShape 1500", "noAbsOver03")
  #rows_to_keep <- !file_meta$comment %in% bad_comments
  #file_meta <- file_meta[rows_to_keep, ]
  #data_mat <- data_mat[rows_to_keep, ]
  
  #------------
  # divide file_meta by locations Umpalá and La Fuente and
  # calc age of tree of respective spectrum if a year is provided
  # (division by location is necessary because samples have same labels)
  
  # filter meta data for spectra that contain valid year information
  year_ix <- which(file_meta$year != "noYear")
  file_meta <- file_meta[year_ix, ]
  data_mat <- data_mat[year_ix, ]
  
  ### calc age
  # Umpalá
  ump_ix <- which(file_meta$loc == "Ump")
  file_meta_ump <- file_meta[ump_ix, ]
  data_mat_ump <- data_mat[ump_ix, ]
  
  rw <- read.rwl("Bursera_Umpala.rwl")
  # for each tree get year of growth start
  first_years <- apply(rw, 2, function(col) {
    years <- rownames(rw)
    first_non_na <- which(!is.na(col))[1]  # ix of first non-NA
    if (!is.na(first_non_na)) {
      return(years[first_non_na])
    } else {
      return(NA)  # if the whole column is NA
    }
  })
  
  names(first_years) <- substring(names(first_years), 4)
  first_years
  # create col in meta data to match chronology data -> join
  file_meta_ump$sample_ind <- paste(file_meta_ump$sample, file_meta_ump$ind, sep = "_")
  file_meta_ump$first_year <- first_years[file_meta_ump$sample_ind]
  # create new col containing years after growth of tree
  file_meta_ump$year <- as.numeric(file_meta_ump$year)
  file_meta_ump$first_year <- as.numeric(file_meta_ump$first_year)
  file_meta_ump$growth_years <- file_meta_ump$year - file_meta_ump$first_year
  
  # filter for heartwood spectra (<= 10 years of growth)
  #keep_rows <- file_meta_ump$growth_years > 14
  #keep_rows <- !(file_meta_ump$ring == "EW" & file_meta_ump$growth_years < 15)
  #file_meta_ump <- file_meta_ump[keep_rows, ]
  #data_mat_ump  <- data_mat_ump[keep_rows, ]
  
  
  # La Fuente
  laf_ix <- which(file_meta$loc == "LaF")
  file_meta_laf <- file_meta[laf_ix, ]
  data_mat_laf <- data_mat[laf_ix, ]
  
  rw <- read.rwl("Bursera_La_Fuente.rwl")
  # for each tree get year of growth start
  first_years <- apply(rw, 2, function(col) {
    years <- rownames(rw)
    first_non_na <- which(!is.na(col))[1]  # index of first non-NA
    if (!is.na(first_non_na)) {
      return(years[first_non_na])
    } else {
      return(NA)  # if the whole column is NA
    }
  })
  
  names(first_years) <- substring(names(first_years), 4)
  first_years
  # create col in meta data to match chronology data -> join
  file_meta_laf$sample_ind <- paste0("A", file_meta_laf$sample, "_", file_meta_laf$ind)#, "_") # according to .rwl
  file_meta_laf$first_year <- first_years[file_meta_laf$sample_ind]
  # create new col containing years after growth of tree
  file_meta_laf$year <- as.numeric(file_meta_laf$year)
  file_meta_laf$first_year <- as.numeric(file_meta_laf$first_year)
  file_meta_laf$growth_years <- file_meta_laf$year - file_meta_laf$first_year
  
  # filter for heartwood spectra (<= 10 years of growth)
  #keep_rows <- file_meta_laf$growth_years > 14
  #keep_rows <- !(file_meta_laf$ring == "EW" & file_meta_laf$growth_years < 15)
  
  #file_meta_laf <- file_meta_laf[keep_rows, ]
  #data_mat_laf  <- data_mat_laf[keep_rows, ]
  #------
  # recombine
  file_meta_ump$year <- as.character(file_meta_ump$year)
  file_meta_laf$year <- as.character(file_meta_laf$year)
  
  # use bind_rows to automatically fill missing columns with NA
  file_meta <- bind_rows(file_meta_ump, file_meta_laf)
  data_mat <- bind_rows(as.data.frame(data_mat_ump), as.data.frame(data_mat_laf))
  
  # check dimensions
  print(nrow(file_meta_ump)    == nrow(data_mat_ump))
  print(nrow(file_meta_laf)    == nrow(data_mat_laf))
  print(file_meta)
  #-----------
  # write to file
  cat("Writing data...\n")
  write_sampled_data(file_meta, data_mat, "data/FTIR_Spectra/Bands_over_time", "Bursera_simaruba", method = "valid_years")
  
  cat("Data preparation for Bursera simaruba data successful. \nReady for plotting band values over time")
}


### wrapper
preprocess_bursera <- function(path_data, path_out, model, heartwood_years = NA) {
  if (model == "PLS-DA") {
    preprocess_bursera_PLSDA(path_data, path_out, heartwood_years)
  } else if (model == "PLSR") {
    preprocess_bursera_PLSR(path_data, path_out)
  }  else if (model == "valid_years") {
    preprocess_bursera_valid_years(path_data, path_out)
  } else {
    stop("Must choose PLS-DA or PLSR.")
  }
}