#==============================================================================
# Supplementary R scripts for mater's thesis
# Title: Potential of FTIR spectroscopy for tropical dendrochronology
# Author: Lukas Erzfeld
# Insitution: University of Leipzig
# Date: August 26, 2025
#===============================================================================

remove_mixed_spectra <- function(metadata, spectra_matrix) {
  metadata$comment[is.na(metadata$comment)] <- "valid"
  
  # if applicable:
  # remove Mixed spectra
  mixed_ix <- which(metadata$comment == "Mixed")
  if (length(mixed_ix) > 0) {
    metadata <- metadata[-mixed_ix, ]
    spectra_matrix <- spectra_matrix[-mixed_ix, ]
  }
  # remove non-clear spectra
  non_clear_ix <- which(metadata$comment == "notClear")
  if (length(non_clear_ix) > 0) {
    metadata <- metadata[-non_clear_ix, ]
    spectra_matrix <- spectra_matrix[-non_clear_ix, ]
  }
  # remove slightly mixed EW spectra
  pot_mixed_EW_ix <- which(metadata$comment == "sliMixed" & metadata$ring == "EW")
  if (length(pot_mixed_EW_ix) > 0) {
    metadata <- metadata[-pot_mixed_EW_ix, ]
    spectra_matrix <- spectra_matrix[-pot_mixed_EW_ix, ]
  }
  # remove potentially mixed EW spectra
  pot_mixed_EW_ix <- which(metadata$comment == "potslilMixed" & metadata$ring == "EW")
  if (length(pot_mixed_EW_ix) > 0) {
    metadata <- metadata[-pot_mixed_EW_ix, ]
    spectra_matrix <- spectra_matrix[-pot_mixed_EW_ix, ]
  }
  
  return(list(metadata = metadata, spectra_matrix = spectra_matrix))
}

filter_heartwood <- function(metadata, spectra_matrix, first_year_dict, heartwood_years=NA) {
  # create sample_ind column to match keys in dictionary
  metadata$sample_ind <- paste(metadata$sample, metadata$ind, sep = "_")
  metadata$year <- as.numeric(metadata$year)
  
  # assign first_year and compute growth_years
  metadata$first_year <- as.numeric(first_year_dict[metadata$sample_ind])
  metadata$growth_years <- metadata$year - metadata$first_year
  
  # filter: keep EW with growth_years >= 6, and all LW
  #keep_rows <- !(metadata$ring == "EW" & metadata$growth_years <= heartwood_years)
  keep_rows <- !(metadata$growth_years <= heartwood_years)
  
  metadata <- metadata[keep_rows, ]
  spectra_matrix  <- spectra_matrix[keep_rows, ]
  
  return(list(metadata = metadata, spectra_matrix = spectra_matrix))
}

calc_age_from_dict <- function(metadata, first_year_dict) {
  # create sample_ind column to match keys in dictionary
  metadata$sample_ind <- paste(metadata$sample, metadata$ind, sep = "_")
  metadata$year <- as.numeric(metadata$year)
  
  # assign first_year and compute growth_years
  metadata$first_year <- as.numeric(first_year_dict[metadata$sample_ind])
  metadata$growth_years <- metadata$year - metadata$first_year
  return(metadata)
}

filter_heartwood <- function(metadata, spectra_matrix, heartwood_years = NA) {
  # if heartwood_years is not NA:
  # filter: keep EW with growth_years >= heartwood_years, and all LW
  if (!is.na(heartwood_years)) {
    keep_rows <- !(metadata$ring == "EW" & metadata$growth_years <= heartwood_years)
    metadata <- metadata[keep_rows, ]
    spectra_matrix <- spectra_matrix[keep_rows, ]
  }
  
  return(list(metadata = metadata, spectra_matrix = spectra_matrix))
}

balance_EW_LW <- function(metafile, spectra_matrix) {
  # indices of each class
  ix_LW <- which(metafile$ring == "LW")
  ix_EW <- which(metafile$ring == "EW")
  
  # randomly sample as many "EW" entries as there are "LW" to obtain a balanced dataset for classification
  set.seed(127)  # for reproducibility
  ix_EW_sampled <- sample(ix_EW, size = length(ix_LW))
  
  # combine to create balanced index
  ix_balanced <- c(ix_LW, ix_EW_sampled)
  
  # subset metadata and spectra
  metafile <- metafile[ix_balanced, ]
  spectra_matrix  <- spectra_matrix[ix_balanced, ]
  
  return(list(metafile = metafile, spectra_matrix = spectra_matrix))
}

write_sampled_data <- function(metafile, spectra_matrix, out_path, species_name, method) {
  # write sampled data to file for easier handling in the future
  file_meta_write <- metafile[, c(1:7, 9:10)]
  if (method == "PLSR") {
    file_meta_write <- metafile[, c(1:7, 9:10)]
    write.csv(file_meta_write, file.path(out_path, paste0(species_name, "_file_meta_PLSR.csv")), row.names = FALSE)
    write.csv(spectra_matrix, file.path(out_path, paste0(species_name, "_data_mat_PLSR.csv")), row.names = FALSE)
  } else if (method == "PLS-DA") {
    file_meta_write <- metafile[, 1:7]
    write.csv(file_meta_write, file.path(out_path, paste0(species_name, "_file_meta_PLS-DA.csv")), row.names = FALSE)
    write.csv(spectra_matrix, file.path(out_path, paste0(species_name, "_data_mat_PLS-DA.csv")), row.names = FALSE)
  } else if (method == "valid_years") {
    file_meta_write <- metafile[,]
    write.csv(file_meta_write, file.path(out_path, paste0(species_name, "_file_meta_valid_years.csv")), row.names = FALSE)
    write.csv(spectra_matrix, file.path(out_path, paste0(species_name, "_data_mat_valid_years.csv")), row.names = FALSE)
  } else {
    cat("Method must be provided.")
  }
}
