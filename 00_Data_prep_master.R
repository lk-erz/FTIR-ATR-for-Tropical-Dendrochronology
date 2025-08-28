#==============================================================================
# Supplementary R scripts for mater's thesis
# Title: Potential of FTIR spectroscopy for tropical dendrochronology
# Author: Lukas Erzfeld
# Insitution: University of Leipzig
# Date: August 26, 2025
#===============================================================================
source("00_PLS-DA_data_prep_Bursera.R")
source("00_PLS-DA_data_prep_Maclura_Pinus.R")
#-------------------------------------------------------------------------------
# specify paths
path_out <- "data/FTIR_Spectra/PLS-DA/"
path_bursera <- "data/FTIR_Spectra/Bursera_simaruba/"
path_maclura <- "data/FTIR_Spectra/Maclura_tinctoria/"
path_pinus <- "data/FTIR_Spectra/Pinus_patula/"

# data preparation PLS-DA
preprocess_bursera(path_data = path_bursera, 
                         path_out = path_out,
                         heartwood_years = 15,
                         model = "PLS-DA")

preprocess_maclura_pinus(path_data = path_maclura, 
                         path_out = path_out, 
                         first_year_dict = c("1_1" = 1954, "2_1" = 1959, "3_1" = 1961, "3_2" = 1959), # based on ring count
                         heartwood_years = 5,
                         species_name = "Maclura_tinctoria",
                         model = "PLS-DA")

preprocess_maclura_pinus(path_data = path_pinus, 
                         path_out = path_out, 
                         first_year_dict = c("1_1" = 2005, "1_2" = 2005, "1_3" = 2005), # based on ring count
                         species_name = "Pinus_patula",
                         model = "PLS-DA")


# data preparation for PLSR
path_out <- "data/FTIR_Spectra/PLSR/"

# data preparation PLSR
preprocess_bursera(path_data = path_bursera, 
                        path_out = path_out,
                        model = "PLSR")

preprocess_maclura_pinus(path_data = path_maclura, 
                         path_out = path_out, 
                         first_year_dict = c("1_1" = 1954, "2_1" = 1959, "3_1" = 1961, "3_2" = 1959), # based on ring count
                         species_name = "Maclura_tinctoria",
                         model = "PLSR")

preprocess_maclura_pinus(path_data = path_pinus, 
                         path_out = path_out, 
                         first_year_dict = c("1_1" = 2005, "1_2" = 2005, "1_3" = 2005), # based on ring count
                         species_name = "Pinus_patula",
                         model = "PLSR")


# data preparation for plots of band values over time
path_out <- "data/FTIR_Spectra/Bands_over_time/"

preprocess_bursera(path_data = path_bursera, 
                   path_out = path_out,
                   model = "valid_years")

preprocess_maclura_pinus_valid_years(path_data = path_maclura,
                                     path_out = path_out, 
                                     first_year_dict = c("1_1" = 1954, "2_1" = 1959, "3_1" = 1961, "3_2" = 1959), # based on ring count
                                     species_name = "Maclura_tinctoria")

preprocess_maclura_pinus_valid_years(path_data = path_pinus, 
                                     path_out = path_out, 
                                     first_year_dict = c("1_1" = 2005, "1_2" = 2005, "1_3" = 2005), # based on ring count
                                     species_name = "Pinus_patula")
