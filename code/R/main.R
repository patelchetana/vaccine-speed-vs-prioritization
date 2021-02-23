# *** Header ***************************************************************************
#
# Master script to build inputs for SEIR model.

# *** Load Packages ********************************************************************
library(dplyr)
library(stringr)
library(feather)
library(assertr)
library(tidyr)
library(data.table)
library(zoo)
library(openxlsx)
library(survey)
library(SASxport)

# Clear workspace
rm(list = ls())
# Turn off scientific notation
options(scipen = 999)

# *** Set Directories ******************************************************************                                     v
base_dir              <- "."
raw_data_dir          <- stringr::str_c(base_dir, "/raw_data")
build_data_dir        <- stringr::str_c(base_dir, "/inputs")
exhibit_data_dir      <- stringr::str_c(base_dir, "/exhibits")

# *** Switches *************************************************************************
build_nasem_tagged_pums_dataset                                                    <- 1
build_contact_matrix                                                               <- 1
build_national_baseline_file                                                       <- 1
build_seir_initial_conditions                                                      <- 1
run_nasem_allocation_simulations                                                   <- 1
generate_supplementary_materials_table_s1                                          <- 1
generate_supplementary_materials_table_s2                                          <- 1
generate_supplementary_materials_table_s3                                          <- 1

# *** Execution ************************************************************************
if (build_nasem_tagged_pums_dataset == 1) {
  print("Beginning script: b1_create_nasem_tagged_pums_dataset.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/input_building_scripts/b1_create_nasem_tagged_pums_dataset.R"
    )
  )
  print("Completed script: b1_create_nasem_tagged_pums_dataset.R")
}

if (build_contact_matrix == 1) {
  print("Beginning script: b2_create_contact_matrix.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/input_building_scripts/b2_create_contact_matrix.R"
    )
  )
  print("Completed script: b2_create_contact_matrix.R")
}

if (build_national_baseline_file == 1) {
  print("Beginning script: b3_create_national_baseline_file.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/input_building_scripts/b3_create_national_baseline_file.R"
    )
  )
  print("Completed script: b3_create_contact_matrix.R")
}

if (build_seir_initial_conditions == 1) {
  print("Beginning script: b4_build_seir_initial_conditions.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/input_building_scripts/b4_create_seir_initial_conditions.R"
    )
  ) 
  print("Completed script: b4_build_seir_initial_conditions.R")
}

if (run_nasem_allocation_simulations == 1) {
  print("Beginning script: b5_run_unreserved_nasem_allocation_simulations_with_uptake.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/input_building_scripts/b5_run_unreserved_nasem_allocation_simulations_with_uptake.R"
    )
  )     
  print("Completed script: b5_run_unreserved_nasem_allocation_simulations_with_uptake.R")
}

if (generate_supplementary_materials_table_s1 == 1) {
  print("Beginning script: b1_supplementary_materials_table_s1.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/supplementary_materials_building_scripts/b1_supplementary_materials_table_s1.R"
    )
  )     
  print("Completed script: b1_supplementary_materials_table_s1.R")
}

if (generate_supplementary_materials_table_s2 == 1) {
  print("Beginning script: b2_supplementary_materials_table_s2.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/supplementary_materials_building_scripts/b2_supplementary_materials_table_s2.R"
    )
  )     
  print("Completed script: b2_supplementary_materials_table_s2.R")
}

if (generate_supplementary_materials_table_s3 == 1) {
  print("Beginning script: b3_supplementary_materials_table_s3.R")
  source(
    stringr::str_c(
      base_dir,
      "/code/R/supplementary_materials_building_scripts/b3_supplementary_materials_table_s3.R"
    )
  )     
  print("Completed script: b3_supplementary_materials_table_s3.R")
}