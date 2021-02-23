# *** Header **************************************************************************
#
# Create age-based contact matrix.
#

# *** Read in 5-Year ACS PUMS Data ****************************************************
raw_pums <- dplyr::as_tibble(data.table::fread(
  file = stringr::str_c(
    build_data_dir,
    "/pums_data_with_nasem_phases.csv"
  ),
  sep = ","
))
# Population share by 5-yr age-group from PUMS data
pop_shs <- raw_pums %>%
  select(
    PWGTP,
    AGEP
  ) %>%
  mutate(
    age_bin = case_when(
      AGEP <= 4               ~ "0 to 4",
      AGEP >  4  & AGEP <= 9  ~ "5 to 9",
      AGEP >  9  & AGEP <= 14 ~ "10 to 14",
      AGEP >  14 & AGEP <= 19 ~ "15 to 19",
      AGEP >  19 & AGEP <= 24 ~ "20 to 24",
      AGEP >  24 & AGEP <= 29 ~ "25 to 29",
      AGEP >  29 & AGEP <= 34 ~ "30 to 34",
      AGEP >  34 & AGEP <= 39 ~ "35 to 39",
      AGEP >  39 & AGEP <= 44 ~ "40 to 44",
      AGEP >  44 & AGEP <= 49 ~ "45 to 49",
      AGEP >  49 & AGEP <= 54 ~ "50 to 54",
      AGEP >  54 & AGEP <= 59 ~ "55 to 59",
      AGEP >  59 & AGEP <= 64 ~ "60 to 64",
      AGEP >  64 & AGEP <= 69 ~ "65 to 69",
      AGEP >  69 & AGEP <= 74 ~ "70 to 74",
      AGEP >  74              ~ "75+",
      TRUE ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(age_bin)
  ) %>%
  group_by(age_bin) %>%
  summarize(
    n = sum(PWGTP)
  ) %>%
  mutate(
    sh = n / sum(n)
  ) %>%
  select(
    -n
  )

# *** Prem et al. (2020) Contact Matrix ***********************************************
prem_raw_age_based_contact_mat <- read.csv(
  stringr::str_c(
    raw_data_dir,
    "/age_based_contact_matrix/synthetic_contacts_2020.csv"
  )
)

clean_prem_contact_mat <- as.matrix(prem_raw_age_based_contact_mat %>%
  filter(
    iso3c            == "USA",
    setting          == "overall",
    location_contact == "all"
  ) %>%
  select(
    age_contactor, age_cotactee, mean_number_of_contacts
  ) %>%
  rename(
    age_bin_i = age_contactor,
    age_bin_j = age_cotactee,
    mean_daily_contacts = mean_number_of_contacts
  ) %>%
  select(
    age_bin_i, age_bin_j, mean_daily_contacts
  ) %>%
  tidyr::pivot_wider(
    id_cols     = age_bin_i,
    names_from  = age_bin_j,
    values_from = mean_daily_contacts
  ) %>%
  tibble::column_to_rownames("age_bin_i"))

# *** Collapse from 5-yr to 10-yr, Extrapolate to 80+ *********************************
pop_sh_mat <- as.matrix(pop_shs %>%
  mutate(
    age_bin = factor(age_bin, levels = rownames(clean_prem_contact_mat))
  ) %>%
  arrange(
    age_bin
  ) %>%
  mutate(
    age_bin = as.character(age_bin)
  ) %>%
  tibble::column_to_rownames("age_bin"))

#' Convert to 10-yr per Bubar et al. (2021)
#' 
#' @param: five_yr_contact_mat<mat<dbl>> | 2n x 2n numeric matrix with row/col sorted by
#'                                         five year age bin
#' @param: pop_sh_mat<mat<dbl>>          | 2n x 1 numeric matrix with rows sorted by five
#'                                         year age bin
#' 
#' @return: out_10_yr_mat<mat<dbl>>      | n x n numeric mat of 10-yr collapsed contacts
bubar_method_to_convert_from_5_yr_to_10_yr <- function(five_yr_contact_mat, pop_sh_mat) {
  num_5_yr_bins = NROW(five_yr_contact_mat)
  out_10_yr_mat = matrix(nrow = num_5_yr_bins / 2, ncol = num_5_yr_bins / 2)
  i_hat = 1
  for (i in seq(1, num_5_yr_bins, by = 2)) {
    j_hat = 1
    for (j in seq(1, num_5_yr_bins, by = 2)) {
      p1 = pop_sh_mat[i]
      p2 = pop_sh_mat[i + 1]
      five_yr_i_j       = five_yr_contact_mat[i, j]
      five_yr_ipl1_j    = five_yr_contact_mat[i + 1, j]
      five_yr_i_jpl1    = five_yr_contact_mat[i, j + 1]
      five_yr_ipl1_jpl1 = five_yr_contact_mat[i + 1, j + 1]
      out_10_yr_mat[i_hat, j_hat] = (
        ((five_yr_i_j + five_yr_i_jpl1) * p1 + (five_yr_ipl1_j + five_yr_ipl1_jpl1) * p2) / (p1 + p2) 
      )
      j_hat = j_hat + 1
    }
    i_hat = i_hat + 1
  }
  return(out_10_yr_mat)
}

#' Extrapolate to 80+ per Bubar et al. (2021)
#' 
#' @param: ten_yr_contact_mat<mat<dbl>> | n x n numeric matrix with row/col sorted by
#'                                        10 year age bin
#' @param: scale_factor<dbl>            | amount to scale down 80pl contacts with young ppl [0,1]
#' 
#' @return: out_extrapolated<mat<dbl>>  | n + 1 x n + 1 numeric matrix with extrapolated 80+ row/col
bubar_method_to_extrapoloate_80_pl <- function(ten_yr_contact_mat, scale_factor) {
  num_bins = NROW(ten_yr_contact_mat) + 1
  i80pl  = num_bins
  i70_79 = num_bins - 1
  i60_69 = num_bins - 2
  i50_59 = num_bins - 3
  i40_49 = num_bins - 4
  i30_39 = num_bins - 5
  i20_29 = num_bins - 6
  i10_19 = num_bins - 7
  i0_9   = num_bins - 8
  stopifnot(i0_9 == 1)
  out_extrapolated = cbind(rbind(ten_yr_contact_mat, 0), 0)
  # 80pl col copy
  out_extrapolated[i10_19:i80pl, i80pl] = out_extrapolated[i0_9:i70_79, i70_79]
  out_extrapolated[i0_9, i80pl]         = out_extrapolated[i0_9, i70_79]
  # 80pl row copy
  out_extrapolated[i80pl, i10_19:i80pl] = out_extrapolated[i70_79, i0_9:i70_79]
  out_extrapolated[i80pl, i0_9]         = out_extrapolated[i70_79, i0_9]
  # Redistribute for LTCF (col)
  alpha = scale_factor * sum(out_extrapolated[i0_9:i60_69, i80pl])
  out_extrapolated[i0_9:i60_69, i80pl]  = (1 - scale_factor) * out_extrapolated[i0_9:i60_69, i80pl] 
  out_extrapolated[i70_79:i80pl, i80pl] = out_extrapolated[i70_79:i80pl, i80pl] + (alpha / 2)
  # Redistribute for LTCF (row)
  alpha = scale_factor * sum(out_extrapolated[i80pl, i0_9:i60_69])
  out_extrapolated[i80pl, i0_9:i60_69]  = (1 - scale_factor) * out_extrapolated[i80pl, i0_9:i60_69] 
  out_extrapolated[i80pl, i70_79:i80pl] = out_extrapolated[i80pl, i70_79:i80pl] + (alpha / 2)
  
  new_age_bins = c(
    "agebin_0-9", "agebin_10-19", "agebin_20-29", "agebin_30-39",
    "agebin_40-49", "agebin_50-59", "agebin_60-69", "agebin_70-79", "agebin_80+"
  )
  colnames(out_extrapolated) = new_age_bins
  rownames(out_extrapolated) = new_age_bins

  return(out_extrapolated)
}

# *** Write to CSV ********************************************************************
contact_mat <- bubar_method_to_extrapoloate_80_pl(
  bubar_method_to_convert_from_5_yr_to_10_yr(
    clean_prem_contact_mat, pop_sh_mat
  ),
  scale_factor = .10
)

write.csv(
  contact_mat,
  stringr::str_c(
    build_data_dir,
    "/age_group_contact_matrix.csv"
  )
)