# *** Header **************************************************************************
#
# Create national baseline file.
#

# *** Population Shares from 5-yr ACS PUMS ********************************************
raw_pums <- dplyr::as_tibble(data.table::fread(
  file = stringr::str_c(
    build_data_dir,
    "/pums_data_with_nasem_phases.csv"
  ),
  sep = ","
))

pop_sh_by_age_bin <- raw_pums %>%
  select(
    PWGTP,
    AGEP
  ) %>%
  mutate(
    group = case_when(
      AGEP < 10              ~ "agebin_0-9",
      AGEP < 20 & AGEP >= 10 ~ "agebin_10-19",
      AGEP < 30 & AGEP >= 20 ~ "agebin_20-29",
      AGEP < 40 & AGEP >= 30 ~ "agebin_30-39",
      AGEP < 50 & AGEP >= 40 ~ "agebin_40-49",
      AGEP < 60 & AGEP >= 50 ~ "agebin_50-59",
      AGEP < 70 & AGEP >= 60 ~ "agebin_60-69",
      AGEP < 80 & AGEP >= 70 ~ "agebin_70-79",
      AGEP >= 80             ~ "agebin_80+",
      TRUE                   ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(group)
  ) %>%
  select(-AGEP) %>%
  group_by(group) %>%
  summarize(
    n_group = sum(PWGTP)
  ) %>%
  mutate(
    sh_group = n_group / sum(n_group),
    nat_population = sum(n_group)
  )

# *** Infection Fatality Rates (IFRS) *************************************************
ifr_by_age_bin <- pop_sh_by_age_bin %>%
  select(group) %>%
  mutate(
    # Levin et al. (2020) IFRs from Supplementary Appendix Q
    # (https://static-content.springer.com/esm/art%3A10.1007%2Fs10654-020-00698-1/MediaObjects/10654_2020_698_MOESM1_ESM.pdf)
    c19_ifr_group = case_when(
      group == "agebin_0-9"    ~ .00001,
      group == "agebin_10-19"  ~ .00003,
      group == "agebin_20-29"  ~ .00011,
      group == "agebin_30-39"  ~ .00037,
      group == "agebin_40-49"  ~ .00123,
      group == "agebin_50-59"  ~ .00413,
      group == "agebin_60-69"  ~ .01380,
      group == "agebin_70-79"  ~ .04620,
      group == "agebin_80+"    ~ .15460,
      TRUE                     ~ NA_real_
    )
  ) %>%
  assertr::verify(
    ! is.na(c19_ifr_group)
  ) %>%
  mutate(
    c19_ifr_group = c19_ifr_group * 100
  )

# *** Add Susceptibility (All + Clinical) to Infection ********************************
sus_to_inf_by_age_bin <- ifr_by_age_bin %>%
  select(group) %>%
  mutate(
    # Davies et al. (2020) from Extended Data Fig. 4
    # https://www.nature.com/articles/s41591-020-0962-9/figures/8
    sus_to_inf = case_when(
      group == "agebin_0-9"    ~ .4,
      group == "agebin_10-19"  ~ .38,
      group == "agebin_20-29"  ~ .79,
      group == "agebin_30-39"  ~ .86,
      group == "agebin_40-49"  ~ .8,
      group == "agebin_50-59"  ~ .82,
      group == "agebin_60-69"  ~ .88,
      group == "agebin_70-79"  ~ .74,
      group == "agebin_80+"    ~ .74,
      TRUE                     ~ NA_real_
    ),
    clinical_sus = case_when(
      group == "agebin_0-9"    ~ .29,
      group == "agebin_10-19"  ~ .21,
      group == "agebin_20-29"  ~ .27,
      group == "agebin_30-39"  ~ .33,
      group == "agebin_40-49"  ~ .4,
      group == "agebin_50-59"  ~ .49,
      group == "agebin_60-69"  ~ .63,
      group == "agebin_70-79"  ~ .69,
      group == "agebin_80+"    ~ .69,
      TRUE                     ~ NA_real_
    )
  ) %>%
  assertr::verify(
    ! is.na(sus_to_inf) & ! is.na(clinical_sus)
  )

# *** Vaccine Uptake ******************************************************************
# US Census Bureau
raw_census_uptake_survey <- openxlsx::read.xlsx(
  stringr::str_c(
    raw_data_dir,
    "/vaccine_hesitancy/health5_week22.xlsx"
  ), sheet = "US"
)

clean_census_uptake_survey <- raw_census_uptake_survey %>%
  slice(
    4:13
  ) %>%
  select(
    `Health.Table.5..COVID-19.Vaccinations.by.Select.Characteristics:.United.States`,
    X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12
  ) %>%
  rename(
    agebin = `Health.Table.5..COVID-19.Vaccinations.by.Select.Characteristics:.United.States`,
    total_adult_pop                = X2,
    got_vax_pop_total              = X3,
    got_vax_pop_both_dose          = X4,
    got_vax_pop_not_both_dose      = X5,
    got_vax_pop_missing_data       = X6,
    no_vax_pop_total               = X7,
    no_vax_pop_def_get_vax         = X8,
    no_vax_pop_prob_get_vax        = X9,
    no_vax_pop_prob_not_get_vax    = X10,
    no_vax_pop_def_not_get_vax     = X11,
    no_vax_pop_missing_data        = X12,
  ) %>%
  filter(
    ! is.na(agebin) & trimws(agebin) != "Age" & agebin != "Select characteristics" & trimws(agebin) != "Total"
  ) %>%
  mutate(
    group = case_when(
      trimws(agebin) == "18 - 24"      ~ "18_24",
      trimws(agebin) == "25 - 39"      ~ "25_39",
      trimws(agebin) == "40 - 54"      ~ "40_54",
      trimws(agebin) == "55 - 64"      ~ "55_64",
      trimws(agebin) == "65 and above" ~ "65+"
    ),
    # One-off non-numeric entry
    got_vax_pop_missing_data = case_when(
      got_vax_pop_missing_data == "-" ~ "0",
      TRUE ~ got_vax_pop_missing_data
    ),
    total_adult_pop               = as.numeric(total_adult_pop),
    got_vax_pop_total             = as.numeric(got_vax_pop_total),
    got_vax_pop_both_dose         = as.numeric(got_vax_pop_both_dose),
    got_vax_pop_not_both_dose     = as.numeric(got_vax_pop_not_both_dose),
    got_vax_pop_missing_data      = as.numeric(got_vax_pop_missing_data),
    no_vax_pop_total              = as.numeric(no_vax_pop_total),
    no_vax_pop_def_get_vax        = as.numeric(no_vax_pop_def_get_vax),
    no_vax_pop_prob_get_vax       = as.numeric(no_vax_pop_prob_get_vax),
    no_vax_pop_prob_not_get_vax   = as.numeric(no_vax_pop_prob_not_get_vax),
    no_vax_pop_def_not_get_vax    = as.numeric(no_vax_pop_def_not_get_vax),
    no_vax_pop_missing_data       = as.numeric(no_vax_pop_missing_data)
  ) %>%
  select(-agebin) %>%
  assertr::verify(
    no_vax_pop_total == (
      no_vax_pop_def_get_vax +
      no_vax_pop_prob_get_vax +
      no_vax_pop_prob_not_get_vax +
      no_vax_pop_def_not_get_vax +
      no_vax_pop_missing_data
    )
  ) %>%
  mutate(
    sh_census_vax_uptake = (
      (no_vax_pop_def_get_vax + no_vax_pop_prob_get_vax + got_vax_pop_total) 
      / (no_vax_pop_total - no_vax_pop_missing_data + got_vax_pop_total))
  ) %>%
  select(
    group, sh_census_vax_uptake
  )

vax_uptake_by_age_bin <- raw_pums %>%
  select(
    PWGTP,
    AGEP
  ) %>%
  mutate(
    fine_age_bin = case_when(
      AGEP <= 4               ~ "0 to 4",
      AGEP  > 4  & AGEP <= 9  ~ "5 to 9",
      AGEP  > 9  & AGEP <= 15 ~ "10 to 15",
      AGEP  > 15 & AGEP <= 19 ~ "16 to 19",
      AGEP  > 19 & AGEP <= 24 ~ "20 to 24",
      AGEP  > 24 & AGEP <= 29 ~ "25 to 29",
      AGEP  > 29 & AGEP <= 34 ~ "30 to 34",
      AGEP  > 34 & AGEP <= 39 ~ "35 to 39",
      AGEP  > 39 & AGEP <= 44 ~ "40 to 44",
      AGEP  > 44 & AGEP <= 49 ~ "45 to 49",
      AGEP  > 49 & AGEP <= 54 ~ "50 to 54",
      AGEP  > 54 & AGEP <= 59 ~ "55 to 59",
      AGEP  > 59 & AGEP <= 64 ~ "60 to 64",
      AGEP  > 64 & AGEP <= 69 ~ "65 to 69",
      AGEP  > 69 & AGEP <= 74 ~ "70 to 74",
      AGEP  > 74              ~ "75+",
      TRUE  ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(fine_age_bin)
  ) %>%
  select(-AGEP) %>%
  group_by(fine_age_bin) %>%
  summarize(
    n = sum(PWGTP)
  ) %>%
  mutate(
    sh = n / sum(n)
  ) %>%
  mutate(
    # Table 1 of Malik et al (2020)
    # https://www.thelancet.com/journals/eclinm/article/PIIS2589-5370(20)30239-X/fulltext
    # 'all_ages' projects youngest surveyed age-group to 0-18 year olds
    all_ages_malik_et_al_survey_vax_uptake = case_when(
      fine_age_bin == "0 to 4"   ~ .59,
      fine_age_bin == "5 to 9"   ~ .59,
      fine_age_bin == "10 to 15" ~ .59,
      fine_age_bin == "16 to 19" ~ .59,
      fine_age_bin == "20 to 24" ~ .59,
      fine_age_bin == "25 to 29" ~ .60,
      fine_age_bin == "30 to 34" ~ .60,
      fine_age_bin == "35 to 39" ~ .64,
      fine_age_bin == "40 to 44" ~ .64,
      fine_age_bin == "45 to 49" ~ .56,
      fine_age_bin == "50 to 54" ~ .56,
      fine_age_bin == "55 to 59" ~ .78,
      fine_age_bin == "60 to 64" ~ .78,
      fine_age_bin == "65 to 69" ~ .78,
      fine_age_bin == "70 to 74" ~ .78,
      fine_age_bin == "75+"      ~ .78,
      TRUE ~ NA_real_
    ),
    # US Census Household Survey Week 22
    # 'all_ages' projects youngest surveyed age-group to 0-18 year olds
    all_ages_census_uptake_survey = case_when(
      fine_age_bin == "0 to 4"   ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "5 to 9"   ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "10 to 15" ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "16 to 19" ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "20 to 24" ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "25 to 29" ~ filter(clean_census_uptake_survey, group == "25_39")[1,2],
      fine_age_bin == "30 to 34" ~ filter(clean_census_uptake_survey, group == "25_39")[1,2],
      fine_age_bin == "35 to 39" ~ filter(clean_census_uptake_survey, group == "25_39")[1,2],
      fine_age_bin == "40 to 44" ~ filter(clean_census_uptake_survey, group == "40_54")[1,2],
      fine_age_bin == "45 to 49" ~ filter(clean_census_uptake_survey, group == "40_54")[1,2],
      fine_age_bin == "50 to 54" ~ filter(clean_census_uptake_survey, group == "40_54")[1,2],
      fine_age_bin == "55 to 59" ~ filter(clean_census_uptake_survey, group == "55_64")[1,2],
      fine_age_bin == "60 to 64" ~ filter(clean_census_uptake_survey, group == "55_64")[1,2],
      fine_age_bin == "65 to 69" ~ filter(clean_census_uptake_survey, group == "65+")[1,2],
      fine_age_bin == "70 to 74" ~ filter(clean_census_uptake_survey, group == "65+")[1,2],
      fine_age_bin == "75+"      ~ filter(clean_census_uptake_survey, group == "65+")[1,2],
    ),
    # Table 1 of Malik et al (2020). Assume 0-15 group uptake is 0 and 16-17 uptake is same as 18-24 age-group
    # https://www.thelancet.com/journals/eclinm/article/PIIS2589-5370(20)30239-X/fulltext
    malik_et_al_survey_vax_uptake = case_when(
      fine_age_bin == "0 to 4"   ~ 0,
      fine_age_bin == "5 to 9"   ~ 0,
      fine_age_bin == "10 to 15" ~ 0,
      fine_age_bin == "16 to 19" ~ .59,
      fine_age_bin == "20 to 24" ~ .59,
      fine_age_bin == "25 to 29" ~ .60,
      fine_age_bin == "30 to 34" ~ .60,
      fine_age_bin == "35 to 39" ~ .64,
      fine_age_bin == "40 to 44" ~ .64,
      fine_age_bin == "45 to 49" ~ .56,
      fine_age_bin == "50 to 54" ~ .56,
      fine_age_bin == "55 to 59" ~ .78,
      fine_age_bin == "60 to 64" ~ .78,
      fine_age_bin == "65 to 69" ~ .78,
      fine_age_bin == "70 to 74" ~ .78,
      fine_age_bin == "75+"      ~ .78,
      TRUE ~ NA_real_
    ),
    # US Census Burea Household Survey Week 22
    # Assume 0-15 year olds have 0 uptake
    census_uptake_survey = case_when(
      fine_age_bin == "0 to 4"   ~ 0,
      fine_age_bin == "5 to 9"   ~ 0,
      fine_age_bin == "10 to 15" ~ 0,
      fine_age_bin == "16 to 19" ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "20 to 24" ~ filter(clean_census_uptake_survey, group == "18_24")[1,2],
      fine_age_bin == "25 to 29" ~ filter(clean_census_uptake_survey, group == "25_39")[1,2],
      fine_age_bin == "30 to 34" ~ filter(clean_census_uptake_survey, group == "25_39")[1,2],
      fine_age_bin == "35 to 39" ~ filter(clean_census_uptake_survey, group == "25_39")[1,2],
      fine_age_bin == "40 to 44" ~ filter(clean_census_uptake_survey, group == "40_54")[1,2],
      fine_age_bin == "45 to 49" ~ filter(clean_census_uptake_survey, group == "40_54")[1,2],
      fine_age_bin == "50 to 54" ~ filter(clean_census_uptake_survey, group == "40_54")[1,2],
      fine_age_bin == "55 to 59" ~ filter(clean_census_uptake_survey, group == "55_64")[1,2],
      fine_age_bin == "60 to 64" ~ filter(clean_census_uptake_survey, group == "55_64")[1,2],
      fine_age_bin == "65 to 69" ~ filter(clean_census_uptake_survey, group == "65+")[1,2],
      fine_age_bin == "70 to 74" ~ filter(clean_census_uptake_survey, group == "65+")[1,2],
      fine_age_bin == "75+"      ~ filter(clean_census_uptake_survey, group == "65+")[1,2],
    )
  ) %>%
  assertr::verify(
    (
      ! is.na(all_ages_malik_et_al_survey_vax_uptake) & ! is.na(all_ages_census_uptake_survey)
      & ! is.na(malik_et_al_survey_vax_uptake) & ! is.na(census_uptake_survey)
    )
  ) %>%
  # ------------ Convert to SEIR age-groups ------------
  mutate(
    filler_col = "Filler"
  ) %>%
  tidyr::pivot_wider(
    id_cols = filler_col,
    names_from = fine_age_bin,
    values_from = c(n, sh,
      all_ages_malik_et_al_survey_vax_uptake,
      all_ages_census_uptake_survey,
      malik_et_al_survey_vax_uptake,
      census_uptake_survey
    ),
    names_glue = "{.value}_{fine_age_bin}"
  ) %>%
  mutate(
    # Malik et al. (2020) All Ages assumption
    `agebin_0-9__vax_uptake_malik_all_ages`   = `all_ages_malik_et_al_survey_vax_uptake_0 to 4`,
    `agebin_10-19__vax_uptake_malik_all_ages` = `all_ages_malik_et_al_survey_vax_uptake_10 to 15`,
    `agebin_20-29__vax_uptake_malik_all_ages` = (
      (`all_ages_malik_et_al_survey_vax_uptake_20 to 24` * `sh_20 to 24` +
       `all_ages_malik_et_al_survey_vax_uptake_25 to 29` * `sh_25 to 29`) / (`sh_20 to 24` + `sh_25 to 29`) 
    ),
    `agebin_30-39__vax_uptake_malik_all_ages` = (
      (`all_ages_malik_et_al_survey_vax_uptake_30 to 34` * `sh_30 to 34` +
       `all_ages_malik_et_al_survey_vax_uptake_35 to 39` * `sh_35 to 39`) / (`sh_30 to 34` + `sh_35 to 39`)
    ),
    `agebin_40-49__vax_uptake_malik_all_ages` = (
      (`all_ages_malik_et_al_survey_vax_uptake_40 to 44` * `sh_40 to 44` +
       `all_ages_malik_et_al_survey_vax_uptake_45 to 49` * `sh_45 to 49`) / (`sh_40 to 44` + `sh_45 to 49`)
    ),
    `agebin_50-59__vax_uptake_malik_all_ages` = (
      (`all_ages_malik_et_al_survey_vax_uptake_50 to 54` * `sh_50 to 54` +
       `all_ages_malik_et_al_survey_vax_uptake_55 to 59` * `sh_55 to 59`) / (`sh_50 to 54` + `sh_55 to 59`)
    ),
    `agebin_60-69__vax_uptake_malik_all_ages` = `all_ages_malik_et_al_survey_vax_uptake_60 to 64`,
    `agebin_70-79__vax_uptake_malik_all_ages` = `all_ages_malik_et_al_survey_vax_uptake_70 to 74`,
    `agebin_80+__vax_uptake_malik_all_ages`   = `all_ages_malik_et_al_survey_vax_uptake_75+`,
    # Malik et al. (2020) 0 uptake for 0-15 year olds
    `agebin_0-9__vax_uptake_malik`            = `malik_et_al_survey_vax_uptake_0 to 4`,
    `agebin_10-19__vax_uptake_malik` = (
      (`malik_et_al_survey_vax_uptake_10 to 15` * `sh_10 to 15` +
       `malik_et_al_survey_vax_uptake_16 to 19` * `sh_16 to 19`) / (`sh_10 to 15` + `sh_16 to 19`)
    ),
    `agebin_20-29__vax_uptake_malik` = (
      (`malik_et_al_survey_vax_uptake_20 to 24` * `sh_20 to 24` +
       `malik_et_al_survey_vax_uptake_25 to 29` * `sh_25 to 29`) / (`sh_20 to 24` + `sh_25 to 29`) 
    ),
    `agebin_30-39__vax_uptake_malik` = (
      (`malik_et_al_survey_vax_uptake_30 to 34` * `sh_30 to 34` +
       `malik_et_al_survey_vax_uptake_35 to 39` * `sh_35 to 39`) / (`sh_30 to 34` + `sh_35 to 39`)
    ),
    `agebin_40-49__vax_uptake_malik` = (
      (`malik_et_al_survey_vax_uptake_40 to 44` * `sh_40 to 44` +
       `malik_et_al_survey_vax_uptake_45 to 49` * `sh_45 to 49`) / (`sh_40 to 44` + `sh_45 to 49`)
    ),
    `agebin_50-59__vax_uptake_malik` = (
      (`malik_et_al_survey_vax_uptake_50 to 54` * `sh_50 to 54` +
       `malik_et_al_survey_vax_uptake_55 to 59` * `sh_55 to 59`) / (`sh_50 to 54` + `sh_55 to 59`)
    ),
    `agebin_60-69__vax_uptake_malik`           = `malik_et_al_survey_vax_uptake_60 to 64`,
    `agebin_70-79__vax_uptake_malik`           = `malik_et_al_survey_vax_uptake_70 to 74`,
    `agebin_80+__vax_uptake_malik`             = `malik_et_al_survey_vax_uptake_75+`,
    # US Census Household Survey Week 22 All Ages assumption
    `agebin_0-9__vax_uptake_census_all_ages`   = `all_ages_census_uptake_survey_0 to 4`,
    `agebin_10-19__vax_uptake_census_all_ages` = `all_ages_census_uptake_survey_10 to 15`,
    `agebin_20-29__vax_uptake_census_all_ages` = (
      (`all_ages_census_uptake_survey_20 to 24` * `sh_20 to 24` +
       `all_ages_census_uptake_survey_25 to 29` * `sh_25 to 29`) / (`sh_20 to 24` + `sh_25 to 29`)
    ),
    `agebin_30-39__vax_uptake_census_all_ages` = `all_ages_census_uptake_survey_30 to 34`,
    `agebin_40-49__vax_uptake_census_all_ages` = `all_ages_census_uptake_survey_40 to 44`,
    `agebin_50-59__vax_uptake_census_all_ages` = (
      (`all_ages_census_uptake_survey_50 to 54` * `sh_50 to 54` +
       `all_ages_census_uptake_survey_55 to 59` * `sh_55 to 59`) / (`sh_50 to 54` + `sh_55 to 59`)
    ),
    `agebin_60-69__vax_uptake_census_all_ages` = (
      (`all_ages_census_uptake_survey_60 to 64` * `sh_60 to 64` +
       `all_ages_census_uptake_survey_65 to 69` * `sh_65 to 69`) / (`sh_60 to 64` + `sh_65 to 69`)
    ),
    `agebin_70-79__vax_uptake_census_all_ages` = `all_ages_census_uptake_survey_70 to 74`,
    `agebin_80+__vax_uptake_census_all_ages`   = `all_ages_census_uptake_survey_75+`,
    # US Census Household Survey Week 22 0 uptake for 0-15 year olds
    `agebin_0-9__vax_uptake_census`            = `census_uptake_survey_0 to 4`,
    `agebin_10-19__vax_uptake_census` = (
      (`census_uptake_survey_10 to 15` * `sh_10 to 15` +
       `census_uptake_survey_16 to 19` * `sh_16 to 19`) / (`sh_10 to 15` + `sh_16 to 19`)
    ),
    `agebin_20-29__vax_uptake_census` = (
      (`census_uptake_survey_20 to 24` * `sh_20 to 24` +
       `census_uptake_survey_25 to 29` * `sh_25 to 29`) / (`sh_20 to 24` + `sh_25 to 29`)
    ),
    `agebin_30-39__vax_uptake_census`          = `census_uptake_survey_30 to 34`,
    `agebin_40-49__vax_uptake_census`          = `census_uptake_survey_40 to 44`,
    `agebin_50-59__vax_uptake_census` = (
      (`census_uptake_survey_50 to 54` * `sh_50 to 54` +
       `census_uptake_survey_55 to 59` * `sh_55 to 59`) / (`sh_50 to 54` + `sh_55 to 59`)
    ),
    `agebin_60-69__vax_uptake_census` = (
      (`census_uptake_survey_60 to 64` * `sh_60 to 64` +
       `census_uptake_survey_65 to 69` * `sh_65 to 69`) / (`sh_60 to 64` + `sh_65 to 69`)
    ),
    `agebin_70-79__vax_uptake_census`          = `census_uptake_survey_70 to 74`,
    `agebin_80+__vax_uptake_census`            = `census_uptake_survey_75+`,
  ) %>%
  select(
    starts_with("agebin")
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("agebin"),
    names_to = c("group", ".value"),
    names_sep = "__"
  ) 

# *** Vaccine Efficacy ****************************************************************
# Manually compute Pfizer efficacy for age-groups not explicitly shown
# See Table 3 (https://www.nejm.org/doi/full/10.1056/NEJMoa2034577) for Pfizer trial data and
# CDC (https://www.cdc.gov/csels/dsepd/ss1978/lesson3/section6.html) for vax efficacy equation
pfizer_56pl_placebo_positives            <- 48
pfizer_56pl_treatment_positives          <- 3
pfizer_65pl_placebo_positives            <- 19
pfizer_65pl_treatment_positives          <- 1
pfizer_75pl_placebo_positives            <- 5
pfizer_75pl_treatment_positives          <- 0 
implied_pfizer_56_64_placebo_positives   <- (pfizer_56pl_placebo_positives - pfizer_65pl_placebo_positives)
implied_pfizer_56_64_treatment_positives <- (pfizer_56pl_treatment_positives - pfizer_65pl_treatment_positives)
implied_pfizer_65_74_placebo_positives   <- (pfizer_65pl_placebo_positives - pfizer_75pl_placebo_positives)
implied_pfizer_65_74_treatment_positives <- (pfizer_65pl_treatment_positives - pfizer_75pl_treatment_positives)
pfizer_56_64_vax_efficacy                <- (
  (implied_pfizer_56_64_placebo_positives - implied_pfizer_56_64_treatment_positives)
  / implied_pfizer_56_64_placebo_positives
)
pfizer_65_74_vax_efficacy                <- (
  (implied_pfizer_65_74_placebo_positives - implied_pfizer_65_74_treatment_positives) 
  / implied_pfizer_65_74_placebo_positives
)

vax_efficacy <- raw_pums %>%
  select(
    PWGTP,
    AGEP
  ) %>%
  mutate(
    fine_age_bin = case_when(
      AGEP <= 4              ~ "0_4",
      AGEP > 4  & AGEP <= 9  ~ "5_9",
      AGEP > 9  & AGEP <= 14 ~ "10_14",
      AGEP > 14 & AGEP <= 19 ~ "15_19",
      AGEP > 19 & AGEP <= 24 ~ "20_24",
      AGEP > 24 & AGEP <= 29 ~ "25_29",
      AGEP > 29 & AGEP <= 34 ~ "30_34",
      AGEP > 34 & AGEP <= 39 ~ "35_39",
      AGEP > 39 & AGEP <= 44 ~ "40_44",
      AGEP > 44 & AGEP <= 49 ~ "45_49",
      AGEP > 49 & AGEP <= 55 ~ "50_55",
      AGEP > 55 & AGEP <= 59 ~ "56_59",
      AGEP > 59 & AGEP <= 64 ~ "60_64",
      AGEP > 64 & AGEP <= 69 ~ "65_69",
      AGEP > 69 & AGEP <= 74 ~ "70_74",
      AGEP > 74 & AGEP <= 79 ~ "75_79",
      AGEP > 79              ~ "80+",
      TRUE ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(fine_age_bin)
  ) %>%
  select(-AGEP) %>%
  group_by(fine_age_bin) %>%
  summarize(
    n = sum(PWGTP)
  ) %>%
  mutate(
    sh = n / sum(n)
  ) %>%
  mutate(
    # Moderna
    # Figure 4 of NEJM Baden et al. (2020) - efficacy at preventing COVID-19 (assumed to be symptom onset)
    # https://www.nejm.org/doi/full/10.1056/NEJMoa2035389
    moderna_efficacy = case_when(
      # For younger ages, assumed to equal 18-65 efficacy
      fine_age_bin == "0_4"     ~ .956,
      fine_age_bin == "5_9"     ~ .956,
      fine_age_bin == "10_14"   ~ .956,
      fine_age_bin == "15_19"   ~ .956,
      fine_age_bin == "20_24"   ~ .956,
      fine_age_bin == "25_29"   ~ .956,
      fine_age_bin == "30_34"   ~ .956,
      fine_age_bin == "35_39"   ~ .956,
      fine_age_bin == "40_44"   ~ .956,
      fine_age_bin == "45_49"   ~ .956,
      fine_age_bin == "50_55"   ~ .956,
      fine_age_bin == "56_59"   ~ .956,
      fine_age_bin == "60_64"   ~ .956,
      fine_age_bin == "65_69"   ~ .864,
      fine_age_bin == "70_74"   ~ .864,
      fine_age_bin == "75_79"   ~ .864,
      fine_age_bin == "80+"     ~ .864,
      TRUE ~ NA_real_
    ),
    # Pfizer
    # Table 3 of NEJM Polack et al. (2020) - efficacy at preventing COVID-19 (assumed to be case positive)
    # https://www.nejm.org/doi/full/10.1056/NEJMoa2034577
    pfizer_efficacy = case_when(
      # For younger ages, assumed to equal 16-55 efficacy
      fine_age_bin == "0_4"     ~ .956,
      fine_age_bin == "5_9"     ~ .956,
      fine_age_bin == "10_14"   ~ .956,
      fine_age_bin == "15_19"   ~ .956,
      fine_age_bin == "20_24"   ~ .956,
      fine_age_bin == "25_29"   ~ .956,
      fine_age_bin == "30_34"   ~ .956,
      fine_age_bin == "35_39"   ~ .956,
      fine_age_bin == "40_44"   ~ .956,
      fine_age_bin == "45_49"   ~ .956,
      fine_age_bin == "50_55"   ~ .956,
      fine_age_bin == "56_59"   ~ pfizer_56_64_vax_efficacy,
      fine_age_bin == "60_64"   ~ pfizer_56_64_vax_efficacy,
      fine_age_bin == "65_69"   ~ pfizer_65_74_vax_efficacy,
      fine_age_bin == "70_74"   ~ pfizer_65_74_vax_efficacy,
      fine_age_bin == "75_79"   ~ 1,
      fine_age_bin == "80+"     ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  assertr::verify(
    ! is.na(pfizer_efficacy) & ! is.na(moderna_efficacy)
  ) %>%
  # Convert to SEIR age-groups
  mutate(
    filler_col = "Filler"
  ) %>%
  tidyr::pivot_wider(
    id_cols     = filler_col,
    names_from  = fine_age_bin,
    values_from = c(n, sh, moderna_efficacy, pfizer_efficacy),
    names_glue  = "{.value}_{fine_age_bin}"
  ) %>%
  mutate(
    # Align Moderna
    `agebin_0-9__moderna_efficacy`   = `moderna_efficacy_0_4`,
    `agebin_10-19__moderna_efficacy` = `moderna_efficacy_10_14`,
    `agebin_20-29__moderna_efficacy` = `moderna_efficacy_20_24`,
    `agebin_30-39__moderna_efficacy` = `moderna_efficacy_30_34`,
    `agebin_40-49__moderna_efficacy` = `moderna_efficacy_40_44`,
    `agebin_50-59__moderna_efficacy` = `moderna_efficacy_50_55`,
    `agebin_60-69__moderna_efficacy` = (
      (`moderna_efficacy_60_64` * `sh_60_64` +
       `moderna_efficacy_65_69` * `sh_65_69`) / (`sh_60_64` + `sh_65_69`)
    ),
    `agebin_70-79__moderna_efficacy` = `moderna_efficacy_70_74`,
    `agebin_80+__moderna_efficacy`   = `moderna_efficacy_80+`,
    # Align Pfizer
    `agebin_0-9__pfizer_efficacy`    = `pfizer_efficacy_0_4`,
    `agebin_10-19__pfizer_efficacy`  = `pfizer_efficacy_10_14`,
    `agebin_20-29__pfizer_efficacy`  = `pfizer_efficacy_20_24`,
    `agebin_30-39__pfizer_efficacy`  = `pfizer_efficacy_30_34`,
    `agebin_40-49__pfizer_efficacy`  = `pfizer_efficacy_40_44`,
    `agebin_50-59__pfizer_efficacy`  = (
      (`pfizer_efficacy_50_55` * `sh_50_55` +
       `pfizer_efficacy_56_59` * `sh_56_59`) / (`sh_50_55` + `sh_56_59`)
    ),
    `agebin_60-69__pfizer_efficacy`  = (
      (`pfizer_efficacy_60_64` * `sh_60_64` +
       `pfizer_efficacy_65_69` * `sh_65_69`) / (`sh_60_64` + `sh_65_69`)
    ),
    `agebin_70-79__pfizer_efficacy`  = (
      (`pfizer_efficacy_70_74` * `sh_70_74` +
       `pfizer_efficacy_75_79` * `sh_75_79`) / (`sh_70_74` + `sh_75_79`)
    ),
    `agebin_80+__pfizer_efficacy`    = `pfizer_efficacy_80+`,
  ) %>%
  select(
    starts_with("agebin")
  ) %>%
  tidyr::pivot_longer(
    cols      = starts_with("agebin"),
    names_to  = c("group", ".value"),
    names_sep = "__"
  ) %>%
  mutate(
    average_2vax_efficacy = (moderna_efficacy + pfizer_efficacy) / 2 
  )

# *** Years Life Lost *****************************************************************
nchs_2017_yll <- dplyr::as_tibble(openxlsx::read.xlsx(
  stringr::str_c(
    raw_data_dir,
    "/nchs_years_life_lost/Table01.xlsx"
  ), sheet = "Table 1"
)) %>%
select(
  `Table.1..Life.table.for.the.total.population:.United.States,.2017`,
  X7
) %>%
rename(
  age_bin = `Table.1..Life.table.for.the.total.population:.United.States,.2017`,
  ex = X7
) %>%
filter(
  ! is.na(age_bin) & age_bin != "Age (years)" & age_bin != "SOURCE: NCHS, National Vital Statistics System, Mortality."
) %>%
mutate(
  lb_yr_range = case_when(
    age_bin == "100 and over" ~ "100",
    stringr::str_length(age_bin) < 5 ~ substr(age_bin, 1, 1),
    stringr::str_length(age_bin) < 7 ~ substr(age_bin, 1, 2),
    TRUE ~ NA_character_
  ),
  ub_yr_range = case_when(
    age_bin == "100 and over" ~ NA_character_,
    stringr::str_length(age_bin) < 5 ~ trimws(substr(age_bin, 3, stringr::str_length(age_bin))),
    stringr::str_length(age_bin) < 7 ~ trimws(substr(age_bin, 4, stringr::str_length(age_bin))),
    TRUE ~ NA_character_
  )
) %>%
assertr::verify(
  ! is.na(lb_yr_range) & (! is.na(ub_yr_range) | age_bin == "100 and over")
) %>%
select(
  -age_bin
) %>%
mutate(
  lb_yr_range = as.numeric(lb_yr_range),
  ub_yr_range = as.numeric(ub_yr_range),
  ex = as.numeric(ex)
) %>%
rename(
  AGEP = lb_yr_range
) %>%
select(
  -ub_yr_range
) %>%
dplyr::right_join(
  select(raw_pums, PWGTP, AGEP),
  by = c("AGEP")
) %>%
assertr::verify(
  ! is.na(ex)
) %>%
mutate(
    group = case_when(
      AGEP < 10              ~ "agebin_0-9",
      AGEP < 20 & AGEP >= 10 ~ "agebin_10-19",
      AGEP < 30 & AGEP >= 20 ~ "agebin_20-29",
      AGEP < 40 & AGEP >= 30 ~ "agebin_30-39",
      AGEP < 50 & AGEP >= 40 ~ "agebin_40-49",
      AGEP < 60 & AGEP >= 50 ~ "agebin_50-59",
      AGEP < 70 & AGEP >= 60 ~ "agebin_60-69",
      AGEP < 80 & AGEP >= 70 ~ "agebin_70-79",
      AGEP >= 80             ~ "agebin_80+",
      TRUE                   ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(group)
  ) %>%
  select(-AGEP) %>%
  group_by(group) %>%
  summarize(
    n = sum(PWGTP),
    ex = weighted.mean(ex, PWGTP)
  ) %>%
  ungroup() %>%
  rename(
    yll = ex
  ) %>%
  select(
    -n
  )

# *** Create Master Table *************************************************************
master_group_table <- pop_sh_by_age_bin %>%
  dplyr::left_join(
    ifr_by_age_bin,
    by = c("group")
  ) %>%
  dplyr::left_join(
    sus_to_inf_by_age_bin,
    by = c("group")
  ) %>%
  dplyr::left_join(
    vax_uptake_by_age_bin,
    by = c("group")
  ) %>%
  dplyr::left_join(
    vax_efficacy,
    by = c("group")
  ) %>%
  dplyr::left_join(
    nchs_2017_yll,
    by = c("group")
  ) %>%
  mutate(
    entity = "USA (50 states + DC)"
  ) %>%
  relocate(
    entity, nat_population, group, n_group, c19_ifr_group, sh_group,
    sus_to_inf, clinical_sus, vax_uptake_malik, vax_uptake_census,
    moderna_efficacy, pfizer_efficacy, average_2vax_efficacy, vax_uptake_malik_all_ages,
    vax_uptake_census_all_ages, yll
  )

write.csv(
  master_group_table,
  stringr::str_c(
    build_data_dir,
    "/national_baseline.csv"
  ),
  row.names = FALSE
)