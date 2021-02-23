# *** Header **************************************************************************
#
# Create S-E-I-R initial conditions as of December 14, 2020
#

# *** Read in 5-Year ACS PUMS Data ****************************************************
raw_pums <- dplyr::as_tibble(data.table::fread(
  file = stringr::str_c(
    build_data_dir,
    "/pums_data_with_nasem_phases.csv"
  ),
  sep = ","
))

# Population share by fine age bin
pop_shares_by_fine_age_bin <- raw_pums %>%
  select(
    PWGTP,
    AGEP
  ) %>%
  mutate(
    fine_age_bin = case_when(
      AGEP <  5              ~ "0_4",
      AGEP >= 5  & AGEP < 10 ~ "5_9",
      AGEP >= 10 & AGEP < 15 ~ "10_14",
      AGEP >= 15 & AGEP < 20 ~ "15_19",
      AGEP >= 20 & AGEP < 25 ~ "20_24",
      AGEP >= 25 & AGEP < 30 ~ "25_29",
      AGEP >= 30 & AGEP < 35 ~ "30_34",
      AGEP >= 35 & AGEP < 40 ~ "35_39",
      AGEP >= 40 & AGEP < 45 ~ "40_44",
      AGEP >= 45 & AGEP < 50 ~ "45_49",
      AGEP >= 50 & AGEP < 55 ~ "50_54",
      AGEP >= 55 & AGEP < 60 ~ "55_59",
      AGEP >= 60 & AGEP < 65 ~ "60_64",
      AGEP >= 65 & AGEP < 70 ~ "65_69",
      AGEP >= 70 & AGEP < 75 ~ "70_74",
      AGEP >= 75 & AGEP < 80 ~ "75_79",
      AGEP >= 80 & AGEP < 85 ~ "80_84",
      AGEP >= 85             ~ "85+",
      TRUE ~ NA_character_
    )
   ) %>%
   assertr::verify(
     ! is.na(fine_age_bin)
   ) %>%
   group_by(fine_age_bin) %>%
   summarize(
     n = sum(PWGTP)
   ) %>%
   mutate(
     sh = n / sum(n),
     entity = "US"
   ) %>%
   ungroup() %>%
   select(
     -n
   ) %>%
   tidyr::pivot_wider(
     id_cols     = entity,
     names_from  = fine_age_bin,
     values_from = sh,
     names_glue  = "{.value}_{fine_age_bin}"
   ) %>%
   assertr::verify(
     round(1 - (
       sh_0_4 +
       sh_5_9 +
       sh_10_14 +
       sh_15_19 +
       sh_20_24 +
       sh_25_29 +
       sh_30_34 +
       sh_35_39 +
       sh_40_44 +
       sh_45_49 +
       sh_50_54 +
       sh_55_59 +
       sh_60_64 +
       sh_65_69 +
       sh_70_74 +
       sh_75_79 +
       sh_80_84 +
       `sh_85+`
       ), 5
     ) == 0
   )

# *** Clean CDC Deaths by Age by Week *************************************************
raw_cdc_by_age_group_by_week <- dplyr::as_tibble(
  read.csv(
    stringr::str_c(
      raw_data_dir,
      "/cdc_death_count/21-02-07/Provisional_COVID-19_Death_Counts_by_Sex__Age__and_Week.csv"
    )
  )
)

clean_cdc_age_over_time <- raw_cdc_by_age_group_by_week %>%
  select(
    -'Data.as.of', -'State', -'Total.Deaths'
  ) %>%
  rename(
    week_id          = MMWR.Week,
    sex              = Sex,
    age_bin          = Age.Group,
    end_of_week_date = End.Week,
    incr_deaths      = COVID.19.Deaths
  ) %>%
  filter(
    age_bin != "All Ages",
    sex == "All Sex"
  ) %>%
  select(
    -sex, -week_id
  ) %>%
  mutate(
    age_bin = case_when(
      age_bin == "Under 1 year"      ~ "-0",
      age_bin == "1-4 Years"         ~ "1_4",
      age_bin == "5-14 Years"        ~ "5_14",
      age_bin == "15-24 Years"       ~ "15_24",
      age_bin == "25-34 Years"       ~ "25_34",
      age_bin == "35-44 Years"       ~ "35_44",
      age_bin == "45-54 Years"       ~ "45_54",
      age_bin == "55-64 Years"       ~ "55_64",
      age_bin == "65-74 Years"       ~ "65_74",
      age_bin == "75-84 Years"       ~ "75_84",
      age_bin == "85 Years and Over" ~ "85+",
      TRUE ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(age_bin)
  ) %>%
  assertr::verify(
    stringr::str_length(end_of_week_date) == 10
  ) %>%
  mutate(
    month = substr(end_of_week_date, 1, 2),
    day   = substr(end_of_week_date, 4, 5),
    year  = substr(end_of_week_date, 7, 10),
    end_of_week_date = stringr::str_c(year, "-", month, "-", day)
  ) %>%
  select(
    -month, -day, -year
  ) %>%
  filter(
    # Last week needed
    end_of_week_date <= "2021-01-09"
  ) %>%
  mutate(
    entity = "US"
  ) %>%
  tidyr::pivot_wider(
    id_cols     = c(end_of_week_date, entity),
    names_from  = age_bin,
    values_from = incr_deaths,
    names_glue  = "incr_deaths__{age_bin}"
  ) %>%
  mutate(
    incr_deaths__0_4 = (`incr_deaths__-0` + incr_deaths__1_4)
  ) %>%
  select(
    -`incr_deaths__-0`, -incr_deaths__1_4
  ) %>%
  left_join(
    pop_shares_by_fine_age_bin,
    by = c("entity")
  ) %>%
  # Re-align age-bins to SEIR groups using PUMS population shares
  mutate(
    incr_deaths__0_9 = (
      (
        incr_deaths__0_4
      ) +
      (
        incr_deaths__5_14 * (sh_5_9 / (sh_5_9 + sh_10_14))
      )
    ),
    incr_deaths__10_19 = (
      (
        incr_deaths__5_14 * (sh_10_14 / (sh_5_9 + sh_10_14))
      ) +
      (
        incr_deaths__15_24 * (sh_15_19 / (sh_15_19 + sh_20_24))
      )
    ),
    incr_deaths__20_29 = (
      (
        incr_deaths__15_24 * (sh_20_24 / (sh_15_19 + sh_20_24))
      ) +
      (
        incr_deaths__25_34 * (sh_25_29 / (sh_25_29 + sh_30_34))
      )
    ),
    incr_deaths__30_39 = (
      (
        incr_deaths__25_34 * (sh_30_34 / (sh_25_29 + sh_30_34))
      ) +
      (
        incr_deaths__35_44 * (sh_35_39 / (sh_35_39 + sh_40_44))
      )
    ),
    incr_deaths__40_49 = (
      (
        incr_deaths__35_44 * (sh_40_44 / (sh_35_39 + sh_40_44))
      ) +
      (
        incr_deaths__45_54 * (sh_45_49 / (sh_45_49 + sh_50_54))
      )
    ),
    incr_deaths__50_59 = (
      (
        incr_deaths__45_54 * (sh_50_54 / (sh_45_49 + sh_50_54))
      ) +
      (
        incr_deaths__55_64 * (sh_55_59 / (sh_55_59 + sh_60_64))
      )
    ),
    incr_deaths__60_69 = (
      (
        incr_deaths__55_64 * (sh_60_64 / (sh_55_59 + sh_60_64))
      ) +
      (
        incr_deaths__65_74 * (sh_65_69 / (sh_65_69 + sh_70_74))
      )
    ),
    incr_deaths__70_79 = (
      (
        incr_deaths__65_74 * (sh_70_74 / (sh_65_69 + sh_70_74))
      ) +
      (
        incr_deaths__75_84 * (sh_75_79 / (sh_75_79 + sh_80_84))
      )
    ),
    `incr_deaths__80+` = (
      (
        incr_deaths__75_84 * (sh_80_84 / (sh_75_79 + sh_80_84))
      ) +
      (
        `incr_deaths__85+`
      )
    )
  ) %>%
  assertr::verify(
    round(
      (
        incr_deaths__0_4 +
        incr_deaths__5_14 +
        incr_deaths__15_24 +
        incr_deaths__25_34 +
        incr_deaths__35_44 +
        incr_deaths__45_54 +
        incr_deaths__55_64 +
        incr_deaths__65_74 +
        incr_deaths__75_84 +
        `incr_deaths__85+`
      ) - (
        incr_deaths__0_9 +
        incr_deaths__10_19 +
        incr_deaths__20_29 +
        incr_deaths__30_39 +
        incr_deaths__40_49 +
        incr_deaths__50_59 +
        incr_deaths__60_69 +
        incr_deaths__70_79 +
        `incr_deaths__80+`
      ), 5
    ) == 0
  ) %>%
  select(
    end_of_week_date,
    incr_deaths__0_9, incr_deaths__10_19, incr_deaths__20_29, incr_deaths__30_39,
    incr_deaths__40_49, incr_deaths__50_59, incr_deaths__60_69, incr_deaths__70_79,
    `incr_deaths__80+`
  ) %>%
  tidyr::pivot_longer(
    cols      = -c(end_of_week_date),
    names_to  = c(".value", "group"),
    names_sep = "__"
  ) %>%
  group_by(end_of_week_date) %>%
  mutate(
    sh_incr_deaths = case_when(
      sum(incr_deaths) == 0 ~ 0,
      TRUE ~ incr_deaths / sum(incr_deaths)
    )
  ) %>%
  ungroup() %>%
  select(
    -incr_deaths
  )

# Day -> week end date x-walk
day_to_end_of_week_date_xwalk <- clean_cdc_age_over_time %>%
  select(
    end_of_week_date
  ) %>%
  distinct(
    end_of_week_date
  ) %>%
  mutate(
    day_1 = as.character(as.Date(end_of_week_date) - 6),
    day_2 = as.character(as.Date(end_of_week_date) - 5),
    day_3 = as.character(as.Date(end_of_week_date) - 4),
    day_4 = as.character(as.Date(end_of_week_date) - 3),
    day_5 = as.character(as.Date(end_of_week_date) - 2),
    day_6 = as.character(as.Date(end_of_week_date) - 1),
    day_7 = as.character(as.Date(end_of_week_date) - 0),
  ) %>%
  tidyr::pivot_longer(
    cols       = -c(end_of_week_date),
    names_to   = c(".value", "day_num"),
    names_sep  = "_"
  ) %>%
  select(
    -day_num
  ) %>%
  rename(
    date = day
  )

# *** Clean New York Times Deaths by State by Day *************************************
raw_nyt_case_death_count_st_lvl <- read.csv(
  file = stringr::str_c(
    raw_data_dir,
    "/nyt_death_infection_count/",
    "21-01-27",
    "/us-states.csv"
  )
)
clean_nyt <- raw_nyt_case_death_count_st_lvl %>%
  filter(
    # 50 states + DC
    ! state %in% c(
      "Guam",
      "Northern Mariana Islands",
      "Puerto Rico",
      "Virgin Islands"
    ) 
  ) %>%
  group_by(date) %>%
  summarize(
    tot_cases = sum(cases),
    tot_death = sum(deaths)
  ) %>%
  arrange(date) %>%
  mutate(
    incr_cases = tot_cases - lag(tot_cases),
    incr_cases = case_when(
      is.na(incr_cases) ~ tot_cases,
      TRUE ~ incr_cases
    ),
    incr_deaths = tot_death - lag(tot_death),
    incr_deaths = case_when(
      is.na(incr_deaths) ~ tot_death,
      TRUE ~ incr_deaths
    )
  ) %>%
  filter(
    # 24 days after 12/14/20
    date <= "2021-01-07"
  ) %>%
  dplyr::left_join(
    day_to_end_of_week_date_xwalk,
    by = c("date")
  ) %>%
  assertr::verify(
    ! is.na(end_of_week_date)
  ) %>%
  mutate(
    # Smooth via 7 day rolling average
    smoothed_incr_deaths = zoo::rollapplyr(incr_deaths, 7, mean, partial = TRUE)
  )

clean_nyt_by_age_group <- clean_nyt %>%
  select(
    -tot_cases, -tot_death, -incr_cases, -incr_deaths
  ) %>%
  full_join(
    clean_cdc_age_over_time,
    by = c("end_of_week_date")
  ) %>%
  filter(
    # Early pandemic days w/o daily death records
    ! is.na(smoothed_incr_deaths)
  ) %>%
  mutate(
    smoothed_incr_deaths = smoothed_incr_deaths * sh_incr_deaths
  ) %>%
  select(
    -end_of_week_date, -sh_incr_deaths
  ) %>%
  relocate(
    date, group, smoothed_incr_deaths
  ) %>%
  arrange(
    date, group
  )

# *** Compute New True Infections by Day **********************************************
# Assume discretely uniform distributed infection-to-death lag on 18 to 24 day interval
lag_floor      <- 18
lag_ceil       <- 24
lag_period_len <- (lag_ceil - lag_floor) + 1

# Assert the NYT is not missing any days in their data
stopifnot(
  length(unique(clean_nyt_by_age_group$date)) ==
  as.Date(unique(clean_nyt_by_age_group$date)[length(unique(clean_nyt_by_age_group$date))]) - as.Date(unique(clean_nyt_by_age_group$date[1])) + 1
)

nyt_clean_true_infections <- clean_nyt_by_age_group %>%
  mutate(
    # Levin et al. (2020) IFRs (Supplementary Appendix Q)
    # https://static-content.springer.com/esm/art%3A10.1007%2Fs10654-020-00698-1/MediaObjects/10654_2020_698_MOESM1_ESM.pdf
    levin_ifr = case_when(
      group == "0_9"   ~ .00001,
      group == "10_19" ~ .00003,
      group == "20_29" ~ .00011,
      group == "30_39" ~ .00037,
      group == "40_49" ~ .00123,
      group == "50_59" ~ .00413,
      group == "60_69" ~ .01380,
      group == "70_79" ~ .04620,
      group == "80+"   ~ .15460,
      TRUE             ~ NA_real_
    )
  ) %>%
  assertr::verify(
    ! is.na(levin_ifr)
  ) %>%
  mutate(
    naive_levin_infections       = smoothed_incr_deaths / levin_ifr,
    distributed_levin_infections = naive_levin_infections / lag_period_len,
  ) %>%
  arrange(
    group, desc(date)
  ) %>%
  group_by(group) %>%
  mutate(
    # Docs on rollapplyr (https://stackoverflow.com/questions/30153835/r-dplyr-rolling-sum)
    incr_levin_infections = (
      zoo::rollapplyr(distributed_levin_infections, 25, sum, partial = TRUE) - 
      zoo::rollapplyr(distributed_levin_infections, 18, sum, partial = TRUE)
    ),
  ) %>%
  ungroup() %>%
  arrange(group, date) %>%
  group_by(group) %>%
  mutate(
    # 4 day latent period + 9 day infectious period
    # Assumes infections take place at EOD. So a new infection
    # today is not included in the latent period because it kicks
    # in at midnight. However, a new infection four days ago is in its
    # fourth day of latency and is still considered latent today.
    infectious_levin_infections = (
      lag(incr_levin_infections, 13, default = 0) +
      lag(incr_levin_infections, 12, default = 0) +
      lag(incr_levin_infections, 11, default = 0) +
      lag(incr_levin_infections, 10, default = 0) +
      lag(incr_levin_infections, 9, default = 0)  +
      lag(incr_levin_infections, 8, default = 0)  +
      lag(incr_levin_infections, 7, default = 0)  +
      lag(incr_levin_infections, 6, default = 0)  +
      lag(incr_levin_infections, 5, default = 0)
    ),
    latent_levin_infections = (
      lag(incr_levin_infections, 4, default = 0) +
      lag(incr_levin_infections, 3, default = 0) +
      lag(incr_levin_infections, 2, default = 0) +
      lag(incr_levin_infections, 1, default = 0)
    )
  ) %>%
  ungroup() %>%
  filter(
      date <= "2020-12-14"
  ) %>%
  select(
    -starts_with("naive"), -starts_with("distributed")
  )

# *** Format SEIR Inputs **************************************************************
national_baseline_file <- dplyr::as_tibble(data.table::fread(
  file = stringr::str_c(
    build_data_dir,
    "/national_baseline.csv"
  ),
  sep = ","
))

# *** Generate Supplementary Tables ***************************************************
supplementary_tables_df <- nyt_clean_true_infections %>%
  arrange(group, date) %>%
  group_by(group) %>%
  mutate(
    cumul_smoothed_deaths  = cumsum(smoothed_incr_deaths),
    cumul_levin_infections = cumsum(incr_levin_infections)
  ) %>%
  ungroup() %>%
  select(
    -levin_ifr
  ) %>%
  mutate(
    group = stringr::str_c("agebin_", stringr::str_replace(group, "_", "-"))
  ) %>%
  left_join(
    select(national_baseline_file, group, n_group),
    by = c("group")
  ) %>%
  select(
    date, group, cumul_levin_infections, n_group
  )

# Get benchmark values
icl <- data.table::fread(
    stringr::str_c(
      raw_data_dir,
      "/icl/2020-12-20_v6.csv"
    )
  ) %>%
  filter(
    compartment == "cumulative_infections",
    iso3c == "USA",
    scenario == "Maintain Status Quo",
    date == "2020-12-14"
  ) %>%
  select(
    y_median
  )
icl <- as.numeric(icl[1,1])

gu <- data.table::fread(
    stringr::str_c(
      raw_data_dir,
      "/covid19projections/latest_all_estimates_us.csv"
    )
  ) %>%
  select(
    date, total_infected_mean
  ) %>%
  filter(
    date == "2020-12-14"
  ) %>%
  select(
    -date
  )

gu <- as.numeric(gu[1,1])

table_s4 <- supplementary_tables_df %>%
  group_by(
    date
  ) %>%
  summarize(
    our_imputation_cumul_infections = sum(cumul_levin_infections),
    n = sum(n_group)
  ) %>% 
  ungroup() %>% 
  filter(
    date %in% c("2020-12-14", "2020-11-15")
  ) %>%
  mutate(
    share_of_total_population_infected = round(our_imputation_cumul_infections / n, 3)*100,
    our_equivalent_imputation = case_when(
      date == "2020-11-15" ~ round(our_imputation_cumul_infections, 0),
      TRUE ~ share_of_total_population_infected
    ),
    benchmark = case_when(
      date == "2020-11-15" ~ "Angulo (2021)",
      TRUE ~ "Noh (2021)"
    )
  ) %>%
  add_row(
    date = "2020-12-14",
    our_imputation_cumul_infections = NA_integer_,
    n = NA_integer_,
    share_of_total_population_infected = NA_real_,
    our_equivalent_imputation = sum(.$our_imputation_cumul_infections[.$benchmark == "Noh (2021)"]),
    benchmark = "Gu (2021)"
  ) %>%
  add_row(
    date = "2020-12-14",
    our_imputation_cumul_infections = NA_integer_,
    n = NA_integer_,
    share_of_total_population_infected = NA_real_,
    our_equivalent_imputation = sum(.$our_imputation_cumul_infections[.$benchmark == "Noh (2021)"]),
    benchmark = "Walker (2020)"
  ) %>%
  mutate(
    our_equivalent_imputation = case_when(
      benchmark == "Noh (2021)" ~ share_of_total_population_infected,
      TRUE ~ our_equivalent_imputation
    )
  ) %>%
  select(
    date, our_equivalent_imputation, benchmark
  ) %>%
  mutate(
    benchmark_val = case_when(
      # (Accessed Feburary 20, 2021) - https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2774584
      benchmark == "Angulo (2021)" ~ 46910006.0,
      # (Accessed February 20, 2021) - https://github.com/JungsikNoh/COVID19_Estimated-Size-of-Infectious-Population/blob/main/output/countries/US/2020-12-14/US_estCumIncidence.png
      benchmark == "Noh (2021)" ~ 15.8,
      # (Accessed February 16, 2021) - https://github.com/mrc-ide/global-lmic-reports/tree/master/data
      benchmark == "Walker (2020)" ~ icl,
      # (Accessed February 16, 2021) - https://github.com/youyanggu/covid19-infection-estimates-latest
      benchmark == "Gu (2021)" ~ gu,
      TRUE ~ NA_real_
    )
  ) %>%
  assertr::verify(
    ! is.na(benchmark_val)
  ) %>%
  relocate(
    benchmark, date, benchmark_val, our_equivalent_imputation
  ) %>%
  rename(
    `Date as of` = date
  )

# Write the table
write.csv(
  table_s4,
  stringr::str_c(
    exhibit_data_dir,
    "/table_s4.csv"
  ),
  row.names = FALSE
)

table_s5 <- supplementary_tables_df %>%
  filter(
    date == "2020-12-14"
  ) %>%
  mutate(
    our_imputation_age_group = case_when(
      group == "agebin_0-9"   ~ "0-19",
      group == "agebin_10-19" ~ "0-19",
      group == "agebin_20-29" ~ "20-49",
      group == "agebin_30-39" ~ "20-49",
      group == "agebin_40-49" ~ "20-49",
      group == "agebin_50-59" ~ "50+",
      group == "agebin_60-69" ~ "50+",
      group == "agebin_70-79" ~ "50+",
      group == "agebin_80+"   ~ "50+",
      TRUE ~ NA_character_
    )
  ) %>%
  assertr::verify(
    ! is.na(our_imputation_age_group)
  ) %>%
  group_by(our_imputation_age_group) %>%
  summarize(
    our_imputation_cumul_infections_on_dec_14 = sum(cumul_levin_infections)
  ) %>%
  ungroup() %>%
  mutate(
    sh_of_our_imputation_cumul_infections_on_dec_14 = round(100*(our_imputation_cumul_infections_on_dec_14 / sum(our_imputation_cumul_infections_on_dec_14)), 1)
  ) %>%
  add_row(
    our_imputation_age_group                         = "All ages",
    our_imputation_cumul_infections_on_dec_14        = sum(.$our_imputation_cumul_infections_on_dec_14),
    sh_of_our_imputation_cumul_infections_on_dec_14  = 100
  ) %>%
  mutate(
    our_imputation_cumul_infections_on_dec_14 = prettyNum(our_imputation_cumul_infections_on_dec_14, big.mark = ","),
    cdc_imputation_age_group = case_when(
      our_imputation_age_group == "0-19" ~ "0-17",
      our_imputation_age_group == "20-49" ~ "18-49",
      TRUE ~ our_imputation_age_group
    ),
    # https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/burden.html (accessed February 20, 2021)
    cdc_imputation_cumul_infections_on_dec_31 = case_when(
      cdc_imputation_age_group == "0-17"     ~ (3001623 + 14550829),
      cdc_imputation_age_group == "18-49"    ~ (41940215),
      cdc_imputation_age_group == "50+"      ~ (14447134 + 9039683),
      cdc_imputation_age_group == "All ages" ~ (83111629),
    ),
    sh_of_cdc_imputation_cumul_infections_on_dec_31 = case_when(
      cdc_imputation_age_group == "All ages" ~ 100,
      TRUE ~ round(100*(cdc_imputation_cumul_infections_on_dec_31 / sum(cdc_imputation_cumul_infections_on_dec_31[cdc_imputation_age_group != "All ages"])), 1)
    ),
    cdc_imputation_cumul_infections_on_dec_31 = prettyNum(cdc_imputation_cumul_infections_on_dec_31, big.mark = ",")
  ) %>%
  relocate(
    cdc_imputation_age_group,
    cdc_imputation_cumul_infections_on_dec_31,
    sh_of_cdc_imputation_cumul_infections_on_dec_31,
    our_imputation_age_group,
    our_imputation_cumul_infections_on_dec_14,
    sh_of_our_imputation_cumul_infections_on_dec_14
  )

# Write the table
write.csv(
  table_s5,
  stringr::str_c(
    exhibit_data_dir,
    "/table_s5.csv"
  ),
  row.names = FALSE
)

seir_inputs_for_dec14_national_model <- nyt_clean_true_infections %>%
  arrange(group, date) %>%
  group_by(group) %>%
  mutate(
    cumul_smoothed_deaths  = cumsum(smoothed_incr_deaths),
    cumul_levin_infections = cumsum(incr_levin_infections)
  ) %>%
  ungroup() %>%
  select(
    -levin_ifr
  ) %>%
  mutate(
    group = stringr::str_c("agebin_", stringr::str_replace(group, "_", "-"))
  ) %>%
  left_join(
    select(national_baseline_file, group, n_group),
    by = c("group")
  ) %>%
  filter(
    date == "2020-12-14"
  ) %>%
  mutate(
    s0 = pmax(0, n_group - cumul_levin_infections) / n_group,
    e0 = latent_levin_infections / n_group,
    i0 = infectious_levin_infections / n_group,
    d0 = cumul_smoothed_deaths / n_group,
    r0 = (cumul_levin_infections - infectious_levin_infections - latent_levin_infections) / n_group
  ) %>%
  # Checks
  assertr::verify(
    (round(s0 + e0 + i0 + r0, 5) == 1 &
    ! is.na(s0) &
    ! is.na(e0) &
    ! is.na(i0) &
    ! is.na(d0) &
    ! is.na(r0) &
    r0 >= d0)
  ) %>%
  select(
    group,
    s0, e0, i0, d0, r0
  )

write.csv(
  seir_inputs_for_dec14_national_model,
  stringr::str_c(
    build_data_dir,
    "/seir_model_nat_lvl_infection_starting_points_by_group_on_dec14.csv"
  ),
  row.names = FALSE
)