# *** Header **************************************************************************
#
# This script assigns a NASEM phase to each person-level observation in the 5-year ACS
# PUMS data.

# ************************************************************************************* #
# -------------------------------- Read in PUMS Data ---------------------------------- #
# ************************************************************************************* #
pums_files <- list()
pums_dims  <- list()
file_keys <- c(
  "data_a",
  "data_b",
  "data_c",
  "data_d"
)
for (i in c("a", "b", "c", "d")) {
  person_lvl_pums <- data.table::fread(
    file = stringr::str_c(
      raw_data_dir,
      "/csv_pus/psam_pus",
      i,
      ".csv"
    ),
    select = c("SERIALNO", "PWGTP", "AGEP", "ST", "SEX", "HISP", "RAC1P", "MIL", "NAICSP", "SOCP", "SCHG")
  )
  person_lvl_pums$id <- stringr::str_c(rownames(person_lvl_pums), i)
  hh_lvl_pums <- data.table::fread(
    file = stringr::str_c(
      raw_data_dir,
      "/csv_hus/psam_hus",
      i,
      ".csv"
    ),
    select = c("SERIALNO", "MULTG", "TYPE")
  )
  
  pums_files[[stringr::str_c("data_", i)]] <- merge(
    person_lvl_pums,
    hh_lvl_pums,
    by = c("SERIALNO"),
    all.x = TRUE,
    all.y = FALSE
  )
  pums_dims[[stringr::str_c("data_", i)]] <- nrow(pums_files[[stringr::str_c("data_", i)]])
}
# Save raw
raw_pums_files <- pums_files


# ************************************************************************************* #
# ------------------------------- BRFSS Risk Factors ---------------------------------- #
# ************************************************************************************* #

# *** Clean BRFSS *********************************************************************
raw_BRFSS <- dplyr::as_tibble(SASxport::read.xport(
  stringr::str_c(
    raw_data_dir,
    '/brfss/LLCP2018.XPT'
  )
))

clean_BRFSS <- raw_BRFSS %>%
  mutate(
    num_covid_risk_factors = 0,
    one_of_key_vars_na = (
      is.na(CHCSCNCR)  | 
      is.na(CHCOCNCR)  | 
      is.na(CHCKDNY1)  | 
      is.na(CHCCOPD1)  | 
      is.na(X.BMI5CAT) | 
      is.na(CVDCRHD4)  | 
      is.na(DIABETE3)
    ),
    # Cancer
    num_covid_risk_factors = case_when(
      (CHCSCNCR == 1) | (CHCOCNCR == 1) ~ num_covid_risk_factors + 1,
      TRUE ~ num_covid_risk_factors
    ),
    # Chronic Kidney
    num_covid_risk_factors = case_when(
      (CHCKDNY1 == 1) ~ num_covid_risk_factors + 1,
      TRUE ~ num_covid_risk_factors
    ),
    # COPD
    num_covid_risk_factors = case_when(
      (CHCCOPD1 == 1) ~ num_covid_risk_factors + 1,
      TRUE ~ num_covid_risk_factors
    ),
    # BMI >= 30
    num_covid_risk_factors = case_when(
      (X.BMI5CAT == 4) ~ num_covid_risk_factors + 1,
      TRUE ~ num_covid_risk_factors
    ),
    # Heart Conditions
    num_covid_risk_factors = case_when(
      (CVDCRHD4 == 1) ~ num_covid_risk_factors + 1,
      TRUE ~ num_covid_risk_factors
    ),
    # Diabetes
    num_covid_risk_factors = case_when(
      (DIABETE3 == 1) ~ num_covid_risk_factors + 1,
      TRUE ~ num_covid_risk_factors
    ),
    num_covid_risk_factors = case_when(
      one_of_key_vars_na ~ NA_real_,
      TRUE ~ num_covid_risk_factors
    ),
    # High risk if >= 2 risk factors
    covidhr = (num_covid_risk_factors >= 2),
    covidmr = (num_covid_risk_factors >= 1),
    # If NA, one question was not answered and the obs
    # is not considered "known"
    covidhrknown = case_when(
      is.na(covidhr) ~ FALSE,
      TRUE ~ covidhr
    ),
    covidmrknown = case_when(
      is.na(covidmr) ~ FALSE,
      TRUE ~ covidmr
    )
  ) %>%
  # Drop those with missing age, sex, and race data as
  # there is no such missing data in the PUMS
  filter(
    X.AGEG5YR != 14,
    SEX1 != 7 & SEX1 != 9,
    X.HISPANC != 9,
    ! is.na(X.MRACE1),
    X.MRACE1 != 77 & X.MRACE1 != 99
  )

# *** Compute Risk Shares *************************************************************
# Compute proportion by age, sex, race
options(survey.lonely.psu = "adjust")
brfssdsgn <- survey::svydesign(id = ~1, strata = ~X.STSTR, weights = ~X.LLCPWT, data = clean_BRFSS)
# High risk
hrknown_agesexrace <- survey::svyby(~covidhrknown, by = ~interaction(X.AGEG5YR, SEX1, X.MRACE1, X.HISPANC), FUN = svymean, design = brfssdsgn)
rownames(hrknown_agesexrace) <- hrknown_agesexrace$`interaction(X.AGEG5YR, SEX1, X.MRACE1, X.HISPANC)`
# Medium risk
mrknown_agesexrace <- survey::svyby(~covidmrknown, by = ~interaction(X.AGEG5YR, SEX1, X.MRACE1, X.HISPANC), FUN = svymean, design = brfssdsgn)
rownames(mrknown_agesexrace) <- mrknown_agesexrace$`interaction(X.AGEG5YR, SEX1, X.MRACE1, X.HISPANC)`

#' Gives probability of having a COVID-19 risky
#' medical condition based on age, sex, and race/ethnicity using BRFSS data
#' 
#' @param: age<int>                  | AGEP in PUMS
#' @param: sex<int>                  | SEX in PUMS
#' @param: hisp<int>                 | HISP in PUMS (1 : Not Hispanic, >1 : Hispanic)
#' @param: race_not_incl_hisp<int>   | RAC1P in PUMS
#' @param: risk_type<chr>            | {"high", "medium"} (BRFSS risk profile)
#' 
#' @return: probability<dbl>         | (0,1)
find_risk_probability <- function(age, sex, hisp, race_not_incl_hisp, risk_type) {
  # BRFSS age category
  if (age < 25) {
    brfss_age = 1
  } else if (age >= 80) {
    brfss_age = 13
  } else {
    brfss_age = floor(age/5 - 3)
  }
  # BRFSS Non-Hispanic race category
  if (race_not_incl_hisp == 2) {
    brfss_race = 2   # Black only
  } else if (race_not_incl_hisp == 1) {
    brfss_race = 1   # White only
  } else if (race_not_incl_hisp %in% c(3, 4, 5)) {
    brfss_race = 3   # American Indian or Alaskan Native only
  } else if (race_not_incl_hisp == 6) {
    brfss_race = 4   # Asian only
  } else if (race_not_incl_hisp == 7) { 
    brfss_race = 5   # Native Hawaiian or other Pacific Islander only
  } else if (race_not_incl_hisp == 8) {
    brfss_race = 6   # Some other race only, non-Hispanic
  } else {
    brfss_race = 7   # Two or more races
  }
  # BRFSS Hispanic category
  if (hisp == 1) {
    brfss_hisp = 2   # Not Spanish/Hispanic/Latino
  } else {
    brfss_hisp = 1
  }

  brfss_row = brfss_age + 13 * (sex - 1) + 13 * 2 * (brfss_race - 1) + 13 * 2 * 7 * (brfss_hisp - 1)
  if (risk_type == "high") {
    probability = hrknown_agesexrace$covidhrknownTRUE[brfss_row]
  } else if (risk_type == "medium") {
    probability = mrknown_agesexrace$covidmrknownTRUE[brfss_row]
  } else {
    print("ERROR: Invalid Risk Type")
  }
  return(probability)
}

# ************************************************************************************* #
# ---------------------------------- SOC X-Walking ------------------------------------ #
# ************************************************************************************* #

# *** List of 5-yr ACS PUMS SOC Codes *************************************************
pums_soc_codes <- list()
for (index in file_keys) {
  pums_soc_codes[[index]] <- raw_pums_files[[index]] %>%
    select(
      SOCP
    )
}
pums_soc_codes <- dplyr::bind_rows(pums_soc_codes) %>%
  distinct() %>%
  filter(
    SOCP != ""
  )

# *** Clean X-walks *******************************************************************
# O*NET -> 2018 SOC
onet_to_2018_soc_xwalk <- dplyr::as_tibble(read.csv(
  file = stringr::str_c(
    raw_data_dir,
    "/soc_xwalks/2010_to_2018_SOC_Crosswalk.csv"
    )
  )) %>%
  rename(
    onet_code       = `O.NET.SOC.2010.Code`,
    onet_title      = `O.NET.SOC.2010.Title`,
    soc_2018_code   = `X2018.SOC.Code`,
    soc_2018_title  = `X2018.SOC.Title`
  )

# Disease Exposure by O*NET
disease_exposure_by_onet <- dplyr::as_tibble(read.csv(
  file = stringr::str_c(
    raw_data_dir,
    "/soc_xwalks/Exposed_to_Disease_or_Infections.csv"
    )
  )) %>%
  rename(
    disease_exposure = Context,
    onet_code        = Code,
    onet_title       = Occupation,
  ) %>%
  select(
    disease_exposure, onet_code
  ) %>%
  mutate(
    disease_exposure = case_when(
      disease_exposure == "Not available" ~ NA_character_,
      TRUE ~ disease_exposure
    )
  )

# Critical Infrastructure by 2018 SOC
crit_infr_status_by_2018_soc <- dplyr::as_tibble(openxlsx::read.xlsx(
  stringr::str_c(
    raw_data_dir,
    "/soc_xwalks/SOC-Codes-CISA-Critical-Infrastructure-Workers-with-OES-Data-Rev-1.xlsx"
    ),
    sheet = "2018 SOC Codes"
  )) %>%
  select(
    `SOC.Occupation`, Critical
  ) %>%
  rename(
    soc_2018_code           = `SOC.Occupation`,
    critical_worker_status  = Critical
  ) %>%
  mutate(
    critical_worker_status = case_when(
      (critical_worker_status == "X") ~ 1,
      is.na(critical_worker_status)   ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  assertr::verify(
    ! is.na(critical_worker_status)
  )

# 2018 SOC -> Collapsed SOC
soc_2018_to_collapsed_soc <- dplyr::as_tibble(openxlsx::read.xlsx(
  stringr::str_c(
    raw_data_dir,
    "/soc_xwalks/2018_acs_soc_to_collapsed_soc.xlsx"
    ),
    sheet = "2018_acs_soc_to_collapsed_soc"
  )) %>%
  rename(
    SOCP                 = `2018_soc`,
    SOCP_manual_xwalked  = pums_soc,
  ) %>%
  mutate(
    SOCP = stringr::str_replace(SOCP, "-", ""),
    SOCP_manual_xwalked = stringr::str_replace(SOCP_manual_xwalked, "-", "")
  ) %>%
  mutate(
    pre_str_length = stringr::str_length(trimws(SOCP)),
    SOCP = as.character(as.numeric(trimws(SOCP))),
    post_str_length = stringr::str_length(SOCP)
  ) %>%
  assertr::verify(
    pre_str_length == post_str_length
  ) %>%
  select(
    -pre_str_length, -post_str_length
  )

# *** Complete Merge ******************************************************************
master_soc_xwalk <- onet_to_2018_soc_xwalk %>%
  dplyr::left_join(
    disease_exposure_by_onet,
    by = c("onet_code")
  ) 
# Replace missing with 2-dig average
two_dig_soc_disease_exposure <- master_soc_xwalk %>%
  select(
    soc_2018_code,
    disease_exposure
  ) %>%
  filter(
    ! is.na(disease_exposure)
  ) %>%
  mutate(
    disease_exposure = as.numeric(disease_exposure)
  ) %>%
  group_by(soc_2018_code) %>%
  summarize(
    disease_exposure = mean(disease_exposure)
  ) %>%
  ungroup() %>%
  mutate(
    two_dig_soc = substr(soc_2018_code, 1, 2)
  ) %>%
  group_by(two_dig_soc) %>%
  summarize(
    `2_dig_disease_exposure` = mean(disease_exposure)
  ) %>%
  ungroup()

master_soc_xwalk <- master_soc_xwalk %>%
  mutate(
    two_dig_soc = substr(soc_2018_code, 1, 2)
  ) %>%
  dplyr::left_join(
    two_dig_soc_disease_exposure,
    by = c("two_dig_soc")
  ) %>%
  select(
    -two_dig_soc
  ) %>%
  mutate(
    disease_exposure = case_when(
      is.na(disease_exposure) ~ `2_dig_disease_exposure`,
      TRUE ~ as.numeric(disease_exposure)
    )
  ) %>%
  select(
    -`2_dig_disease_exposure`
  ) %>%
  dplyr::left_join(
    crit_infr_status_by_2018_soc,
    by = c("soc_2018_code")
  ) %>%
  mutate(
    SOCP = stringr::str_replace(soc_2018_code, "-", ""),
    SOCP_pums = case_when(
      (SOCP %in% pums_soc_codes$SOCP) ~ SOCP,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::left_join(
    soc_2018_to_collapsed_soc,
    by = c("SOCP")
  ) %>%
  mutate(
    SOCP_pums = case_when(
      is.na(SOCP_pums) ~ SOCP_manual_xwalked,
      TRUE ~ SOCP_pums
    )
  )
# Match broad categories of occupation codes 
unmatched_SOCP_pums <- pums_soc_codes %>%
  filter(
    ! SOCP %in% master_soc_xwalk$SOCP_pums
  )

#' @param: x<chr>
#' @param: n<int>
#' 
#' @return: rightmost nth chars of x
substrRight <- function(x, n) {
  return(substr(x, nchar(x)-n+1, nchar(x)))
}

manual_xwalk <- c(NA, NA)
for (i in 1:nrow(unmatched_SOCP_pums)) {
  if (substrRight(unmatched_SOCP_pums[i, 1], 1) == "0") {
    for (j in 1:9) {
      manual_xwalk <- rbind(
        manual_xwalk,
        c(
          stringr::str_c(
            substr(unmatched_SOCP_pums[i, 1], 1, 5),
            j
          ),
          as.character(unmatched_SOCP_pums[i, 1])
        )
      )
    }
  }
}

manual_xwalk <- dplyr::as_tibble(manual_xwalk) %>%
  rename(
    SOCP = V1,
    SOCP_manually_xwalked = V2
  ) %>%
  filter(
    ! is.na(SOCP)
  )

master_soc_xwalk <- master_soc_xwalk %>%
  dplyr::left_join(
    manual_xwalk,
    by = c("SOCP")
  ) %>%
  mutate(
    SOCP_pums = case_when(
      is.na(SOCP_pums) ~ SOCP_manually_xwalked,
      TRUE ~ SOCP_pums
    ),
    # Complete last few manually
    SOCP_pums = case_when(
      substr(SOCP, 1, 3) == "251" ~ "251000",
      substr(SOCP, 1, 3) == "391" ~ "391000",
      substr(SOCP, 1, 5) %in% c("47501", "47505", "47508", "47509") ~ "4750XX",
      substr(SOCP, 1, 3) == "299" ~ "299000",
      substr(SOCP, 1, 3) == "531" ~ "531000",
      substr(SOCP, 1, 5) %in% c("51406", "51407") ~ "5140XX",
      TRUE ~ SOCP_pums
    )
  ) %>%
  rename(
    pums_soc_2018_code = SOCP_pums
  ) %>%
  filter(
    ! is.na(pums_soc_2018_code)
  ) %>%
  select(
    pums_soc_2018_code,
    disease_exposure,
    critical_worker_status
  ) %>%
  group_by(pums_soc_2018_code) %>%
  summarize(
    disease_exposure        = mean(disease_exposure),
    critical_worker_status  = mean(critical_worker_status)
  ) %>%
  rename(
    SOCP = pums_soc_2018_code
  )

for (i in file_keys) {
  pums_files[[i]] <- pums_files[[i]] %>%
    dplyr::left_join(
      master_soc_xwalk,
      by = c("SOCP")
    )
}

# ************************************************************************************* #
# ----------------------------- NASEM Phase Labeling ---------------------------------- #
# ************************************************************************************* #
# Make replicable
set.seed(42)
# *** Phase 1A ************************************************************************
# High risk health care workers
hc_worker_naics <- c(
  "622M",     # 622M - MED-General Medical And Surgical Hospitals, And Specialty (Except Psychiatric And Substance Abuse) Hospitals
  "6231",     # 6231 - MED-Nursing Care Facilities (Skilled Nursing Facilities)
  "623M",     # 623M - MED-Residential Care Facilities, Except Skilled Nursing Facilities
  "6216"      # 6216 - MED-Home Health Care Services
)
morticians_socs <- c(
  "394031",   # 394031 - PRS-Morticians, Undertakers, And Funeral Arrangers
  "3940XX"    # 3940XX - PRS-Embalmers, Crematory Operators, And Funeral Attendants
)
pharmacists_socs <- c(
  "292052",   # 292052 - MED-Pharmacy Technicians
  "291051"    # 291051 - MED-Pharmacists
)
pub_health_worker_socs <- c(
  "21109X",   # 21109X - CMS-Other Community And Social Service Specialists (includes Community Health Workers)
  "211022"    # 211022 - CMS-Healthcare Social Workers
)
dentist_socs <- c(
  "291020",   # 291020 - MED-Dentist
  "319091"    # 319091 - HLS-Dental Assistants
)
# From Baker et al., 96.1% and 91.5% of SOC Codes corresponding to HC professions (Two digit SOC in {31, 29} respectively)
# are exposed >1 time/month to infection / disease in the workplace. Make a very crude assumption and randomly sample
# the average of the two. Another crude assumption here is that morticians/pharmacists/public health workers/dentists
# face similar exposure levels.
fraction_exposed_to_infection = (.961 + .915) / 2
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_hc_worker_in_risky_env = (NAICSP %in% hc_worker_naics),
      is_mortician              = (SOCP %in% morticians_socs),
      is_pharma                 = (SOCP %in% pharmacists_socs),
      is_pub_health             = (SOCP %in% pub_health_worker_socs),
      is_dentist                = (SOCP %in% dentist_socs)
    ) %>%
    mutate(
      is_ph1a_high_risk_hc_worker = (
        (
          is_hc_worker_in_risky_env | 
          is_mortician | 
          is_pharma | 
          is_pub_health | 
          is_dentist
        ) & sample(
          c(TRUE, FALSE),
          size = pums_dims[[index]],
          replace = TRUE,
          prob = c(fraction_exposed_to_infection, 1 - fraction_exposed_to_infection)
        )
      )
    )
}
# First responder
first_responders <- c(
  "292042",    # 292042 - MED-Emergency Medical Technicians
  "292043",    # 292043 - MED-Paramedics
  "533011",    # 533011 - TRN-Ambulance Drivers And Attendants, Except Emergency Medical Technicians
  "331012",    # 331012 - PRT-First-Line Supervisors Of Police And Detectives
  "333050",    # 333050 - PRT-Police Officers
  "331021",    # 331021 - RT-First-Line Supervisors Of Fire Fighting And Prevention Workers
  "332011"     # 332011 - PRT-Firefighters
)
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_ph1a_first_responder = (SOCP %in% first_responders),
    ) 
}

for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      phase_1a = is_ph1a_high_risk_hc_worker | is_ph1a_first_responder
    ) 
}

# *** Phase 1B ************************************************************************
# High risk
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      risk_probabilities_high = mapply(find_risk_probability, AGEP, SEX, HISP, RAC1P, "high"),
      # Choose a uniform random number from 0 to 1. Assign someone as high risk if their random number 
      # is less than their risk probability. Assume 18-24 risk probabilities for people <18 years old
      # dut to lack of data in BRFSS
      is_high_risk = (runif(pums_dims[[index]]) < risk_probabilities_high)
    ) 
}
# Older adult in crowded housing
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      # MULTG == 2 -> multi-gernational housing | TYPE == 2 -> institutional group quarters
      is_ph1b_older_adult_crowd_hous = (AGEP >= 65 & (MULTG == 2 | TYPE == 2)),
      is_ph1b_older_adult_crowd_hous = case_when(
        is.na(is_ph1b_older_adult_crowd_hous) ~ FALSE,
        TRUE ~ is_ph1b_older_adult_crowd_hous
      ) 
    ) 
}
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      phase_1b = (is_high_risk | is_ph1b_older_adult_crowd_hous) & ! phase_1a
    ) 
}
# *** Phase 2 *************************************************************************
# Essential workers and exposed at work
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_critical_worker_and_exposed = (critical_worker_status > 0) & (disease_exposure >= 50),
      is_critical_worker_and_exposed = case_when(
        is.na(is_critical_worker_and_exposed) ~ FALSE,
        TRUE ~ is_critical_worker_and_exposed
      )
    ) 
}
# Teachers and school staff
teachers_and_school_staff <- c(
  "252010",   # 252010 - EDU-Preschool And Kindergarten Teachers
  "252020",   # 252020 - EDU-Elementary And Middle School Teachers
  "252030",   # 252030 - EDU-Secondary School Teachers
  "252050",   # 252050 - EDU-Special Education Teachers
  "2530XX",   # 2530XX - EDU-Other Teachers And Instructors
  "259040",   # 259040 - EDU-Teaching Assistants
  "2590XX",   # 2590XX - EDU-Other Educational Instruction and Library Workers
  "193034",   # 193034 - SCI-School Psychologists
  "339094",   # 339094 - PRT-School Bus Monitors
  "533051",   # 533051 - TRN-Bus Drivers, School
  "119030",   # 119030 - MGR-Education And Childcare Administrators
  "211021",   # 211021 - CMS-Child, Family, And School Social Workers
  "211012",   # 211012 - CMS-Educational, Guidance, And Career Counselors And Advisors
  "399011"    # 399011 - PRS-Childcare Workers
)
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_teacher_and_school_staff = (SOCP %in% teachers_and_school_staff)
    ) 
}
# Medium risk
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      # Choose a uniform random number from 0 to 1. Assign someone as medium risk if their random number 
      # is less than their risk probability. Assume 18-24 risk probabilities for people <18 years old
      # dut to lack of data in BRFSS
      risk_probabilities_medium = mapply(find_risk_probability, AGEP, SEX, HISP, RAC1P, "medium"),
      is_medium_risk = (runif(pums_dims[[index]]) < risk_probabilities_medium)
    ) 
}
# Older adults
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_older_adult = (AGEP >= 65)
    ) 
}
# Homeless or group home
military_naics <- c(
  "928110P1",    # 928110P1 .MIL-U.S. Army
  "928110P2",    # 928110P2 .MIL-U.S. Air Force
  "928110P3",    # 928110P3 .MIL-U.S. Navy
  "928110P4",    # 928110P4 .MIL-U.S. Marines
  "928110P5",    # 928110P5 .MIL-U.S. Coast Guard
  "928110P6",    # 928110P6 .MIL-U.S. Armed Forces, Branch Not Specified
  "928110P7"     # 928110P7 .MIL-Military Reserves Or National Guard
)
# (Imperfect approximation) Social and Community Service Managers
shelter_group_home_staff <- c("119151")

for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_active_duty = (MIL == 1),
      is_active_duty = case_when(
        is.na(is_active_duty) ~ FALSE,
        TRUE ~ is_active_duty
      ),
      is_undergrad = (SCHG == 15),
      is_undergrad = case_when(
        is.na(is_undergrad) ~ FALSE,
        TRUE ~ is_undergrad
      ),
      # Housing type is not granular enough to indicate
      # whether or not someone is in a homeless shelter. We can approximate this
      # by including non-institutional group quarters residents. According to
      # the census (https://www.census.gov/topics/income-poverty/poverty/guidance/group-quarters.html)
      # this includes: (1) College dormitories, (2) Military baracks, (3) Group homes,
      # (4) Missions, and (5) Shelters. To mitigate this , we don't consider active 
      # military members, people with a military NAICs code and college undergraduate students. 
      is_homeless_or_gr_home = (
        ! NAICSP %in% military_naics & 
        ! is_active_duty & 
        (TYPE == 3 | SOCP %in% shelter_group_home_staff) & 
        ! is_undergrad
      )
    ) 
}
# Prisoners or staff
prison_staff <- c(
  "331011",    # 331011 - PRT-First-Line Supervisors Of Correctional Officers
  "333012"     # 333012 - PRT-First-Line Supervisors Of Police And Detectives
)
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      # Define prisoners as institional group quartered population. Note that
      # this pool was already dipped into for Phase 1B. Note further that this
      # will have the unfortunate side effect of including nursing home residents
      # and mental hospital residents along with prisoners.
      is_prisoner_or_staff = (SOCP %in% prison_staff) | (TYPE == 2)
    )
}
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      phase_2 = (
        (
          is_critical_worker_and_exposed | 
          is_teacher_and_school_staff | 
          is_medium_risk | 
          is_older_adult | 
          is_homeless_or_gr_home | 
          is_prisoner_or_staff
        ) & (! phase_1a) & (! phase_1b)
      )
    )
}

# *** Phase 3 *************************************************************************
# Young adults
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_young_adult = (AGEP >= 18 & AGEP <= 30)
    )
}
# Children
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_child = (AGEP < 18)
    )
}
# Essential worker and moderately exposed (liberal definition taken)
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      is_critical_worker_and_med_exposed = (critical_worker_status > 0),
      is_critical_worker_and_med_exposed = case_when(
        is.na(is_critical_worker_and_med_exposed) ~ FALSE,
        TRUE ~ is_critical_worker_and_med_exposed
      )
    )
}

for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      phase_3 = (
        (
          is_young_adult | 
          is_child | 
          is_critical_worker_and_med_exposed
        ) & (! phase_1a) & (! phase_1b) & (! phase_2)
      )
    )
}

# *** Phase 4 *************************************************************************
# Everyone else
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    mutate(
      phase_4 = (
        (! phase_1a) & 
        (! phase_1b) & 
        (! phase_2)  & 
        (! phase_3)
      )
    )
}

# ************************************************************************************* #
# ------------------------------------- Save ------------------------------------------ #
# ************************************************************************************* #
for (index in file_keys) {
  pums_files[[index]] <- pums_files[[index]] %>%
    select(
      PWGTP,
      AGEP,
      ST,
      id,
      phase_1a,
      phase_1b,
      phase_2,
      phase_3,
      phase_4,
    )
}
pums_nasem_phase_tagged_df <- dplyr::bind_rows(pums_files)

data.table::fwrite(
  pums_nasem_phase_tagged_df, 
  stringr::str_c(
    build_data_dir,
    "/pums_data_with_nasem_phases.csv"
  )
)