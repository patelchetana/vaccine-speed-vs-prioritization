# *** Header **************************************************************************
#
# Create .csv version of Table S3
#

# Read in national baseline
national_baseline <- read.csv(
  stringr::str_c(
    build_data_dir,
    "/national_baseline.csv"
  )
)

# Format the table
table_s3 <- national_baseline %>%
  select(
    group,
    n_group,
    c19_ifr_group,
    sus_to_inf,
    vax_uptake_census,
    average_2vax_efficacy,
    yll,
  ) %>%
  mutate(
    group = stringr::str_replace(
      group,
      "agebin_",
      ""
    ),
    d_E                   = "4 days (age invariant)",
    d_I                   = "9 days (age invariant)",
    sus_to_inf            = round(sus_to_inf, 2),
    c19_ifr_group         = round(c19_ifr_group, 3),
    n_group               = prettyNum(n_group, big.mark = ","),
    vax_uptake_census     = round(vax_uptake_census, 3),
    c_ij                  = "See Table S1",
    average_2vax_efficacy = round(average_2vax_efficacy, 3),
    yll                   = round(yll, 1)
  ) %>%
  rename(
    `Age group` = group,
    beta_i      = sus_to_inf,
    IFR_i       = c19_ifr_group,
    N_i         = n_group,
    vu_i        = vax_uptake_census,
    ve_i        = average_2vax_efficacy,
    YLL_i       = yll
  ) %>%
  relocate(
    `Age group`,
    d_E,
    d_I,
    beta_i,
    IFR_i,
    N_i,
    vu_i,
    c_ij,
    ve_i,
    YLL_i
  )

# Write the table
write.csv(
  table_s3,
  stringr::str_c(
    exhibit_data_dir,
    "/table_s3.csv"
  ),
  row.names = FALSE
)