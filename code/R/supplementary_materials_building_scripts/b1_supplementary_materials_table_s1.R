# *** Header **************************************************************************
#
# Create .csv version of Table S1
#

# Read in contact matrix
contact_mat <- read.csv(
  stringr::str_c(
    build_data_dir,
    "/age_group_contact_matrix.csv"
  )
)

# Format table
table_s1 <- contact_mat %>%
  mutate(
    X = stringr::str_replace(
      X,
      "agebin_",
      ""
    )
  ) %>%
  rename(
    `0-9`   = agebin_0.9,
    `10-19` = agebin_10.19,
    `20-29` = agebin_20.29,
    `30-39` = agebin_30.39,
    `40-49` = agebin_40.49,
    `50-59` = agebin_50.59,
    `60-69` = agebin_60.69,
    `70-79` = agebin_70.79,
    `80+`   = agebin_80.
  ) %>%
  mutate(
    `0-9`    = round(`0-9`, 2),
    `10-19`  = round(`10-19`, 2),
    `20-29`  = round(`20-29`, 2),
    `30-39`  = round(`30-39`, 2),
    `40-49`  = round(`40-49`, 2),
    `50-59`  = round(`50-59`, 2),
    `60-69`  = round(`60-69`, 2),
    `70-79`  = round(`70-79`, 2),
    `80+`    = round(`80+`, 2),
  ) %>%
  rename(
    `Contactor` = X
  )

# Write the table
write.csv(
  table_s1,
  stringr::str_c(
    exhibit_data_dir,
    "/table_s1.csv"
  ),
  row.names = FALSE
)