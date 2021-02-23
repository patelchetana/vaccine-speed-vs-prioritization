# *** Header **************************************************************************
#
# Create .csv version of Table S2
#

# Read in initial conditions
initial_conditions <- read.csv(
  stringr::str_c(
    build_data_dir,
    "/seir_model_nat_lvl_infection_starting_points_by_group_on_dec14.csv"
  )
)

# Format table
table_s2 <- initial_conditions %>%
  mutate(
    group = stringr::str_replace(
      group,
      "agebin_",
      ""
    ),
    `R(0)` = round(r0 - d0, 4),
    d0 = as.character(round(d0, 4)),
    d0 = case_when(
      d0 == "0" ~ "<0.0001",
      TRUE      ~ d0
    ),
    s0 = round(s0, 4),
    e0 = round(e0, 4),
    i0 = round(i0, 4),
    r0 = round(r0, 4),
  ) %>%
  relocate(
    group, s0, e0, i0, `R(0)`, d0, r0
  ) %>%
  rename(
    `Age group`   = group,
    `S(0)`        = s0,
    `E(0)`        = e0,
    `I(0)`        = i0,
    `D(0)`        = d0,
    `R(0) + D(0)` = r0,
  )
  
# Write the table
write.csv(
  table_s2,
  stringr::str_c(
    exhibit_data_dir,
    "/table_s2.csv"
  ),
  row.names = FALSE
)