# *** Header **************************************************************************
#
# This script allocates vaccines from federal -> state -> individual
# under NASEM guidelines given no reserve while incorporating parameter
# assumptions around vaccine uptake stratified by age.
#

# *** Define Functions to Help Out with Allocation ************************************
#' Initialize stateful state-level matrix to mantain allocation state over time
#' 
#' @param: pums_nasem_phase_labeled<df>   | At least 8 fields: PWGTP<int>, AGEP<int>, ST<int>, phase_XX<logical>
#' @param: uptake_by_age_bin<df>          | Two fields: group<chr>, uptake<dbl>
#' 
#' @return: state_lvl_stateful_mat<df>    | unique by ST<int>, phase<chr>, group<chr> with fields as n<int>, uptake<dbl>, noff<dbl>
#'                                                            naccept<dbl>, ndecline<dbl>, sh<dbl>, nunoff<dbl>, st_n<dbl>
create_state_level_stateful_mat <- function(pums_nasem_phase_labeled, uptake_by_age_bin) {
  state_lvl_stateful_mat = pums_nasem_phase_labeled %>%
    select(
      PWGTP,
      ST,
      AGEP,
      phase_1a, phase_1b, phase_2, phase_3, phase_4
    ) %>%
    mutate(
      phase = case_when(
        phase_1a ~ "p1a",
        phase_1b ~ "p1b",
        phase_2  ~ "p2",
        phase_3  ~ "p3",
        phase_4  ~ "p4",
        TRUE ~ NA_character_
      )
    ) %>%
    assertr::verify(
      ! is.na(phase)
    ) %>%
    select(-starts_with("phase_")) %>%
    mutate(
      group = case_when(
        AGEP < 10              ~ "0_9",
        AGEP < 16 & AGEP >= 10 ~ "10_15",
        AGEP < 20 & AGEP >= 16 ~ "16_19",
        AGEP < 30 & AGEP >= 20 ~ "20_29",
        AGEP < 40 & AGEP >= 30 ~ "30_39",
        AGEP < 50 & AGEP >= 40 ~ "40_49",
        AGEP < 60 & AGEP >= 50 ~ "50_59",
        AGEP < 70 & AGEP >= 60 ~ "60_69",
        AGEP < 80 & AGEP >= 70 ~ "70_79",
        AGEP >= 80             ~ "80pl",
        TRUE                   ~ NA_character_
      )
    ) %>%
    assertr::verify(
      ! is.na(group)
    ) %>%
    select(-AGEP) %>%
    mutate(
      group = as.factor(group)
    ) %>%
    group_by(ST, phase, group, .drop = FALSE) %>%
    summarize(
      n = sum(PWGTP)
    ) %>%
    ungroup() %>%
    mutate(
      group = as.character(group)
    ) %>%
    group_by(ST, phase) %>%
    mutate(
      sh = n / sum(n)
    ) %>%
    ungroup() %>%
    dplyr::left_join(
      uptake_by_age_bin,
      by = c("group")
    ) %>%
    assertr::verify(
      ! is.na(uptake)
    ) %>%
    mutate(
      noff          = 0,
      naccept       = 0,
      ndecline      = 0,
      nunoff        = n,
      neventaccept  = n * uptake,
    ) %>%
    group_by(ST, phase) %>%
    mutate(
      sh_eventaccept = neventaccept / sum(neventaccept),
      ph_uptake = weighted.mean(uptake, n)
    ) %>%
    ungroup() %>%
    group_by(ST) %>%
    mutate(
      st_n = sum(n)
    )
  return(state_lvl_stateful_mat)
}

#' Compute single round of vax allocation by updating param stateful matrix
#' 
#' @param: state_lvl_stateful_mat<df>           | At least 8 fields: PWGTP<int>, AGEP<int>, ST<int>, phase_XX<logical>, uptake<dbl>
#' @param: num_incremental_units<int>           | 100,000
#' 
#' @return: updated_state_lvl_stateful_mat<df>  | State of @param state_lvl_stateful_matrix after allocation @param num_incremental_units
compute_round_fed_to_st_to_phase_unreserved_nasem <- function(state_lvl_stateful_mat, num_incremental_units) {
  updated_state_lvl_stateful_mat = state_lvl_stateful_mat %>%
    group_by(ST) %>%
    mutate(
      # S' : states with some who have not yet been offered (Phase 4 still has people yet to be offered)
      is_s_prime = round(sum(nunoff[phase == "p4"]), 5) > 0,
    ) %>%
    ungroup() %>%
    mutate(
      # lambda_s : state population share among S'
      sum_N_s_prime = sum(n[is_s_prime]),
      lambda_s = case_when(
        is_s_prime ~ st_n / sum_N_s_prime,
        TRUE ~ 0
      ),
      # V_s : supply each state receives from @param: num_incremental_units
      V_s = lambda_s * num_incremental_units
    ) %>%
    group_by(ST) %>%
    mutate(
      # --- Phase 1A ---
      V_s_working = case_when(
        (phase == "p1a") ~ V_s * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p1a_accepts = pmin(V_s_working, neventaccept - naccept),
      V_s = pmax(V_s - sum(incr_p1a_accepts), 0),
      # --- Phase 1B ---
      V_s_working = case_when(
        (phase == "p1b") ~ V_s * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p1b_accepts = pmin(V_s_working, neventaccept - naccept),
      V_s = pmax(V_s - sum(incr_p1b_accepts), 0),
      # --- Phase 2 ---
      V_s_working = case_when(
        (phase == "p2") ~ V_s * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p2_accepts = pmin(V_s_working, neventaccept - naccept),
      V_s = pmax(V_s - sum(incr_p2_accepts), 0),
      # --- Phase 3 ---
      V_s_working = case_when(
        (phase == "p3") ~ V_s * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p3_accepts = pmin(V_s_working, neventaccept - naccept),
      V_s = pmax(V_s - sum(incr_p3_accepts), 0),
      # --- Phase 4 ---
      V_s_working = case_when(
        (phase == "p4") ~ V_s * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p4_accepts = pmin(V_s_working, neventaccept - naccept),
      V_s = pmax(V_s - sum(incr_p4_accepts), 0),
      # --- Tally Accepts ----
      incr_accepts = (
        incr_p1a_accepts +
        incr_p1b_accepts +
        incr_p2_accepts  +
        incr_p3_accepts  +
        incr_p4_accepts
      )
    ) %>%
    group_by(ST, phase) %>%
    mutate(
      # --- Offers ----
      ph_incr_offers = case_when(
        phase == "p1a" ~ sum(incr_p1a_accepts) / ph_uptake,
        phase == "p1b" ~ sum(incr_p1b_accepts) / ph_uptake,
        phase == "p2"  ~ sum(incr_p2_accepts)  / ph_uptake,
        phase == "p3"  ~ sum(incr_p3_accepts)  / ph_uptake,
        phase == "p4"  ~ sum(incr_p4_accepts)  / ph_uptake
      ),
      incr_offers = ph_incr_offers * sh,
      incr_declines = incr_offers - incr_accepts
    ) %>%
    ungroup() %>%
    mutate(
      noff     = noff + incr_offers,
      nunoff   = nunoff - incr_offers,
      naccept  = naccept + incr_accepts,
      ndecline = ndecline + incr_declines
    ) %>%
    assertr::verify(
      round(noff - naccept - ndecline, 5) == 0
    ) %>%
    assertr::verify(
      round(n - noff - nunoff, 5) == 0
    ) %>%
    group_by(ST) %>%
    mutate(
      state_is_complete = (round(sum(n) - sum(noff), 5) == 0),
      state_made_offers = sum(incr_offers) > 0,
      state_is_complete_and_made_offers = (state_is_complete & state_made_offers)
    ) %>%
    ungroup() %>%
    assertr::verify(
      round(sum(incr_accepts) - num_incremental_units, 5) == 0 | sum(state_is_complete_and_made_offers) > 0 | sum(noff) > 320000000
    ) %>%
    group_by(ST) %>%
    assertr::verify(
      sum(incr_accepts) > 0 | ! is_s_prime
    ) %>%
    ungroup() %>%
    # Save "wasted" units that accrue because of boundary condition.
    # For lack of better alternative, split the leftover by share of remaining units
    # to be allocated
    mutate(
      leftover = num_incremental_units - sum(incr_accepts)
    ) %>%
    # Clear out variables from heart of allocation
    select(
      -starts_with("incr"), -is_s_prime, -sum_N_s_prime, -lambda_s,
      -V_s, -V_s_working, -ph_incr_offers, -starts_with("state_is"), -starts_with("state_made")
    ) %>%
    group_by(ST) %>%
    mutate(
      amount_remaining        = neventaccept - naccept,
      st_amount_remaining     = sum(amount_remaining),
      has_remaining_to_alloc  = round(st_amount_remaining, 5) > 0,
    ) %>%
    ungroup() %>%
    mutate(
      # lambda_r : state share of remaining people to be allocated
      sum_N_state_remaining = sum(amount_remaining[has_remaining_to_alloc]),
      lambda_r = case_when(
        has_remaining_to_alloc ~ st_amount_remaining / sum_N_state_remaining,
        TRUE ~ 0
      ),
      # Split the leftover to states by share of remaining to be allocated
      V_r = lambda_r * leftover
    ) %>%
    group_by(ST) %>%
    mutate(
      # --- Phase 1A ---
      V_r_working = case_when(
        (phase == "p1a") ~ V_r * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p1a_accepts = pmin(V_r_working, neventaccept - naccept),
      V_r = pmax(V_r - sum(incr_p1a_accepts), 0),
      # --- Phase 1B ---
      V_r_working = case_when(
        (phase == "p1b") ~ V_r * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p1b_accepts = pmin(V_r_working, neventaccept - naccept),
      V_r = pmax(V_r - sum(incr_p1b_accepts), 0),
      # --- Phase 2 ---
      V_r_working = case_when(
        (phase == "p2") ~ V_r * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p2_accepts = pmin(V_r_working, neventaccept - naccept),
      V_r = pmax(V_r - sum(incr_p2_accepts), 0),
      # --- Phase 3 ---
      V_r_working = case_when(
        (phase == "p3") ~ V_r * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p3_accepts = pmin(V_r_working, neventaccept - naccept),
      V_r = pmax(V_r - sum(incr_p3_accepts), 0),
      # --- Phase 4 ---
      V_r_working = case_when(
        (phase == "p4") ~ V_r * sh_eventaccept,
        TRUE ~ 0
      ),
      incr_p4_accepts = pmin(V_r_working, neventaccept - naccept),
      V_r = pmax(V_r - sum(incr_p4_accepts), 0),
      # --- Tally Accepts ----
      incr_accepts = (
        incr_p1a_accepts +
        incr_p1b_accepts +
        incr_p2_accepts  +
        incr_p3_accepts  +
        incr_p4_accepts
      )
    ) %>%
    group_by(ST, phase) %>%
    mutate(
      # --- Offers ----
      ph_incr_offers = case_when(
        phase == "p1a" ~ sum(incr_p1a_accepts) / ph_uptake,
        phase == "p1b" ~ sum(incr_p1b_accepts) / ph_uptake,
        phase == "p2"  ~ sum(incr_p2_accepts)  / ph_uptake,
        phase == "p3"  ~ sum(incr_p3_accepts)  / ph_uptake,
        phase == "p4"  ~ sum(incr_p4_accepts)  / ph_uptake
      ),
      incr_offers = ph_incr_offers * sh,
      incr_declines = incr_offers - incr_accepts
    ) %>%
    ungroup() %>%
    mutate(
      noff     = noff + incr_offers,
      nunoff   = nunoff - incr_offers,
      naccept  = naccept + incr_accepts,
      ndecline = ndecline + incr_declines
    ) %>%
    assertr::verify(
      round(noff - naccept - ndecline, 5) == 0
    ) %>%
    assertr::verify(
      round(n - noff - nunoff, 5) == 0
    ) %>%
    # Clear out the variables
    select(
      -starts_with("incr"), -leftover, -sum_N_state_remaining, -amount_remaining, -st_amount_remaining, -has_remaining_to_alloc, -lambda_r, -V_r, -V_r_working, -ph_incr_offers
    )

  return(updated_state_lvl_stateful_mat)
}

#' Complete allocation where every person is offered once and time is scaled in terms of accepts
#' 
#' @param: pums_nasem_phase_labeled<df>    | At least 8 fields: PWGTP<int>, AGEP<int>, ST<int>, phase_XX<logical>
#' @param: uptake_by_age_bin<df>           | Two fields: group<chr>, uptake<dbl>
#' @param: size_of_increment<int>          | Usually 100,000
#' 
#' @return: sim_results<df>                | Unique by ST, n_units_accepted_nationally showing NASEM allocation over time
compute_unreserved_nasem_with_uptake <- function(pums_nasem_phase_labeled, uptake_by_age_bin, size_of_increment) {
  sim_results = list()
  st_lvl_stateful_mat = create_state_level_stateful_mat(pums_nasem_phase_labeled, uptake_by_age_bin) %>%
    mutate(
      n_units_accepted_nationally = 0
    )

  lower_bound = 1
  upper_bound = (ceiling(sum(st_lvl_stateful_mat$neventaccept) / size_of_increment) * size_of_increment) / size_of_increment
  sim_results[["0"]] = st_lvl_stateful_mat
  for (i in lower_bound:upper_bound) {
    st_lvl_stateful_mat = compute_round_fed_to_st_to_phase_unreserved_nasem(st_lvl_stateful_mat, size_of_increment) %>%
      mutate(
        n_units_accepted_nationally = i * size_of_increment
      )
    sim_results[[as.character(format(i * size_of_increment, scientific = F))]] = st_lvl_stateful_mat
  }
  sim_results = dplyr::bind_rows(sim_results)
  return(sim_results)
}

#'
#' @param: long_sim_result<df>     | unique by ST<int>, n_units_accepted_nationally<int>, phase<chr>, group<chr>
#' 
#' @return: seir_input_form<df>    | group<chr>, entity<chr>, cumul_n_accepted_nationally<dbl>, 
#'                                   cumul_n_offered_actual<dbl>, cumul_n_group_offered<dbl>, 
#'                                   nat_population<dbl>, n_group<dbl>, sh_group<dbl>, cumul_sh_vaccinated<dbl>
#'                                   n_offered_nationally<dbl>, flow_vaccinated<dbl>
convert_to_seir_input_form <- function(long_sim_result) {
  seir_input_form = long_sim_result %>%
    select(
      -st_n, -nunoff
    ) %>%
    filter(
      (n_units_accepted_nationally != 0 & n_units_accepted_nationally %% 5000000 == 0) | n_units_accepted_nationally == max(n_units_accepted_nationally)
    ) %>%
    group_by(n_units_accepted_nationally, group) %>%
    # Sum over states and phases
    summarize(
      n        = sum(n),
      noff     = sum(noff),
      naccept  = sum(naccept),
      ndecline = sum(ndecline)
    ) %>%
    ungroup() %>%
    mutate(
      group = stringr::str_replace(group, "pl", "+"),
      group = stringr::str_replace(group, "_", "-"), 
      group = stringr::str_c("agebin_", group),
    ) %>%
    # Collapse to SEIR group structure from group structure
    # created to ensure uptake among individuals <16 years old is 0
    mutate(
      seir_group = case_when(
        group == "agebin_10-15" ~ "agebin_10-19",
        group == "agebin_16-19" ~ "agebin_10-19",
        TRUE ~ group
      )
    ) %>%
    select(
      -group
    ) %>%
    rename(
      group = seir_group
    ) %>%
    group_by(n_units_accepted_nationally, group) %>%
    summarize(
      n        = sum(n),
      noff     = sum(noff),
      naccept  = sum(naccept),
      ndecline = sum(ndecline)
    ) %>%
    ungroup() %>%
    arrange(
      group, n_units_accepted_nationally
    ) %>%
    mutate(
      entity = "USA (50 states + DC)",
    ) %>%
    rename(
      cumul_n_accepted_nationally = n_units_accepted_nationally,
      n_group                     = n,
      cumul_n_group_offered       = noff,
      cumul_n_group_accepted      = naccept,
      cumul_n_group_declined      = ndecline
    ) %>%
    group_by(cumul_n_accepted_nationally) %>%
    mutate(
      cumul_n_offered_actual   = sum(cumul_n_group_offered),
      cumul_n_accepted_actual  = sum(cumul_n_group_accepted),
      cumul_n_declined_actual  = sum(cumul_n_group_declined),
      nat_population           = sum(n_group),
      sh_group                 = n_group / nat_population
    ) %>%
    ungroup() %>%
    mutate(
      cumul_n_accepted_nationally = cumul_n_accepted_actual
    ) %>%
    mutate(
      cumul_share_vaccinated = cumul_n_group_accepted / n_group
    ) %>%
    group_by(group) %>%
    mutate(
      # Incremental units offered nationally
      n_offered_nationally = cumul_n_offered_actual - lag(cumul_n_offered_actual),
      n_offered_nationally = case_when(
        is.na(n_offered_nationally) ~ cumul_n_offered_actual,
        TRUE ~ n_offered_nationally
      ),
      # Incremental units accepted nationally
      n_accepted_nationally = cumul_n_accepted_actual - lag(cumul_n_accepted_actual),
      n_accepted_nationally = case_when(
        is.na(n_accepted_nationally) ~ cumul_n_accepted_actual,
        TRUE ~ n_accepted_nationally
      ),
      # Incremental share of group accepted vaccination offer
      flow_vaccinated = cumul_share_vaccinated - lag(cumul_share_vaccinated),
      flow_vaccinated = case_when(
        is.na(flow_vaccinated) ~ cumul_share_vaccinated,
        TRUE ~ flow_vaccinated
      )
    ) %>%
    ungroup() %>%
    relocate(
      group, entity, cumul_n_accepted_nationally, cumul_n_offered_actual, cumul_n_accepted_actual, 
      cumul_n_declined_actual, cumul_n_group_offered, cumul_n_group_accepted, cumul_n_group_declined,
      nat_population, n_group, sh_group, cumul_share_vaccinated, n_offered_nationally, n_accepted_nationally, flow_vaccinated
    ) %>%
    # --- Asserts on Final Form ----
    assertr::verify(
      round((cumul_n_accepted_actual + cumul_n_declined_actual) - cumul_n_offered_actual, 5) == 0
    ) %>%
    assertr::verify(
      round((cumul_n_group_accepted + cumul_n_group_declined) - cumul_n_group_offered, 5) == 0
    ) 
  return(seir_input_form)
}

# *** Read in 5-Year ACS PUMS Data ****************************************************
pums_nasem_labeled <- data.table::fread(
  file = stringr::str_c(
    build_data_dir,
    "/pums_data_with_nasem_phases.csv"
  ),
  sep = ","
)
# *** Read in Uptake Assumptions ******************************************************
national_baseline_file <- dplyr::as_tibble(data.table::fread(
  file = stringr::str_c(
    build_data_dir,
    "/national_baseline.csv"
  ),
  sep = ","
))

uptake_scenarios <- national_baseline_file %>%
  select(
    # Must be all-age assumption as zero-ing of 0-15 year olds is in
    # these scripts
    group, vax_uptake_malik_all_ages, vax_uptake_census_all_ages
  ) %>%
  mutate(
    group = case_when(
      group == "agebin_10-19" ~ "agebin_16-19",
      TRUE ~ group
    )
  ) %>%
  add_row(
    group = "agebin_10-15",
    vax_uptake_malik_all_ages = 0,
    vax_uptake_census_all_ages = 0
  ) %>%
  mutate(
    vax_uptake_100p = 1,
    vax_uptake_const_70 = .7
  ) %>%
  mutate_at(
    vars(-group),
    list(~case_when(
      group %in% c("agebin_0-9", "agebin_10-15") ~ 0,
      TRUE ~ .
    ))
  ) %>%
  mutate(
    group = stringr::str_replace(
      group,
      "agebin_",
      ""
    ),
    group = stringr::str_replace(
      group,
      "-",
      "_"
    ),
    group = stringr::str_replace(
      group,
      "\\+",
      "pl"
    )
  )

# *** Run the Sims Parameterized by Uptake ********************************************
sim_runs <- list()
for (uptake_assumption in c("malik_all_ages", "census_all_ages", "100p", "const_70")) {
  print(stringr::str_c("Uptake parameter of NASEM unreserved rollout simulation: ", uptake_assumption))
  uptake_input_to_sim <- uptake_scenarios %>%
    select(
      group, ends_with(uptake_assumption)
    ) %>%
    rename(
      uptake = stringr::str_c("vax_uptake_", uptake_assumption)
    )
  if (uptake_assumption == "malik_all_ages") {
    lab_uptake_assumption <- "malik"  
  } else if (uptake_assumption == "census_all_ages") {
    lab_uptake_assumption <- "census"
  } else {
    lab_uptake_assumption <- uptake_assumption
  }
  
  sim_runs[[lab_uptake_assumption]] <- convert_to_seir_input_form( 
    compute_unreserved_nasem_with_uptake(
      pums_nasem_labeled,
      uptake_input_to_sim,
      size_of_increment = 100000
    )
  )
  write.csv(
    sim_runs[[lab_uptake_assumption]],
    stringr::str_c(
      build_data_dir,
      "/vaccination_flow_strategy_nasem_unreserved_with_",
      lab_uptake_assumption,
      "_uptake_param.csv"
    ),
    row.names = FALSE
  )
}