spim_simulate <- function(i, args, combined) {
  el <- args[[i]]
  message(sprintf("-----\nRunning scenario %d / %d (%s [%s])", i, length(args),
                  el$analysis, el$scenario))
  time <- system.time(
    ret <- spim_simulate_one(el, combined))
  message(sprintf("Finished scenario %d in %2.1f s", i, time[["elapsed"]]))
  ret
}


simulation_central_analysis <- function(full_run = TRUE, multistrain = TRUE) {
  grid <- tibble::tibble(
    RUN = TRUE,
    full_run = full_run,
    analysis = "Central", adherence_to_baseline_npis = "central",
    seasonality = "central",
    vaccine_daily_doses = "co_central", vaccine_efficacy = "central",
    vaccine_uptake = "central",
    vaccine_eligibility = "min_18", vaccine_lag_groups = "no_delay",
    vaccine_lag_days = "no_delay", vaccine_booster_daily_doses = "no_booster",
    vaccine_booster_efficacy = "central"
   )

  if (multistrain) {
    grid <- cbind(
      grid,
      tibble::tibble(
        strain_transmission = "no_voc", strain_seed_rate = "no_seeding",
        strain_vaccine_efficacy = "central",
        strain_initial_proportion = "no_voc",
        strain_vaccine_booster_efficacy = "central"
      )
    )
  }

  grid
}
