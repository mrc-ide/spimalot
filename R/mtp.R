spim_mtp_prepare <- function(mtp_commission, npi_key, n_par, end_date,
                             vaccine_parameters, combined, rt) {
  date <- combined$date
  npi_key <- mtp_npi_key_prepare(npi_key, rt, date)
  ## If date is NULL we get NaN entries for Rt and an obscure error
  ## later.

  ## TODO: might move 'rt' into the combined depending on size?

  ## expand data frame to include all regions
  f <- function(df) {
    msg <- setdiff(names(combined$pars), df$region)
    if (length(msg) > 0) {
      f <- function(r) {
        d <- df[df$region == "combined", ]
        d$region <- r
        d
      }
      df <- dplyr::bind_rows(c(list(df), lapply(msg, f)))
    }
    df[df$region != "combined", ]
  }

  mtp_commission <- split(mtp_commission, mtp_commission$scenario)
  ret <- lapply(mtp_commission, f) %>%
    dplyr::bind_rows()

  run_grid <- ret %>%
    dplyr::select(.data$scenario, .data$daily_doses, .data$spim_name,
                  .data$type_rt) %>%
    dplyr::distinct()

  rt_future <- ret %>%
    dplyr::select(-c(.data$daily_doses, .data$spim_name)) %>%
    tidyr::pivot_longer(c(tidyr::starts_with("npi"),
                          tidyr::starts_with("date")),
                        names_to = c(".value", "changepoint"),
                        names_pattern = "(.+)_(.+)") %>%
    dplyr::mutate(Rt_sd = npi_key[npi, "Rt_sd"],
                  Rt = npi_key[npi, "Rt"],
                  .after = npi)

  rt_future <- split(rt_future, f = rt_future$scenario)
  rt_future <- rt_future[run_grid$scenario]

  ## Things that we want to keep:
  inc_states <- c("deaths", "admitted", "diagnoses", "infections")
  prev_states <- c("icu", "general", "hosp", "react_pos")
  output <- list(incidence = inc_states,
                 prevalence = prev_states,
                 keep = c(inc_states, prev_states))

  vaccine <- list(
    future_daily_doses = lapply(vaccine_parameters$daily_doses, "[[", "mtp"),
    efficacy = vaccine_parameters$vacc_efficacy$central,
    uptake = vaccine_parameters$uptake_by_age$central)

  n_threads <- spim_control_cores()
  message(sprintf("Running on %d threads", n_threads))

  list(length = nrow(run_grid),
       combined = combined,
       n_threads = n_threads,
       n_par = n_par,
       end_date = end_date,
       run_grid = run_grid,
       rt_future = rt_future,
       vaccine = vaccine,
       output = output)
}


mtp_npi_key_prepare <- function(npi_key, rt, date) {
  curr_eff_Rt <- rt$eff_Rt_general[rt$date == sircovid::sircovid_date(date)]
  curr_Rt <- rt$Rt_general[rt$date == sircovid::sircovid_date(date)]
  npi_key$Rt <- round(npi_key$Reff_t * mean(curr_Rt / curr_eff_Rt), 1)
  npi_key["schools_main", "Rt"] <- round(mean(curr_Rt), 1)

  ## translate 0.5 reduction in Rt_excl_immunity for school closures
  ## into equivalent reduction in eff_Rt, this is likely to be around
  ## 0.3

  npi_key_holidays <- npi_key
  npi_key_holidays$Rt <- npi_key$Rt - 0.5
  rownames(npi_key_holidays) <- gsub("schools", "holidays", rownames(npi_key))

  rbind(npi_key, npi_key_holidays)
}


spim_mtp_simulate <- function(obj, subset = NULL) {
  n_run <- obj$length
  res <- vector("list", n_run)

  index <- subset %||% seq_len(n_run)

  i <- 1L
  for (i in index) {
    message(sprintf("Running scenario %d / %d", i, n_run))
    res[[i]] <- spim_simulate_schedule(
      combined = obj$combined,
      n_par = obj$n_par,
      end_date = obj$end_date,
      n_threads = obj$n_threads,
      keep = obj$output$keep,
      rt_future = obj$rt_future[[i]],
      rt_type = obj$run_grid$type_rt[i],
      future_daily_doses = obj$vaccine$future_daily_doses,
      vaccine_efficacy = obj$vaccine$efficacy,
      uptake_by_age = obj$vaccine$uptake,
      calculate_rt = TRUE,
      n_strain = 1)
  }
  names(res) <- obj$run_grid$scenario

  res
}
