##' Process fit data for `carehomes` model
##'
##' @title Process a fit for `carehomes` model
##' @param samples The pmcmc_samples object from [mcstate::pmcmc]
##'
##' @param parameters The parameter specification
##'   ([spimalot::spim_pars_pmcmc_load])
##'
##' @param data Data sets used in fitting, via
##'   [spimalot::spim_fit_process_data]
##'
##' @param control The forecast control from [spimalot::spim_control]
##'
##' @param random_sample Logical parameter, if `TRUE` will obtain the
##'   posterior samples via random sampling, otherwise thinning will
##'   be used
##'
##' @export
spim_fit_process <- function(samples, parameters, data, control,
                             random_sample = TRUE) {
  region <- samples$info$region

  message("Computing restart information")
  restart <- fit_process_restart(samples, parameters, data, control)
  samples$restart <- NULL

  message("Running forecasts")
  incidence_states <- "deaths"
  ## Add 1 to burnin to account for removal of initial parameters
  forecast <- sircovid::carehomes_forecast(samples,
                                           control$n_sample,
                                           control$burnin + 1L,
                                           control$forecast_days,
                                           incidence_states,
                                           random_sample = random_sample,
                                           thin = control$thin)

  message("Computing Rt")
  rt <- calculate_carehomes_Rt(forecast, samples$info$multistrain, TRUE)
  # TODO: very slow
  if (samples$info$multistrain) {
    variant_rt <-
      calculate_carehomes_Rt(forecast, samples$info$multistrain, FALSE)
  } else {
    variant_rt <- NULL
  }
  message("Computing IFR")
  ifr_t <-
    calculate_carehomes_ifr_t(forecast, samples$info$multistrain) # a bit slow

  if (is.null(data$admissions)) {
    admissions <- NULL
  } else {
    message("Summarising admissions")
    admissions <- extract_outputs_by_age(forecast, "cum_admit") # slow
    admissions[["data"]] <- data$admissions
  }

  message("Summarising deaths")
  deaths <- extract_outputs_by_age(forecast, "D_hosp") # slow
  i_deaths_data <- colnames(deaths$output_t)
  deaths$data <- data$rtm[data$rtm$region == region,
                          c("date", "region", i_deaths_data)]
  deaths$data[is.na(deaths$data)] <- 0

  ## TODO: someone needs to document what this date is for (appears to
  ## filter trajectories to start at this date) and when we might
  ## change it.
  message("Preparing onward simulation object")
  start_date_sim <- "2021-01-01"
  simulate <- create_simulate_object(
    forecast, samples$vaccine$efficacy, start_date_sim, samples$info$date)

  ## Reduce trajectories in forecast before saving
  message("Reducing trajectories")
  forecast <- reduce_trajectories(forecast)

  if (!is.null(restart)) {
    ## When adding the trajectories, we might as well strip them down
    ## to the last date in the restart
    i <- forecast$trajectories$date <= max(restart$state$time)
    restart$trajectories <- trajectories_filter_time(forecast$trajectories, i)
  }

  message("Computing parameter MLE and covariance matrix")
  parameters <- spim_fit_parameters(samples, parameters)

  if (!is.null(restart)) {
    ## When adding the trajectories, we might as well strip them down
    ## to the last date in the restart
    restart_date <- max(restart$state$time)
    i <- forecast$trajectories$date <= restart_date

    restart$parent <- list(
      trajectories = trajectories_filter_time(forecast$trajectories, i),
      rt = rt_filter_time(rt, i),
      ifr_t = rt_filter_time(ifr_t, i),
      deaths = deaths_filter_time(deaths, restart_date),
      admissions = deaths_filter_time(deaths, restart_date),
      ## TODO: check to make sure that this is just the one region's
      ## parameters at this point (see the region column)
      prior = parameters$prior)
  }

  ## Drop the big objects from the output
  samples[c("state", "trajectories", "predict")] <- list(NULL)

  list(samples = forecast, # note complicated naming change here
       pmcmc = samples,
       rt = rt,
       variant_rt = variant_rt,
       ifr_t = ifr_t,
       admissions = admissions,
       deaths = deaths,
       simulate = simulate,
       parameters = parameters,
       restart = restart,
       vaccination = data$vaccination,
       data = list(fitted = data$fitted, full = data$full))
}


##' Process fit data for `lancelot` model
##'
##' @title Process a fit for `lancelot` model
##' @param samples The pmcmc_samples object from [mcstate::pmcmc]
##'
##' @param parameters The parameter specification
##'   ([spimalot::spim_pars_pmcmc_load])
##'
##' @param data Data sets used in fitting, via
##'   [spimalot::spim_fit_process_data]
##'
##' @param control The forecast control from [spimalot::spim_control]
##'
##' @param random_sample Logical parameter, if `TRUE` will obtain the
##'   posterior samples via random sampling, otherwise thinning will
##'   be used
##'
##' @export
spim_lancelot_fit_process <- function(samples, parameters, data, control,
                                      random_sample = TRUE) {
  region <- samples$info$region

  message("Computing restart information")
  restart <- fit_process_restart(samples, parameters, data, control)
  samples$restart <- NULL

  message("Running forecasts")
  incidence_states <- "deaths"
  ## Add 1 to burnin to account for removal of initial parameters
  forecast <- sircovid::lancelot_forecast(samples,
                                          control$n_sample,
                                          control$burnin + 1L,
                                          control$forecast_days,
                                          incidence_states,
                                          random_sample = random_sample,
                                          thin = control$thin)

  message("Computing Rt")
  rt <- calculate_lancelot_Rt(forecast, samples$info$multistrain, TRUE)
  # TODO: very slow
  if (samples$info$multistrain) {
    variant_rt <-
      calculate_lancelot_Rt(forecast, samples$info$multistrain, FALSE)
  } else {
    variant_rt <- NULL
  }
  message("Computing IFR")
  ifr_t <-
    calculate_lancelot_ifr_t(forecast, samples$info$multistrain) # a bit slow

  if (is.null(data$admissions)) {
    admissions <- NULL
  } else {
    message("Summarising admissions")
    admissions <- extract_outputs_by_age(forecast, "cum_admit") # slow
    admissions[["data"]] <- data$admissions
  }

  message("Summarising deaths")
  deaths <- extract_outputs_by_age(forecast, "D_hosp") # slow
  i_deaths_data <- colnames(deaths$output_t)
  deaths$data <- data$rtm[data$rtm$region == region,
                          c("date", "region", i_deaths_data)]
  deaths$data[is.na(deaths$data)] <- 0

  ## TODO: someone needs to document what this date is for (appears to
  ## filter trajectories to start at this date) and when we might
  ## change it.
  message("Preparing onward simulation object")
  start_date_sim <- "2021-01-01"
  simulate <- create_simulate_object(
    forecast, samples$vaccine$efficacy, start_date_sim, samples$info$date)

  ## Reduce trajectories in forecast before saving
  message("Reducing trajectories")
  forecast <- reduce_trajectories(forecast)

  if (!is.null(restart)) {
    ## When adding the trajectories, we might as well strip them down
    ## to the last date in the restart
    i <- forecast$trajectories$date <= max(restart$state$time)
    restart$trajectories <- trajectories_filter_time(forecast$trajectories, i)
  }

  message("Computing parameter MLE and covariance matrix")
  parameters <- spim_fit_parameters(samples, parameters)

  if (!is.null(restart)) {
    ## When adding the trajectories, we might as well strip them down
    ## to the last date in the restart
    restart_date <- max(restart$state$time)
    i <- forecast$trajectories$date <= restart_date

    restart$parent <- list(
      trajectories = trajectories_filter_time(forecast$trajectories, i),
      rt = rt_filter_time(rt, i),
      ifr_t = rt_filter_time(ifr_t, i),
      deaths = deaths_filter_time(deaths, restart_date),
      admissions = deaths_filter_time(deaths, restart_date),
      ## TODO: check to make sure that this is just the one region's
      ## parameters at this point (see the region column)
      prior = parameters$prior)
  }

  ## Drop the big objects from the output
  samples[c("state", "trajectories", "predict")] <- list(NULL)

  list(samples = forecast, # note complicated naming change here
       pmcmc = samples,
       rt = rt,
       variant_rt = variant_rt,
       ifr_t = ifr_t,
       admissions = admissions,
       deaths = deaths,
       simulate = simulate,
       parameters = parameters,
       restart = restart,
       vaccination = data$vaccination,
       data = list(fitted = data$fitted, full = data$full))
}


##' Collect data sets for use with [spimalot::spim_fit_process]
##'
##' @title Collect data sets
##' @param admissions The admissions data set from
##'   [spimalot::spim_data_admissions]. Set to NULL if not fitting or plotting
##'   age-specific data
##'
##' @param rtm The rtm data set
##'
##' @param fitted The data set as passed to
##'   [spimalot::spim_particle_filter]
##'
##' @param full Full data set, before any right-censoring
##'
##' @param vaccination The vaccination data set as passed to
##'   [spimalot::spim_pars]
##'
##' @export
spim_fit_process_data <- function(admissions, rtm, fitted, full, vaccination) {
  list(admissions = admissions,
       rtm = rtm,
       full = full,
       fitted = fitted,
       vaccination = vaccination)
}


create_simulate_object <- function(samples, vaccine_efficacy, start_date_sim,
                                   date) {
  start_date_sim <- sircovid::sircovid_date(start_date_sim)
  fit_dates <- samples$trajectories$date
  idx_dates <- (fit_dates >= start_date_sim) &
    (fit_dates <= sircovid::sircovid_date(date))

  # trim dates to only those needed
  ret <- list(date = fit_dates[idx_dates],
              state = samples$trajectories$state[, , idx_dates])
  # add state_by_age
  ret$state_by_age <- extract_age_class_state(ret$state)
  # add n_protected and n_doses2s
  cross_immunity <- samples$predict$transform(samples$pars[1, ])$cross_immunity

  ret <-
    c(ret, calculate_vaccination(ret$state, vaccine_efficacy, cross_immunity))

  # thin trajectories
  ret$state <- ret$state[c("deaths", "deaths_comm", "deaths_hosp", "admitted",
                           "diagnoses", "infections", "hosp", "icu"), , ]

  # reshape to add a regional dimension
  ret$state <- mcstate::array_reshape(ret$state, i = 2, c(ncol(ret$state), 1))

  ret
}


calculate_carehomes_Rt <- function(samples, multistrain, weight_Rt) {
  step <- samples$trajectories$step

  index_S <- grep("^S_", names(samples$predict$index))
  index_R <- grep("^R_", names(samples$predict$index))
  index_ps <- grep("^prob_strain", names(samples$predict$index))

  if (multistrain) {
    R <- samples$trajectories$state[index_R, , , drop = FALSE]
    prob_strain <- samples$trajectories$state[index_ps, , , drop = FALSE]
  } else {
    R <- NULL
    prob_strain <- NULL
  }

  S <- samples$trajectories$state[index_S, , , drop = FALSE]

  pars <- lapply(seq_rows(samples$pars), function(i)
    samples$predict$transform(samples$pars[i, ]))

  sircovid::carehomes_Rt_trajectories(
    step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE, R = R, prob_strain = prob_strain,
    weight_Rt = weight_Rt)
}


calculate_lancelot_Rt <- function(samples, multistrain, weight_Rt) {
  step <- samples$trajectories$step

  index_S <- grep("^S_", names(samples$predict$index))
  index_R <- grep("^R_", names(samples$predict$index))
  index_ps <- grep("^prob_strain", names(samples$predict$index))

  if (multistrain) {
    R <- samples$trajectories$state[index_R, , , drop = FALSE]
    prob_strain <- samples$trajectories$state[index_ps, , , drop = FALSE]
  } else {
    R <- NULL
    prob_strain <- NULL
  }

  S <- samples$trajectories$state[index_S, , , drop = FALSE]

  pars <- lapply(seq_rows(samples$pars), function(i)
    samples$predict$transform(samples$pars[i, ]))

  sircovid::lancelot_Rt_trajectories(
    step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE, R = R, prob_strain = prob_strain,
    weight_Rt = weight_Rt)
}


calculate_carehomes_ifr_t <- function(samples, multistrain) {
  step <- samples$trajectories$step

  index_S <- grep("^S_", names(samples$predict$index))
  index_I_weighted <- grep("^I_weighted_", names(samples$predict$index))
  index_R <- grep("^R_", names(samples$predict$index))

  S <- samples$trajectories$state[index_S, , , drop = FALSE]
  I_weighted <- samples$trajectories$state[index_I_weighted, , , drop = FALSE]
  if (multistrain) {
    R <- samples$trajectories$state[index_R, , , drop = FALSE]
  } else {
    R <- NULL
  }

  pars <- lapply(seq_rows(samples$pars), function(i)
    samples$predict$transform(samples$pars[i, ]))

  sircovid::carehomes_ifr_t_trajectories(
    step, S, I_weighted, pars, R = R,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)
}


calculate_lancelot_ifr_t <- function(samples, multistrain) {
  step <- samples$trajectories$step

  index_S <- grep("^S_", names(samples$predict$index))
  index_I_weighted <- grep("^I_weighted_", names(samples$predict$index))
  index_R <- grep("^R_", names(samples$predict$index))

  S <- samples$trajectories$state[index_S, , , drop = FALSE]
  I_weighted <- samples$trajectories$state[index_I_weighted, , , drop = FALSE]
  if (multistrain) {
    R <- samples$trajectories$state[index_R, , , drop = FALSE]
  } else {
    R <- NULL
  }

  pars <- lapply(seq_rows(samples$pars), function(i)
    samples$predict$transform(samples$pars[i, ]))

  sircovid::lancelot_ifr_t_trajectories(
    step, S, I_weighted, pars, R = R,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)
}


## All the functions below here have awful names, but none are
## exported so we can tidy this up later.
extract_outputs_by_age <- function(sample, what) {
  trajectories <- sample$trajectories$state
  i <- grep(paste0("^", what), rownames(trajectories))
  cum_output <- trajectories[i, , , drop = FALSE]

  total_output <- cum_output[, , dim(cum_output)[3]]
  prop_output <- t(total_output) / colSums(total_output) * 100

  output <- apply(cum_output, 1:2, diff)
  mean_output <- apply(output, 1:2, mean)

  ## TODO: consider tdigest::tquantile (see quantile_digest in
  ## globals)
  out <- list(prop_total_output = prop_output,
              mean_prop_total_output = colMeans(prop_output),
              output_t = mean_output,
              lower_bound = apply(output, 1:2, quantile, 0.025),
              upper_bound = apply(output, 1:2, quantile, 0.975))

  out <- aggregate_outputs_by_age(out, what)

  if (diff(sample$trajectories$date[1:2]) != 1) {
    for (i in seq_along(out)) {
      out[[i]][1, ] <- NA
    }
  }

  out$date <- sample$trajectories$date[-1]

  out
}


aggregate_outputs_by_age <- function(object, what) {

  which_df <- c("output_t", "lower_bound", "upper_bound")
  out <- NULL
  if (what == "cum_admit") {
    for (df in which_df) {
      adm_0 <- rowSums(object[[df]][, 1:5], na.rm = TRUE)
      adm_25 <- rowSums(object[[df]][, 6:12], na.rm = TRUE) +
        (object[[df]][, 18] * 0.75)
      adm_55 <- rowSums(object[[df]][, 12:13], na.rm = TRUE) +
        (object[[df]][, 18] * 0.25)
      adm_65 <- rowSums(object[[df]][, 14:15], na.rm = TRUE) +
        (object[[df]][, 19] * 0.1)
      adm_75 <- rowSums(object[[df]][, 16:17], na.rm = TRUE) +
        (object[[df]][, 19] * 0.9)

      out[[df]] <- cbind(adm_0, adm_25, adm_55, adm_65, adm_75)
    }
  } else if (what == "D_hosp") {
    for (df in which_df) {

      death_0 <- rowSums(object[[df]][, 1:11], na.rm = TRUE) +
        (object[[df]][, 18] * 0.75)
      death_55 <- rowSums(object[[df]][, 12:13], na.rm = TRUE) +
        (object[[df]][, 18] * 0.25)
      death_65 <- rowSums(object[[df]][, 14:15], na.rm = TRUE)
      death_75 <- rowSums(object[[df]][, 16:17], na.rm = TRUE)
      death_chr <- object[[df]][, 19]

      out[[df]] <- cbind(death_0, death_55, death_65, death_75, death_chr)

    }
  } else {
    stop(sprintf("Can't aggregate '%s'", what))
  }
  out
}


extract_age_class_state <- function(state) {
  n_groups <- sircovid:::carehomes_n_groups()
  names_index <- c("cum_infections_disag", "diagnoses_admitted", "D_all")

  ## output cumulative states by
  ## age / vaccine class / sample / region / time
  arrays <- lapply(names_index, function(x) state[grep(x, rownames(state)), , ])
  names(arrays) <- c("infections", "diagnoses_admitted", "deaths")
  strata <- nrow(arrays$deaths) / n_groups

  f <- function(array) {
    x <- mcstate::array_reshape(array, 1L, c(n_groups, strata))

    if (ncol(x) == 5) {
      colnames(x) <- c("unvaccinated", "partial_protection", "full_protection",
                       "waned_protection", "booster_protection")
    } else {
      colnames(x) <- c("unvaccinated", "partial_protection", "full_protection",
                       "waned_protection")
    }


    ## aggregate age groups
    groups <- list(age_0 = 1:6, # 0-4, 5-9, 10-14, 15-19, 20-24, 25-29
                   age_30 = 7:10,  # 30-34, 35-39, 40-44, 45-49
                   age_50 = 11:15, # 50-54, 55-59, 60-64, 65-69, 70-74
                   age_75 = 16:17, # 75-79, 80+
                   chw = 18, chr = 19)

    res <- lapply(groups,
                  function(i) apply(x[i, , , , drop = FALSE], 2:4, sum))

    # distribute CHW between 30-49 and 50-74 age groups
    # distribute CHR between 50-74 and 75+ age groups
    res$age_30 <- res$age_30 + 0.75 * res$chw
    res$age_50 <- res$age_50 + 0.25 * res$chw + 0.1 * res$chr
    res$age_75 <- res$age_75 + 0.9 * res$chr
    res$chw <- NULL
    res$chr <- NULL

    # take mean across particles
    ret <- apply(abind_quiet(res, along = 4), c(1, 3, 4), mean)

    # [age, vaccine status, region, time]
    ret <- round(aperm(ret, c(3, 1, 2)))
    ret <- mcstate::array_reshape(ret, 3, c(1, dim(ret)[3]))

    ret
  }

  lapply(arrays, f)

}


reduce_trajectories <- function(samples) {
  ## Remove unused trajectories for predict function in combined
  remove_strings <- c("prob_strain", "^R_", "I_weighted_", "D_hosp_", "D_all_",
                      "diagnoses_admitted_", "cum_infections_disag_",
                      "cum_n_vaccinated")

  pars <- samples$predict$transform(samples$pars[1, ])
  n_groups <- pars$n_groups
  n_vacc_classes <- pars$n_vacc_classes

  state <- samples$trajectories$state

  index_remove <- lapply(remove_strings, function(s) {
    grep(paste0("^", s), rownames(state))
  })
  state <- state[-unlist(index_remove), , ]


  ## Add index_S
  nms_S <- grep("^S_", rownames(state), value = TRUE)
  S <- state[nms_S, , ]

  ## calculate the effective number of susceptibles
  calc_eff_S <- function(i) {
    p <- samples$predict$transform(samples$pars[i, ])
    rel_susceptibility <- p$rel_susceptibility
    apply(S[, i, ], 2, function(x) sum(x * c(rel_susceptibility)))
  }

  eff_S <- t(vapply(seq_len(dim(S)[2]), calc_eff_S, numeric(dim(S)[3])))
  eff_S <- array(eff_S, c(1, dim(eff_S)))
  row.names(eff_S) <- "eff_S"
  state <- abind1(state, eff_S)

  samples$trajectories$state <- abind1(samples$trajectories$state, eff_S)
  ## sum across vaccine stages for S
  if (n_groups != length(nms_S)) {
    ## sum across vacc classes
    S <- mcstate::array_reshape(S, i = 1, d = c(n_groups, n_vacc_classes))
    S <- apply(S, c(1, 3, 4), sum)
    rownames(S) <- paste0("S_", seq_len(n_groups))
  }

  samples$trajectories$state <-
    abind1(state[setdiff(rownames(state), nms_S), , ], S)

  ## Calculate Pillar 2 positivity and cases
  if (samples$info$model_type == "BB") {
    samples <- calculate_positivity(samples)
  } else {
    samples <- calculate_cases(samples)
  }

  samples
}


trajectories_filter_time <- function(trajectories, i) {
  trajectories$step <- trajectories$step[i]
  trajectories$date <- trajectories$date[i]
  trajectories$predicted <- trajectories$predicted[i]
  trajectories$state <- trajectories$state[, , i, drop = FALSE]
  trajectories
}


rt_filter_time <- function(rt, i) {
  ret <- lapply(rt, function(x) x[i, , drop = FALSE])
  class(ret) <- class(rt)
  ret
}


deaths_filter_time <- function(x, restart_date) {
  i <- x$date < restart_date
  for (v in c("output_t", "lower_bound", "upper_bound")) {
    x[[v]] <- x[[v]][i, , drop = FALSE]
  }
  x$date <- x$date[i]
  x$data <- x$data[sircovid::sircovid_date(x$data$date) < restart_date, ]
  x
}


calculate_vaccination <- function(state, vaccine_efficacy, cross_immunity) {

  de <- dim(vaccine_efficacy[[1]])
  if (length(de) == 2L) {
    ## We've changed the parameter preparation so that this branch
    ## should never be used
    message("WARNING: you are using old versions of tasks")
    n_groups <- de[[1]]
    n_strain <- 1
    n_vacc_classes <- de[[2]]
    multistrain <- FALSE
  } else {
    n_groups <- de[[1]]
    n_strain <- de[[2]]
    n_vacc_classes <- de[[3]]
    multistrain <- n_strain > 1
  }
  n_days <- dim(state)[[3]]

  # extract array of  mean by age / vaccine class / region == 1 / time
  get_mean_avt <- function(nm, state, strain = TRUE) {
    idx <- grep(nm, rownames(state))
    d <- c(n_groups, if (strain) n_strain else 1, n_vacc_classes)
    x <- mcstate::array_reshape(state[idx, , ], i = 1, d = d)
    out <- apply(x, c(1, 2, 3, 5), mean)
    ## TODO: this aperm should be removed, and code below here updated
    aperm(out, c(1, 3, 2, 4))
  }

  ## mean R by vaccine class / region == 1, time
  R <- get_mean_avt("^R_", state, TRUE)
  ## sum out age
  R <- apply(R, c(2, 3, 4), sum)

  ## mean cumulative vaccinations by age / vaccine class / region == 1 / time
  n_vaccinated <- get_mean_avt("^cum_n_vaccinated_", state, FALSE)

  vp <- get_vaccine_protection(vaccine_efficacy)

  # Methodology: calculate incidence of first / second doses,
  # number in each strata in total,
  # number in each strata who are recovered, use these to calculate proportion
  # protected as shown in main fig A / B.
  # to calculate proportion protected we need
  # s = stratum, r = region, t = time
  # let V be number in each vaccination stage V[a, s, r, t]
  # ler R be number recovered R[s, r, t]
  # from top to bottom of figure 1B
  # - number vaccinated sum_sr(V[s, r, t])
  # - number protected after vaccination against severe disease:
  #   i.e., sum_asr(V[a, s, r, t] * eff_severe[a, s])
  # - number protected after vaccination against infection
  #   i.e., sum_asr(V[a, s, r, t] * eff_inf[a, s])
  # - number protected after infection (this is added to all of the above):
  #   i.e., sum_sr(R[s, r, t])
  # - number protected after infection only:
  #   i.e. sum_r(R[1, r, t]


  ## create array of number in each vaccine stratum
  V <- array(0, dim = dim(n_vaccinated))
  V[, n_vacc_classes, , ] <- n_vaccinated[, n_vacc_classes - 1L, , ]
  for (i in seq(2, n_vacc_classes - 1)) {
    V[, i, , ] <- n_vaccinated[, i - 1, , ] - n_vaccinated[, i, , ]
  }

  sum_sr <- function(x) c(apply(x, c(2, 3), sum))
  sum_asr <- function(x) c(apply(x, c(3, 4), sum))

  if (multistrain) {

    ever_vaccinated <- colSums(n_vaccinated[, 1, , ])
    R_strain_1 <- apply(R[, -2, ], c(1, 3), sum) + R[, 2, ] * cross_immunity[2]
    R_strain_2 <-  apply(R[, -1, ], c(1, 3), sum) + R[, 1, ] * cross_immunity[1]

    n_protected <- list(
      strain_1 = rbind(
        ever_vaccinated = ever_vaccinated,
        protected_against_infection = sum_asr(c(vp$infection[, 1, ]) * V),
        protected_against_severe_disease =
          sum_asr(c(vp$severe_disease[, 1, ]) * V),
        protected_against_death = sum_asr(c(vp$death[, 1, ]) * V),
        ever_infected = colSums(R_strain_1),
        ever_infected_unvaccinated = R_strain_1[1, ]),
      strain_2 = rbind(
        ever_vaccinated = ever_vaccinated,
        protected_against_infection = sum_asr(c(vp$infection[, 2, ]) * V),
        protected_against_severe_disease =
          sum_asr(c(vp$severe_disease[, 2, ]) * V),
        protected_against_death = sum_asr(c(vp$death[, 2, ]) * V),
        ever_infected = colSums(R_strain_2),
        ever_infected_unvaccinated = R_strain_2[1, ]))

  } else {

    n_protected <- list(strain_1 = rbind(
      ever_vaccinated = colSums(n_vaccinated[, 1, , ]),
      protected_against_infection = sum_asr(c(vp$infection) * V),
      protected_against_severe_disease = sum_asr(c(vp$severe_disease) * V),
      protected_against_death = sum_asr(c(vp$death) * V),
      ever_infected = sum_sr(R),
      ever_infected_unvaccinated = R[1, , , drop = FALSE]))
  }


  ## calculate n_doses

  # Output number of first and second doses
  doses <- n_vaccinated[, c(1, 3), , , drop = FALSE]
  colnames(doses) <- c("first_dose", "second_dose")
  doses_inc <- aperm(apply(doses, c(1, 2, 3), diff), c(2, 3, 4, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-4], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")

  n_doses <- abind_quiet(doses, doses_inc, along = 2)

  list(n_protected = lapply(n_protected, mcstate::array_reshape, i = 2,
                            d = c(1, ncol(n_protected[[1]]))),
       n_doses = n_doses)
}


get_vaccine_protection <- function(vaccine_efficacy, booster_efficacy = NULL) {
  if (!is.null(booster_efficacy)) {
    stopifnot(identical(names(vaccine_efficacy), names(booster_efficacy)))
    vaccine_efficacy <- Map(cbind, vaccine_efficacy, booster_efficacy)
  }

  efficacy_infection <- 1 - vaccine_efficacy$rel_susceptibility
  efficacy_disease <- efficacy_infection + (1 - efficacy_infection) *
    (1 - vaccine_efficacy$rel_p_sympt)
  efficacy_severe_disease <- efficacy_disease + (1 - efficacy_disease) *
    (1 - vaccine_efficacy$rel_p_hosp_if_sympt)
  efficacy_death <- efficacy_severe_disease + (1 - efficacy_severe_disease) *
    (1 - vaccine_efficacy$rel_p_death)

  list(infection = efficacy_infection,
       disease = efficacy_disease,
       severe_disease = efficacy_severe_disease,
       death = efficacy_death)
}


extract_age_class_outputs <- function(samples) {
  ## Get index of states with age/vacc outputs
  names_index <- c("cum_infections_disag", "diagnoses_admitted", "D_all")
  index <- grep(paste(names_index, collapse = "|"),
                dimnames(samples$trajectories$state)[[1]])

  ## Save object for the age_vacc_class states
  if (length(dim(samples$trajectories$state)) == 4) {
    samples$trajectories$state[index, , , ]
  } else {
    samples$trajectories$state[index, , ]
  }
}


spim_fit_parameters <- function(samples, parameters) {
  info <- parameters$info[parameters$info$region == samples$info$region, ]
  rownames(info) <- NULL
  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]
  info$initial[match(names(initial), info$name)] <- unname(initial)

  prior <- parameters$prior[parameters$prior$region == samples$info$region, ]
  rownames(prior) <- NULL

  covariance <- cov(samples$pars)
  rownames(covariance) <- NULL
  proposal <- data_frame(region = samples$info$region,
                         name = colnames(covariance),
                         covariance)

  parameters$info <- info
  parameters$prior <- prior
  parameters$proposal <- proposal

  parameters
}


calculate_positivity <- function(samples) {

  p_NC_names <- c("p_NC_under15", "p_NC_15_24", "p_NC_25_49",
                  "p_NC_50_64", "p_NC_65_79", "p_NC_80_plus",
                  "p_NC_weekend_under15", "p_NC_weekend_15_24",
                  "p_NC_weekend_25_49", "p_NC_weekend_50_64",
                  "p_NC_weekend_65_79", "p_NC_weekend_80_plus")

  x <- sircovid::sircovid_date_as_date(samples$trajectories$date)

  base_pars <- samples$predict$transform(samples$pars[1, ])

  pars <- t(vapply(seq_len(nrow(samples$pars)),
                 function(i) unlist(samples$predict$transform(
                   samples$pars[i, ])[p_NC_names]),
                 numeric(length(p_NC_names))))

  pos_under15 <-
    samples$trajectories$state[paste0("sympt_cases_under15_inc"), , ]
  pos_15_24 <- samples$trajectories$state[paste0("sympt_cases_15_24_inc"), , ]
  pos_25_49 <- samples$trajectories$state[paste0("sympt_cases_25_49_inc"), , ]
  pos_50_64 <- samples$trajectories$state[paste0("sympt_cases_50_64_inc"), , ]
  pos_65_79 <- samples$trajectories$state[paste0("sympt_cases_65_79_inc"), , ]
  pos_80_plus <-
    samples$trajectories$state[paste0("sympt_cases_80_plus_inc"), , ]

  pos_over25 <- pos_25_49 + pos_50_64 + pos_65_79 + pos_80_plus
  pos_all <- pos_under15 + pos_15_24 + pos_over25

  calc_negs <- function(group) {
    neg <- base_pars[[paste0("N_tot_", group)]] -
      samples$trajectories$state[paste0("sympt_cases_", group, "_inc"), , ]

    neg[, grepl("^S", weekdays(x))] <-
      neg[, grepl("^S", weekdays(x))] * pars[, paste0("p_NC_weekend_", group)]
    neg[, !grepl("^S", weekdays(x))] <-
      neg[, !grepl("^S", weekdays(x))] * pars[, paste0("p_NC_", group)]

    neg
  }

  neg_under15 <- calc_negs("under15")
  neg_15_24 <- calc_negs("15_24")
  neg_25_49 <- calc_negs("25_49")
  neg_50_64 <- calc_negs("50_64")
  neg_65_79 <- calc_negs("65_79")
  neg_80_plus <- calc_negs("80_plus")

  neg_over25 <- neg_25_49 + neg_50_64 + neg_65_79 + neg_80_plus
  neg_all <- neg_under15 + neg_15_24 + neg_over25

  calc_pos <- function(pos, neg) {
    positivity1 <-
      (pos * base_pars$pillar2_sensitivity +
         neg * (1 - base_pars$pillar2_specificity)) / (pos + neg) * 100
    array(positivity1, c(1, dim(positivity1)))
  }

  positivity <- calc_pos(pos_all, neg_all)
  positivity <- abind1(positivity, calc_pos(pos_over25, neg_over25))
  positivity <- abind1(positivity, calc_pos(pos_under15, neg_under15))
  positivity <- abind1(positivity, calc_pos(pos_15_24, neg_15_24))
  positivity <- abind1(positivity, calc_pos(pos_25_49, neg_25_49))
  positivity <- abind1(positivity, calc_pos(pos_50_64, neg_50_64))
  positivity <- abind1(positivity, calc_pos(pos_65_79, neg_65_79))
  positivity <- abind1(positivity, calc_pos(pos_80_plus, neg_80_plus))

  row.names(positivity) <- paste0("pillar2_positivity",
                                  c("", "_over25", "_under15", "_15_24",
                                    "_25_49", "_50_64", "_65_79", "_80_plus"))

  samples$trajectories$state <-
    abind1(samples$trajectories$state, positivity)


  samples
}


calculate_cases <- function(samples) {

  phi_pillar2_cases_names <-
    c("phi_pillar2_cases_under15", "phi_pillar2_cases_15_24",
      "phi_pillar2_cases_25_49", "phi_pillar2_cases_50_64",
      "phi_pillar2_cases_65_79", "phi_pillar2_cases_80_plus",
      "phi_pillar2_cases_weekend_under15", "phi_pillar2_cases_weekend_15_24",
      "phi_pillar2_cases_weekend_25_49", "phi_pillar2_cases_weekend_50_64",
      "phi_pillar2_cases_weekend_65_79", "phi_pillar2_cases_weekend_80_plus")

  x <- sircovid::sircovid_date_as_date(samples$trajectories$date)

  base_pars <- samples$predict$transform(samples$pars[1, ])

  pars <- t(vapply(seq_len(nrow(samples$pars)),
                   function(i) unlist(samples$predict$transform(
                     samples$pars[i, ])[phi_pillar2_cases_names]),
                   numeric(length(phi_pillar2_cases_names))))

  calc_cases <- function(group) {
    cases1 <-
      samples$trajectories$state[paste0("sympt_cases_", group, "_inc"), , ]

    cases1[, grepl("^S", weekdays(x))] <-
      cases1[, grepl("^S", weekdays(x))] *
      pars[, paste0("phi_pillar2_cases_weekend_", group)]
    cases1[, !grepl("^S", weekdays(x))] <-
      cases1[, !grepl("^S", weekdays(x))] *
      pars[, paste0("phi_pillar2_cases_", group)]

    array(cases1, c(1, dim(cases1)))
  }

  cases_under15 <- calc_cases("under15")
  cases_15_24 <- calc_cases("15_24")
  cases_25_49 <- calc_cases("25_49")
  cases_50_64 <- calc_cases("50_64")
  cases_65_79 <- calc_cases("65_79")
  cases_80_plus <- calc_cases("80_plus")

  cases_over25 <- cases_25_49 + cases_50_64 + cases_65_79 + cases_80_plus
  cases_all <- cases_under15 + cases_15_24 + cases_over25

  cases <- abind1(cases_all, cases_over25)
  cases <- abind1(cases, cases_under15)
  cases <- abind1(cases, cases_15_24)
  cases <- abind1(cases, cases_25_49)
  cases <- abind1(cases, cases_50_64)
  cases <- abind1(cases, cases_65_79)
  cases <- abind1(cases, cases_80_plus)

  row.names(cases) <- paste0("pillar2_cases",
                             c("", "_over25", "_under15", "_15_24",
                               "_25_49", "_50_64", "_65_79", "_80_plus"))

  samples$trajectories$state <-
    abind1(samples$trajectories$state, cases)

  samples
}
