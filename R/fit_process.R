##' Process fit data
##'
##' @title Process a fit
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
##'
##' @export
spim_fit_process <- function(samples, parameters, data, control) {
  region <- samples$info$region

  message("Computing restart information")
  restart <- fit_process_restart(samples, parameters, data, control)
  samples$restart <- NULL

  message("Running forecasts")
  incidence_states <- "deaths"
  forecast <- sircovid::carehomes_forecast(samples,
                                           control$n_sample,
                                           control$burnin,
                                           control$forecast_days,
                                           incidence_states)

  message("Computing Rt")
  rt <- calculate_Rt(forecast, samples$info$multistrain, TRUE) # TODO: very slow
  if (samples$info$multistrain) {
    variant_rt <- calculate_Rt(forecast, samples$info$multistrain, FALSE)
  } else {
    variant_rt <- NULL
  }
  message("Computing IFR")
  ifr_t <-
    calculate_ifr_t(forecast, samples$info$multistrain) # TODO: a bit slow

  message("Summarising admissions")
  admissions <- extract_outputs_by_age(forecast, "cum_admit") # slow
  admissions[["data"]] <- data$admissions

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
##' @param data_admissions The admissions data set from
##'   [spimalot::spim_data_admissions]
##'
##' @param rtm The rtm data set
##'
##'
##' @param data The data set as passed to
##'   [spimalot::spim_particle_filter]
##'
##' @param data_full Full data set, before any right-censoring
##'
##' @param data_vaccination The vaccination data set as passed to
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

  ret <- c(ret, calculate_vaccination(ret$state, vaccine_efficacy, cross_immunity))

  # thin trajectories
  ret$state <- ret$state[c("deaths", "deaths_comm", "admitted", "diagnoses",
                           "infections", "hosp", "icu"), , ]

  # reshape to add a regional dimension
  ret$state <- mcstate::array_reshape(ret$state, i = 2, c(ncol(ret$state), 1))

  ret
}


calculate_Rt <- function(samples, multistrain, weight_Rt) {
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


calculate_ifr_t <- function(samples, multistrain) {
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

    ## aggregate partially immunised strata
    x[, 2L, , ] <- x[, 2L, , ] + x[, 3L, , ]
    x <- x[, -3L, , ]
    colnames(x) <- c("unvaccinated", "partial_protection", "full_protection")

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
    ret <- apply(abind::abind(res, along = 4), c(1, 3, 4), mean)

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
    n_groups <- de[[1]]
    n_strain <- 1
    n_vacc_classes <- de[[2]]
    multistrain <- FALSE
  } else {
    n_groups <- de[[1]]
    n_strain <- de[[2]]
    n_vacc_classes <- de[[3]]
    multistrain <- TRUE
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
        protected_against_severe_disease = sum_asr(c(vp$severe_disease[, 1, ]) * V),
        protected_against_death = sum_asr(c(vp$death[, 1, ]) * V),
        ever_infected = colSums(R_strain_1),
        ever_infected_unvaccinated = R_strain_1[1, ]),
      strain_2 = rbind(
        ever_vaccinated = ever_vaccinated,
        protected_against_infection = sum_asr(c(vp$infection[, 2, ]) * V),
        protected_against_severe_disease = sum_asr(c(vp$severe_disease[, 2, ]) * V),
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

  n_doses <- abind::abind(doses, doses_inc, along = 2)

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

  list(info = info,
       prior = prior,
       proposal = proposal)
}
