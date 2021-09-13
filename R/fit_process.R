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
##' @param random_sample Logical parameter, if `TRUE` will obtain the
##'   posterior samples via random sampling, otherwise thinning will
##'   be used
##'
##' @export
spim_fit_process <- function(samples, parameters, data, control,
                             random_sample = TRUE) {

  region <- samples$info$region
  single_region <- length(region) == 1

  ## inflate to multiregion and store original to return later
  original_samples <- samples
  samples <- inflate_samples(samples)

  message("Computing restart information")
  restart <- fit_process_restart(samples, parameters, data, control)
  samples$restart <- NULL

  message("Running forecasts")
  incidence_states <- "deaths"

  ## Add 1 to burnin to account for removal of initial parameters
  ## This will automatically work with original (non-inflated) samples in
  ##  single and multi-region cases
  forecast <- sircovid::carehomes_forecast(original_samples,
                                           control$n_sample,
                                           control$burnin + 1L,
                                           control$forecast_days,
                                           incidence_states,
                                           random_sample = random_sample,
                                           thin = control$thin)

  ## Inflate forecast in same way as samples
  original_forecast <- forecast
  forecast <- inflate_samples(forecast)

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
  deaths$data <- data$rtm[data$rtm$region %in% region,
                          c("date", "region", i_deaths_data)]
  deaths$data[is.na(deaths$data)] <- 0

  ## TODO: someone needs to document what this date is for (appears to
  ## filter trajectories to start at this date) and when we might
  ## change it.
  message("Preparing onward simulation object")
  start_date_sim <- "2021-01-01"

  simulate <- create_simulate_object(
    forecast, samples$vaccine[[1]]$efficacy, start_date_sim, samples$info$date
  )

  message("Computing parameter MLE and covariance matrix")
  parameters <- spim_fit_parameters(samples, parameters)

  ## Reduce trajectories in forecast before saving
  ## Note we revert to non-inflated
  message("Reducing trajectories")
  forecast <- deflate_samples(reduce_trajectories(forecast))

  if (!is.null(restart)) {
    ## When adding the trajectories, we might as well strip them down
    ## to the last date in the restart
    restart_date <- max(restart$state$time)
    i <- forecast$trajectories$date <= restart_date

    if (single_region) {
      restart_rt <- rt_filter_time(rt, i)
      restart_ifr_t <- rt_filter_time(ifr_t, i)
    } else {
      ## FIXME - Is a list really the best method?
      restart_rt <- lapply(rt, rt_filter_time, i)
      restart_ifr_t <- lapply(ifr_t, rt_filter_time, i)
    }

    restart$parent <- list(
      ## this implies restart$parent$trajectories = restart$trajectories
      trajectories = trajectories_filter_time(forecast$trajectories, i),
      rt = restart_rt,
      ifr_t = restart_ifr_t,
      deaths = deaths_filter_time(deaths, restart_date),
      admissions = deaths_filter_time(deaths, restart_date),
      prior = parameters$prior)
  }

  ## Drop the big objects from the output
  ##  (note not using inflated samples here)
  original_samples[c("state", "trajectories", "predict")] <- list(NULL)

  list(samples = forecast, # note complicated naming change here
       pmcmc = original_samples,
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

  regions <- samples$info$region

  # trim dates to only those needed
  ret <- list(date = fit_dates[idx_dates])

  state <- samples$trajectories$state[, , , idx_dates, drop = FALSE]
  if (length(idx_dates) == 1) {
    dim(state) <- dim(state)[1:3]
  }

  # add state_by_age
  ret$state_by_age <- extract_age_class_state(state, regions)
  # add n_protected and n_doses2s
  cross_immunity <- samples$predict$transform(t(samples$pars[, , 1]))
  if (length(regions) > 1) {
    cross_immunity <- cross_immunity[[1]][["cross_immunity"]]
  } else {
    cross_immunity <- cross_immunity[["cross_immunity"]]
  }

  ret <-
    c(ret, calculate_vaccination(state, vaccine_efficacy, cross_immunity,
                                 regions))

  # thin trajectories
  state <- state[c("deaths", "deaths_comm", "deaths_hosp", "admitted",
                   "diagnoses", "infections", "hosp", "icu"), , , ,
                   drop = FALSE]

  ret$state <- state
  ret
}


calculate_Rt_ifr_t <- function(samples, multistrain, weight_Rt, what) {
  step <- samples$trajectories$step

  index_S <- grep("^S_", names(samples$predict$index))
  S <- samples$trajectories$state[index_S, , , , drop = FALSE]

  index_I_weighted <- grep("^I_weighted_", names(samples$predict$index))
  I_weighted <- samples$trajectories$state[index_I_weighted, , , ,
                                           drop = FALSE]

  if (multistrain) {
    index_ps <- grep("^prob_strain", names(samples$predict$index))
    index_R <- grep("^R_", names(samples$predict$index))
    R <- samples$trajectories$state[index_R, , , , drop = FALSE]
    prob_strain <- samples$trajectories$state[index_ps, , , , drop = FALSE]
  } else {
    R <- NULL
    prob_strain <- NULL
  }

  region <- samples$info$region
  single_region <- length(region) == 1

  pars <- lapply(seq_len(nlayers(samples$pars)), function(i) {
    out <- samples$predict$transform(t(samples$pars[, , i]))
    if (!single_region) {
      names(out) <- region
    }
    out
  })

  if (single_region) {
    pars <- list(pars)
  } else {
    pars <- list_transpose(pars)
  }


  ret <- lapply(seq_along(pars), function(i) {
    if (what == "rt") {
      sircovid::carehomes_Rt_trajectories(
        step, S[, , i, ], pars[[i]],
        initial_step_from_parameters = TRUE,
        shared_parameters = FALSE, R = R[, , i, ],
        prob_strain = prob_strain[, , i, ],
        weight_Rt = weight_Rt
      )
    } else if (what == "ifr_t") {
      sircovid::carehomes_ifr_t_trajectories(
        step, S[, , i, ], I_weighted[, , i, ], pars[[i]],
        initial_step_from_parameters = TRUE,
        shared_parameters = FALSE, R = R[, , i, ]
      )
    } else {
      stop(sprintf("'what' should be 'rt' or 'ifr_t', but got '%s'", what))
    }
  })

  if (single_region) {
    ret <- ret[[1]]
  } else {
    names(ret) <- region
  }

  ret
}


calculate_Rt <- function(samples, multistrain, weight_Rt) {
  calculate_Rt_ifr_t(samples, multistrain, weight_Rt, "rt")
}


calculate_ifr_t <- function(samples, multistrain) {
  calculate_Rt_ifr_t(samples, multistrain, NULL, "ifr_t")
}


## All the functions below here have awful names, but none are
## exported so we can tidy this up later.
extract_outputs_by_age <- function(sample, what) {

  trajectories <- sample$trajectories$state

  i <- grep(paste0("^", what), rownames(trajectories))
  cum_output <- trajectories[i, , , , drop = FALSE]

  total_output <- cum_output[, , , dim(cum_output)[4], drop = FALSE]
  dim(total_output) <- dim(total_output)[1:3]
  ## TODO - Confirm correct dimensions used here
  prop_output <- mcstate::array_reshape(
    apply(total_output, 3, function(x) t(proportions(x, 2)) * 100),
    1, c(ncol(total_output), nrow(total_output))
  )

  ## TODO - Confirm correct dimensions used here
  output <- apply(cum_output, 1:3, diff)
  mean_output <- apply(output, c(1, 2, 4), mean)

  ## TODO: consider tdigest::tquantile (see quantile_digest in
  ## globals)
  ## FIXME: Removed `prop_total_output` and `mean_prop_total_output`.
  ##  Someone please confirm
  out <- list(output_t = mean_output,
              lower_bound = apply(output, c(1:2, 4), quantile, 0.025),
              upper_bound = apply(output, c(1:2, 4), quantile, 0.975))

  out <- aggregate_outputs_by_age(out, what, sample$info$region)

  if (diff(sample$trajectories$date[1:2]) != 1) {
    for (i in seq_along(out)) {
      out[[i]][1, , ] <- NA
    }
  }

  out$date <- sample$trajectories$date[-1]

  out
}


aggregate_outputs_by_age <- function(object, what, regions) {

  which_df <- c("output_t", "lower_bound", "upper_bound")
  out <- NULL
  if (what == "cum_admit") {
    for (df in which_df) {
      adm_0 <- apply(object[[df]][, 1:5, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE)
      adm_25 <- apply(object[[df]][, 6:12, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE) + (object[[df]][, 18, ] * 0.75)
      adm_55 <- apply(object[[df]][, 12:13, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE) + (object[[df]][, 18, ] * 0.25)
      adm_65 <- apply(object[[df]][, 14:15, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE) + (object[[df]][, 19, ] * 0.1)
      adm_75 <- apply(object[[df]][, 16:17, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE) + (object[[df]][, 19, ] * 0.9)

      out[[df]] <- cbind(adm_0, adm_25, adm_55, adm_65, adm_75)
      dim(out[[df]]) <- c(nrow(out[[df]]), 5, nlayers(object[[df]]))
      dimnames(out[[df]]) <- list(NULL, paste0("adm_", c(0, 25, 55, 65, 75)),
                                  regions)

    }
  } else if (what == "D_hosp") {
    for (df in which_df) {

      death_0 <- apply(object[[df]][, 1:11, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE) + (object[[df]][, 18, ] * 0.75)
      death_55 <- apply(object[[df]][, 12:13, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE) + (object[[df]][, 18, ] * 0.25)
      death_65 <- apply(object[[df]][, 14:15, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE)
      death_75 <- apply(object[[df]][, 16:17, , drop = FALSE], 3, rowSums,
                      na.rm = TRUE)
      death_chr <- object[[df]][, 19, ]

      out[[df]] <- cbind(death_0, death_55, death_65, death_75, death_chr)
      dim(out[[df]]) <- c(nrow(out[[df]]), 5, nlayers(object[[df]]))
      dimnames(out[[df]]) <- list(NULL,
                                  paste0("death_", c(0, 55, 65, 75, "chr")),
                                  regions)
    }
  } else {
    stop(sprintf("Can't aggregate '%s'", what))
  }
  out
}


extract_age_class_state <- function(state, region) {

  n_groups <- sircovid:::carehomes_n_groups()
  names_index <- c("cum_infections_disag", "diagnoses_admitted", "D_all")

  ## output cumulative states by
  ## age / vaccine class / sample / region / time
  arrays <- lapply(names_index,
                   function(x) state[grep(x, rownames(state)), , , ,
                                     drop = FALSE])
  names(arrays) <- c("infections", "diagnoses_admitted", "deaths")
  strata <- nrow(arrays$deaths) / n_groups

  f <- function(array) {

    x <- mcstate::array_reshape(array, 1L, c(n_groups, strata))

    ## aggregate partially immunised strata
    ## need to preserve dimensions but can't add on unit dimension
    x2 <- x
    x[, 2L, , , ] <- x[, 2L, , , ] + x[, 3L, , , ]
    x <- x[, -3L, , , ]
    if (length(dim(x)) != length(dim(x2))) {
      x <- add_dimension(x, 4)
      stopifnot(all(dim(x) == c(n_groups, 3, nlayers(x2), 1, dim(x2)[[5]])))
    }

    colnames(x) <- c("unvaccinated", "partial_protection", "full_protection")

    ## aggregate age groups
    groups <- list(age_0 = 1:6, # 0-4, 5-9, 10-14, 15-19, 20-24, 25-29
                   age_30 = 7:10,  # 30-34, 35-39, 40-44, 45-49
                   age_50 = 11:15, # 50-54, 55-59, 60-64, 65-69, 70-74
                   age_75 = 16:17, # 75-79, 80+
                   chw = 18, chr = 19)

    res <- lapply(groups,
                  function(i) apply(x[i, , , , , drop = FALSE], 2:5, sum))

    # distribute CHW between 30-49 and 50-74 age groups
    # distribute CHR between 50-74 and 75+ age groups
    res$age_30 <- res$age_30 + 0.75 * res$chw
    res$age_50 <- res$age_50 + 0.25 * res$chw + 0.1 * res$chr
    res$age_75 <- res$age_75 + 0.9 * res$chr
    res$chw <- NULL
    res$chr <- NULL

    # take mean across particles
    ret <- apply(abind::abind(res, along = 5), c(1, 3:5), mean)

    # [age, vaccine status, region, time]
    ret <- round(aperm(ret, c(4, 1, 2, 3)))
    dimnames(ret)[3] <- list(region)

    ret
  }

  lapply(arrays, f)

}


reduce_trajectories <- function(samples) {
  single_region <- length(samples$info$region) == 1
  ## Remove unused trajectories for predict function in combined
  remove_strings <- c("prob_strain", "^R_", "I_weighted_", "D_hosp_", "D_all_",
                      "diagnoses_admitted_", "cum_infections_disag_",
                      "cum_n_vaccinated")

  pars <- samples$predict$transform(t(samples$pars[, , 1]))
  if (!single_region) {
    pars <- pars[[1]]
  }

  n_groups <- pars$n_groups
  n_vacc_classes <- pars$n_vacc_classes

  state <- samples$trajectories$state

  index_remove <- lapply(remove_strings, function(s) {
    grep(paste0("^", s), rownames(state))
  })
  state <- state[-unlist(index_remove), , , , drop = FALSE]


  ## Add index_S
  nms_S <- grep("^S_", rownames(state), value = TRUE)
  S <- state[nms_S, , , , drop = FALSE]

  ## calculate the effective number of susceptibles
  calc_eff_S <- function(i) {

    p <- samples$predict$transform(t(samples$pars[, , i]))

    if (single_region) {
      ## same for all regions
      rel_susceptibility <- p$rel_susceptibility
    } else {
      rel_susceptibility <- p[[1]]$rel_susceptibility
    }

    apply(drop_dimension(S[, i, , , drop = FALSE], 2), 2:3,
          function(x) sum(x * c(rel_susceptibility)))
  }

  eff_S <- vapply(seq_len(dim(S)[2]), calc_eff_S,
                  matrix(0, dim(S)[3L], dim(S)[4L]))
  eff_S <- aperm(eff_S, c(3, 1, 2))
  eff_S <- array(eff_S, c(1, dim(eff_S)))
  row.names(eff_S) <- "eff_S"
  state <- abind::abind(state, eff_S, along = 1)

  samples$trajectories$state <- abind::abind(samples$trajectories$state,
                                             eff_S, along = 1)
  ## sum across vaccine stages for S
  if (n_groups != length(nms_S)) {
    ## sum across vacc classes
    S <- mcstate::array_reshape(S, i = 1, d = c(n_groups, n_vacc_classes))
    S <- apply(S, c(1, 3, 4, 5), sum)
    rownames(S) <- paste0("S_", seq_len(n_groups))
  }

  samples$trajectories$state <-
    abind::abind(state[setdiff(rownames(state), nms_S), , , , drop = FALSE], S,
                 along = 1)

  ## Calculate Pillar 2 positivity
  pillar2_positivity <- calculate_positivity(samples, FALSE)
  pillar2_positivity <-
    array(pillar2_positivity, c(1, dim(pillar2_positivity)))
  pillar2_positivity_over25 <- calculate_positivity(samples, TRUE)
  pillar2_positivity_over25 <-
    array(pillar2_positivity_over25, c(1, dim(pillar2_positivity_over25)))
  pillar2_positivity <- abind::abind(pillar2_positivity,
                                     pillar2_positivity_over25, along = 1)
  row.names(pillar2_positivity) <-
    c("pillar2_positivity", "pillar2_positivity_over25")
  samples$trajectories$state <- abind::abind(samples$trajectories$state,
                                             pillar2_positivity, along = 1)

  samples
}


trajectories_filter_time <- function(trajectories, i) {
  trajectories$step <- trajectories$step[i]
  trajectories$date <- trajectories$date[i]
  trajectories$predicted <- trajectories$predicted[i]
  if (length(dim(trajectories$state)) == 4) {
    trajectories$state <- trajectories$state[, , , i, drop = FALSE]
  } else {
    trajectories$state <- trajectories$state[, , i, drop = FALSE]
  }

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
    x[[v]] <- x[[v]][i, , , drop = FALSE]
  }
  x$date <- x$date[i]
  x$data <- x$data[sircovid::sircovid_date(x$data$date) < restart_date, ]
  x
}


calculate_vaccination <- function(state, vaccine_efficacy, cross_immunity,
                                  regions) {

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
  n_days <- dim(state)[[4]]

  # extract array of mean by age / strain / vacc class / regions / time
  get_mean_avt <- function(nm, state, strain = TRUE) {
    idx <- grep(nm, rownames(state))
    d <- c(n_groups, if (strain) n_strain else 1, n_vacc_classes)
    x <- mcstate::array_reshape(state[idx, , , , drop = FALSE], i = 1, d = d)
    apply(x, c(1, 2, 3, 5, 6), mean)
    ## TODO: this aperm should be removed, and code below here updated
    # aperm(out, c(1, 3, 2, 4))
  }

  ## mean R
  R <- get_mean_avt("^R_", state, TRUE)
  ## sum out age - strain / vacc class / region / time
  R <- apply(R, 2:5, sum)

  ## mean cumulative vaccinations
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
  V[, , n_vacc_classes, , ] <- n_vaccinated[, , n_vacc_classes - 1L, , ]
  for (i in seq(2, n_vacc_classes - 1)) {
    V[, , i, , ] <- n_vaccinated[, , i - 1, , ] - n_vaccinated[, , i, , ]
  }

  new_dim <- c(1, length(regions), dim(V)[5L])
  ## don't differentiate between strains
  ever_vaccinated <- array(apply(n_vaccinated[, , 1, , , drop = FALSE],
                                   4:5, sum), new_dim)


  calculate_n_protected <- function(ever, inf, sd, death, R) {
    sum_asr <- function(x) apply(c(x) * V, c(2, 4, 5), sum)

    ret <- abind::abind(
      ever_vaccinated = ever,
      ## TODO (LW): Confirm calculation below still correct
      protected_against_infection = sum_asr(inf),
      protected_against_severe_disease = sum_asr(sd),
      protected_against_death = sum_asr(death),
      ## sum over vaccine class (everyone infected indep. vacc class)
      ever_infected = apply(R, 2:3, sum),
      ## everyone recovered in vacc class 1
      ever_infected_unvaccinated = array(R[1, , , drop = FALSE], new_dim),
      along = 1
    )
    colnames(ret) <- regions
    ret
  }

  if (multistrain) {
    R_strain_1 <- apply(R[-2, , , , drop = FALSE], 2:4, sum) +
      drop_dimension(R[2, , , , drop = FALSE], 1) * cross_immunity[2]
    R_strain_2 <- apply(R[-1, , , , drop = FALSE], 2:4, sum) +
      drop_dimension(R[1, , , , drop = FALSE], 1) * cross_immunity[1]

    n_protected <- list(
      strain_1 = calculate_n_protected(
        ever_vaccinated, vp$infection[, 1, ], vp$severe_disease[, 1, ],
        vp$death[, 1, ], R_strain_1),
      strain_2 = calculate_n_protected(
        ever_vaccinated, vp$infection[, 2, ], vp$severe_disease[, 2, ],
        vp$death[, 2, ], R_strain_2)
    )

  } else {
    n_protected <- list(
      strain_1 = calculate_n_protected(
        ever_vaccinated, vp$infection, vp$severe_disease, vp$death,
        drop_dimension(R[1, , , , drop = FALSE], 1)),
      strain_2 = NULL
    )
  }

  ## calculate n_doses

  # Output number of first and second doses
  doses <- n_vaccinated[, , c(1, 3), , , drop = FALSE]
  dimnames(doses) <- list(NULL, NULL, c("first_dose", "second_dose"),
                          regions, NULL)
  doses_inc <- aperm(apply(doses, 1:4, diff), c(2:5, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-5], 1)),
                                   doses_inc)
  dimnames(doses_inc) <- dimnames(doses)
  dimnames(doses_inc)[[3L]] <- paste0(dimnames(doses_inc)[[3L]], "_inc")

  n_doses <- abind::abind(doses, doses_inc, along = 3)

  list(n_protected = n_protected,
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
  region <- samples$info$region
  info <- parameters$info[parameters$info$region %in% region, ]
  rownames(info) <- NULL
  for (r in region) {
    i <- which.max(samples$probabilities["log_posterior", r, ])
    initial <- samples$pars[, r, i]
    which <- which(info$region == r)
    info$initial[match(names(initial), info$name[which])] <- unname(initial)
  }

  prior <- parameters$prior[parameters$prior$region %in% region, ]
  rownames(prior) <- NULL

  nr <- nrow(samples$pars)
  covariance <- lapply(
    region, function(x) data.frame(
      region = x, cov(t(samples$pars[, x, ]))))
  covariance <- do.call(rbind, covariance)
  rownames(covariance) <- NULL
  proposal <- data_frame(name = colnames(covariance)[-1],
                          covariance)
  proposal <- proposal[, c("region", "name", colnames(covariance)[-1])]

  list(info = info,
       prior = prior,
       proposal = proposal)
}


calculate_positivity <- function(samples, over25) {

  single_region <- length(samples$info$region) == 1

  model_params <- samples$predict$transform(t(samples$pars[, , 1]))

  if (!single_region) {
    ## FIXME - All params below are currently 'fixed' but what do we do when
    ##  this changes? Ed?
    model_params <- model_params[[1]]
    samples$pars <- t(samples$pars[, 1, ])
  }

  if ("p_NC" %in% colnames(samples$pars)) {
    p_NC <- samples$pars[, "p_NC"]
  } else {
    p_NC <- model_params$p_NC
  }

  if ("p_NC_weekend" %in% colnames(samples$pars)) {
    p_NC_weekend <- samples$pars[, "p_NC_weekend"]
  } else {
    p_NC_weekend <- p_NC
  }

  if (over25) {
    pos <- drop_dimension(
      samples$trajectories$state["sympt_cases_over25_inc", , , , drop = FALSE],
      1)
    neg <- (sum(model_params$N_tot[6:19]) - pos)
  } else {
    pos <- drop_dimension(
      samples$trajectories$state["sympt_cases_inc", , , , drop = FALSE],
      1)
    neg <- (sum(model_params$N_tot) - pos)
  }

  dates <- sircovid::sircovid_date_as_date(samples$trajectories$date)

  neg[, , grepl("^S", weekdays(dates))] <-
    neg[, , grepl("^S", weekdays(dates))] * p_NC_weekend
  neg[, , !grepl("^S", weekdays(dates))] <-
    neg[, , !grepl("^S", weekdays(dates))] * p_NC

  (pos * model_params$pillar2_sensitivity +
    neg * (1 - model_params$pillar2_specificity)) / (pos + neg) * 100
}


inflate_samples <- function(s) {
  region <- s$info$region
  if (length(region) == 1) {
    s$pars <- add_dimension(t(s$pars), 2, region)
    s$probabilities <- add_dimension(t(s$probabilities), 2, region)
    s$state <- add_dimension(s$state, 2, region)
    s$trajectories$state <- add_dimension(s$trajectories$state, 3, region)
    if (!is.null(s$restart)) {
      s$restart$state <- add_dimension(s$restart$state, 3, region)
    }
    s$vaccine <- setNames(list(s$vaccine), region)
  }

  s
}


deflate_samples <- function(s) {
  region <- s$info$region
  if (length(region) == 1) {
    s$pars <- t(drop_dimension(s$pars, 2))
    s$probabilities <- t(drop_dimension(s$probabilities, 2))
    s$state <- drop_dimension(s$state, 2)
    s$trajectories$state <- drop_dimension(s$trajectories$state, 3)
    if (!is.null(s$restart)) {
      s$restart$state <- drop_dimension(s$restart$state, 3)
    }
    s$vaccine <- s$vaccine[[1L]]
  }

  s
}