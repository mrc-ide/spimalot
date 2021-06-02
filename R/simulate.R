spim_simulate_schedule <- function(combined, n_par, end_date, n_threads, keep,
                                   rt_future, rt_type,
                                   future_daily_doses, vaccine_efficacy,
                                   uptake_by_age, n_strain = 1L,
                                   seasonality = 0,
                                   output_time_series = TRUE,
                                   calculate_rt = TRUE,
                                   calculate_vaccination = TRUE,
                                   calculate_state_by_age = TRUE,
                                   vaccine_efficacy_strain_2 = NULL,
                                   strain_transmission = NULL,
                                   strain_cross_immunity = NULL,
                                   strain_seed_date = NULL,
                                   strain_seed_rate = NULL,
                                   inflate_strain = FALSE,
                                   lag_days = NULL, lag_groups = NULL,
                                   rel_strain_modifier = NULL,
                                   prop_voc = NULL,
                                   mean_vacc_delay_multiplier = 1,
                                   seed = NULL) {
  regions <- names(combined$pars)
  transform <- combined$transform[regions]
  vaccine <- combined$vaccine[regions]
  future_daily_doses <- future_daily_doses[regions]
  info <- combined$info[regions]

  strain_seed_rate <- strain_seed_rate[regions]
  strain_seed_rate <- strain_seed_rate %||% list(NULL)

  ## Take a random sample of our parameters without replacement.
  n_regions <- length(regions)
  n_par_combined <- nrow(combined$pars[[1]])
  n_par <- min(n_par, n_par_combined)
  i <- sort(sample(n_par_combined, n_par, replace = FALSE))
  pars_mcmc <- lapply(combined$pars, function(x) x[i, , drop = FALSE])
  state <- lapply(combined$state, function(x) x[, i, drop = FALSE])

  ## First, convert from pmcmc parameters to model parameters; we'll
  ## use and modify these below. This is the only time that the
  ## transform is called.
  ##
  ## TODO: do the transform in parallel? This is quite a timesink
  pars_model <- Map(function(pars_region, transform)
    apply(pars_region, 1, transform),
    pars_mcmc, transform)

  if (inflate_strain) {
    stop("CHECK")
    ## At this point we should have multistrain = FALSE and n_strain > 1
    tmp <- prepare_inflate_strain(pars_model, state, info)
    pars_model <- tmp$pars
    state <- tmp$state
    info <- tmp$info
  }

  priority_population <- lapply(names(pars_model),
                                sircovid::vaccine_priority_population,
                                uptake = uptake_by_age)

  ## Seems odd that this is not doing the future betas too
  pars_future <- Map(
    simulate_prepare_parameters_future,
    pars = pars_model,
    vaccine = vaccine,
    future_daily_doses = future_daily_doses,
    priority_population = priority_population,
    strain_seed_rate = strain_seed_rate,
    MoreArgs = list(
      vaccine_efficacy = vaccine_efficacy,
      vaccine_efficacy_strain_2 = vaccine_efficacy_strain_2,
      n_strain = n_strain,
      end_date = end_date,
      strain_seed_date = strain_seed_date,
      strain_cross_immunity = strain_cross_immunity,
      strain_transmission = strain_transmission,
      lag_days = lag_days,
      lag_groups = lag_groups,
      rel_strain_modifier = rel_strain_modifier,
      mean_vacc_delay_multiplier = mean_vacc_delay_multiplier
    )
  )

  ## For the final object we will use a list-matrix of parameters and
  ## a 3d array of state as these will feed more easily into dust.
  pars_future <- unlist(unname(pars_future), FALSE)
  dim(pars_future) <- c(n_par, n_regions)
  colnames(pars_future) <- regions

  state <- array(unlist(unname(state)),
                 c(nrow(state[[1]]), n_par, n_regions))

  steps_per_day <- combined$steps_per_day

  date_start <- combined$step / steps_per_day
  end_date <- sircovid::sircovid_date(end_date)
  dates <- seq(date_start, end_date)
  steps <- dates * steps_per_day
  step_end <- end_date * steps_per_day

  info <- combined$info[[1]]$info
  index <- simulate_index(info, keep, calculate_vaccination, n_strain > 1)

  S <- mcstate::array_flatten(state[index$S, , , drop = FALSE], 2:3)

  if (n_strain > 1) {
    R <- mcstate::array_flatten(state[index$R, , , drop = FALSE], 2:3)
    prob_strain <- mcstate::array_flatten(
      state[index$prob_strain, , , drop = FALSE], 2:3)
  } else {
    R <- NULL
    prob_strain <- NULL
  }

  pars <- simulate_prepare_future_betas(
    pars_future, rt_future, S, rt_type,
    combined$step, step_end, 1 / steps_per_day,
    seasonality, R, prob_strain)

  if (!is.null(prop_voc)) {
    ## NOTE: not 100% sure that this is correct, was 'dat$info' and
    ## I've renamed to combined$info
    state <- move_strain_compartments(
      state, info, c("E", "I_A", "I_P", "I_C_1"),
      1, 2, prop_voc, names(combined$info)
    )
  }

  obj <- sircovid::carehomes$new(pars, combined$step, NULL,
                                 n_threads = n_threads,
                                 seed = seed,
                                 pars_multi = TRUE)
  obj$set_state(state)
  obj$set_index(index$run)

  trajectories <- obj$simulate(steps)
  dimnames(trajectories)[[3]] <- names(combined$pars)

  ret <- list(date = dates,
              summary_state = simulate_summary_state(trajectories, keep, dates))

  if (output_time_series) {
    ret$state <- trajectories[keep, , , ]
  }

  if (calculate_state_by_age) {
    ret$state_by_age <- simulate_summary_by_age(trajectories, index)
  }

  if (calculate_rt) {
    critical_dates <- unique(sircovid::sircovid_date(rt_future$date))
    critical_dates <- critical_dates[critical_dates > date_start]
    message("Calculating Rt")
    rt <- simulate_rt(steps,
                      trajectories[names(index$S), , , ],
                      pars,
                      sort(critical_dates),
                      trajectories[names(index$R), , , ],
                      trajectories[names(index$prob_strain), , , ],
                      no_seeding = identical(strain_seed_rate[[1]], numeric(2)),
                      prop_voc = prop_voc)
    ret <- c(ret, rt)
  }

  if (calculate_vaccination) {
    ret <- c(ret,
             simulate_summary_vaccination(trajectories, index, vaccine_efficacy,
                                          n_strain))
  }

  ret
}


simulate_prepare <- function(combined, n_par, end_date,
                             ## vaccination
                             future_daily_doses, vaccine_efficacy,
                             uptake_by_age, lag_days = NULL, lag_groups = NULL,
                             mean_vacc_delay_multiplier = 1,
                             ## strain
                             n_strain = 1L,
                             strain_transmission = NULL,
                             strain_seed_date = NULL,
                             strain_seed_rate = NULL,
                             strain_cross_immunity = NULL,
                             inflate_strain = FALSE,
                             vaccine_efficacy_strain_2 = NULL,
                             rel_strain_modifier = NULL) {

  ## Extra things that we might want to change soon
  ret$no_seeding <- identical(strain_seed_rate[[1]], numeric(2))
  ret$vaccine_efficacy <- vaccine_efficacy
  ret$multistrain <- n_strain > 1
  ret
}


simulate_prepare_parameters_future <- function(pars, vaccine,
                                               future_daily_doses, end_date,
                                               vaccine_efficacy,
                                               priority_population,
                                               n_strain = NULL,
                                               vaccine_efficacy_strain_2 = NULL,
                                               strain_seed_rate = NULL,
                                               strain_seed_date = NULL,
                                               strain_cross_immunity = NULL,
                                               strain_transmission = NULL,
                                               lag_days = NULL,
                                               lag_groups = NULL,
                                               rel_strain_modifier = NULL,
                                               mean_vacc_delay_multiplier = 1) {

  if (is.null(vaccine_efficacy)) {
    stop("FIXME: Rewriteto allow not changing efficacy")
  }

  p <- pars ## just a rename for historical reasons, could use pars throughout
  vaccine_index_dose2 <- p[[1]]$index_dose[[2]]
  vaccine_progression_rate <- p[[1]]$vaccine_progression_rate_base

  N_tot <- p[[1]]$N_tot

  vaccine_schedule <- sircovid::vaccine_schedule_scenario(
    vaccine$schedule, future_daily_doses, end_date,
    round(vaccine$mean_days_between_doses * mean_vacc_delay_multiplier),
    priority_population, lag_groups, lag_days)

  n_vacc_strata <- ncol(vaccine_efficacy[[1]])
  n_groups <- 19
  if (!all(lengths(vaccine_efficacy) == n_vacc_strata * n_groups)) {
    stop("Vaccine efficacy parameters must have the same length")
  }
  if (!is.null(vaccine_efficacy_strain_2) &&
      !all(lengths(vaccine_efficacy_strain_2) == n_vacc_strata * n_groups)) {
    stop("Vaccine efficacy strain 2 parameters must have the same length")
  }


  ## vaccine_efficacy parameters come in as vectors of length n_vacc_strata
  ## let's create arrays of dim = c(n_groups, n_strain, n_vacc_strata)
  ## these are required for sircovid:::carehomes_parameters_vaccination
  dim <- c(n_groups, n_strain, n_vacc_strata)
  rel_list <- rep(list(array(rep(NA_integer_), dim = dim)), 4)
  names(rel_list) <- c("rel_p_sympt", "rel_p_hosp_if_sympt",
                       "rel_susceptibility", "rel_infectivity")

  if (!is.null(vaccine_efficacy_strain_2) && n_strain != 4) {
    stop("'n_strain' should be '4' iff 'vaccine_efficacy_strain_2' is non-NULL")
  } else if (is.null(vaccine_efficacy_strain_2) && n_strain != 1) {
    stop("'n_strain' should be '1' iff 'vaccine_efficacy_strain_2' is NULL")
  }

  for (rel in names(rel_list)) {
    for (s in seq_len(n_strain)) {
      if (is.null(rel_strain_modifier)) {
        mod <- 1
      } else {
        mod <- rel_strain_modifier[[s]][[rel]]
      }
      for (v_s in seq_len(n_vacc_strata)) {
        for (g in seq_len(n_groups)) {
          ## If not multistrain then all use same params, otherwise split
          ##  by strains. Strains: 1 (=1), 2(=2), 3(=1->2), 4(=2->1)
          if (is.null(vaccine_efficacy_strain_2) || s %in% c(1, 4)) {
            ## Strains 1 and 2->1 (if multistrain)
            rel_list[[rel]][g, s, v_s] <- vaccine_efficacy[[rel]][g, v_s] * mod
          } else {
            ## Strains 2 and 1->2 (if multistrain)
            rel_list[[rel]][g, s, v_s] <-
              vaccine_efficacy_strain_2[[rel]][g, v_s] * mod
          }
        }
      }
    }
  }

  extra <- sircovid:::carehomes_parameters_vaccination(
    N_tot,
    p[[1]]$dt,
    rel_susceptibility = rel_list$rel_susceptibility,
    rel_p_sympt = rel_list$rel_p_sympt,
    rel_p_hosp_if_sympt = rel_list$rel_p_hosp_if_sympt,
    rel_infectivity = rel_list$rel_infectivity,
    vaccine_schedule = vaccine_schedule,
    vaccine_index_dose2 = vaccine_index_dose2,
    vaccine_progression_rate = vaccine_progression_rate,
    n_strains = n_strain
  )

  if (!is.null(strain_transmission)) {
    strain_params <- sircovid:::carehomes_parameters_strain(
      strain_transmission, sircovid::sircovid_date(strain_seed_date),
      strain_seed_rate, dt = p[[1]]$dt
    )
    strain_params$cross_immunity <- strain_cross_immunity
    extra <- c(extra, strain_params)
  }

  for (i in seq_along(p)) {
    p[[i]][names(extra)] <- extra
  }

  lapply(p, sircovid::carehomes_check_severity)
}


simulate_prepare_future_betas <- function(pars, rt_future, S, rt_type,
                                          step_current, step_end, dt,
                                          seasonality, R, prob_strain) {
  ## For seasonality, we assume a sine wave with a trough at day 228 = 15th Aug
  ## (and a peak 6 months earlier on day 46 = 15th Feb)
  seasonality_date_peak <- sircovid::sircovid_date("2020-02-15")

  ## Past betas, as inferred from the pmcmc
  beta <- t(vapply(pars, "[[", numeric(length(pars[[1]]$beta_step)),
                   "beta_step"))

  n_beta_add <- step_end - ncol(beta)
  beta <- cbind(beta, matrix(rep(beta[, ncol(beta)], n_beta_add), nrow(beta)))
  beta <- array(beta, c(dim(pars), ncol(beta)))

  rt <- vnapply(seq_along(pars), function(i)
    sircovid::carehomes_Rt(step_current, S[, i, drop = FALSE], pars[[i]],
                           type = rt_type, R = R[, i, drop = FALSE],
                           prob_strain = prob_strain[, i, drop = FALSE],
                           weight_Rt = FALSE
    )[[rt_type]][1])
  beta_rt_ratio <- beta[, , step_current] / rt

  for (region_index in seq_len(ncol(pars))) {
    r <- colnames(pars)[[region_index]]
    rt_future_r <- rt_future[rt_future$region == r, ]
    rt_future_r$step_start <- sircovid::sircovid_date(rt_future_r$date) / dt
    rt_future_r$step_end <- c(rt_future_r$step_start[-1L] - 1L, step_end)

    for (i in seq_len(nrow(rt_future_r))) {
      j <- seq(rt_future_r$step_start[[i]], rt_future_r$step_end[[i]])

      if (rt_future_r$Rt_sd[[i]] > 0) {
        dpars <- data.frame(mean = rt_future_r$Rt[[i]],
                            sd = rt_future_r$Rt_sd[[i]])

        # order by increasing Rt to preserve direction of Rt scaling
        rt_i <- sort(distr6::dstr("Lognormal",
                                  mean = rt_future_r$Rt[[i]],
                                  sd = rt_future_r$Rt_sd[[i]])$rand(nrow(beta)))
      } else {
        rt_i <- rt_future_r$Rt[[i]]
      }

      beta[, region_index, j] <- rt_i * beta_rt_ratio[, region_index]
    }
  }

  beta <- mcstate::array_flatten(beta, 1:2)

  ## Apply seasonality to the new betas:

  ## Back-calculate dates here to deal with any truncating of
  ## step_current:step_end caused by simulating with new data.
  date <- seq(to = step_end, length.out = n_beta_add) * dt
  ## Relative distance between us and the peak date (on [0..1])
  beta_mult <- calc_seasonality(date, seasonality_date_peak, seasonality)

  i <- seq(to = ncol(beta), length.out = n_beta_add)
  beta[, i] <- beta[, i] * rep(beta_mult, each = nrow(beta))

  for (i in seq_along(pars)) {
    pars[[i]]$beta_step <- beta[i, ]
  }

  pars
}


calc_seasonality <- function(date, seasonality_date_peak, seasonality) {
  delta <- ((date - seasonality_date_peak) %% 365) / 365
  1 + cos(2 * pi * delta) * seasonality
}


simulate_model_init <- function(pars, step, state, n_threads, seed = NULL) {
  obj <- sircovid::carehomes$new(pars, step, NULL, n_threads = n_threads,
                                 seed = seed, pars_multi = TRUE)
  obj$set_state(state)
  obj
}


simulate_summary_by_age <- function(state, index) {
  n_groups <- sircovid:::carehomes_n_groups()

  ## output cumulative states by
  ## age / vaccine class / sample / region / time
  arrays <- list(
    deaths <- state[names(index$D), , , ],
    infections <- state[names(index$I), , , ],
    admissions <- state[names(index$A), , , ]
  )
  names(arrays) <- c("deaths", "infections", "diagnoses_admitted")
  strata <- nrow(arrays$deaths) / n_groups

  f <- function(array) {

    x <- mcstate::array_reshape(array, 1L, c(n_groups, strata))

    ## aggregate partially immunised strata
    x[, 2L, , , ] <- x[, 2L, , , ] + x[, 3L, , , ]
    x <- x[, -3L, , , ]
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
    ret <- apply(abind::abind(res, along = 5), c(1, 3, 4, 5), mean)

    # [age, vaccine status, region, time]
    ret <- round(aperm(ret, c(4, 1, 2, 3)))

    ret
  }

  lapply(arrays, f)
}


simulate_index <- function(info, keep, calculate_vaccination, multistrain) {
  index_S <- info$index$S
  names(index_S) <- paste0("S_", seq_along(index_S))

  index_D <- info$index$D
  names(index_D) <- paste0("D_", seq_along(index_D))

  index_I <- info$index$cum_infections_disag
  names(index_I) <- paste0("I_", seq_along(index_I))

  index_A <- info$index$diagnoses_admitted
  names(index_A) <- paste0("A_", seq_along(index_A))

  if (calculate_vaccination || multistrain) {
    index_n_vaccinated <- set_names(
      info$index$cum_n_vaccinated,
      paste0("n_vaccinated_", seq_along(index_S))
    )
    index_R <- info$index$R
    names(index_R) <- paste0("R_", seq_along(index_R))
  } else {
    index_n_vaccinated <- NULL
    index_R <- NULL
  }

  if (multistrain) {
    index_prob_strain <- info$index$prob_strain
    names(index_prob_strain) <-
      paste0("prob_strain_", seq_along(index_prob_strain))
  } else {
    index_prob_strain <- NULL
  }

  index <- c(sircovid::carehomes_index(info)$state[keep],
             index_S, index_D, index_I, index_A, index_n_vaccinated,
             index_R, index_prob_strain)

  list(run = index,
       n_vaccinated = index_n_vaccinated,
       S = index_S,
       I = index_I,
       A = index_A,
       R = index_R,
       D = index_D,
       prob_strain = index_prob_strain)
}


simulate_summary_state <- function(state, keep, dates) {
  summary_state <- state[keep, , , ]
  # aggregate over regions to get England peak hosp
  summary_state_england <-
    apply(summary_state[, , sircovid::regions("england"), ], c(1, 2, 4), sum)
  summary_state <- abind_quiet(list(summary_state,
                                    england = summary_state_england),
                               along = 3)

  # calculate number and date of peak hospital bed occupancy
  peak_hosp <- list(
    peak_hosp = apply(summary_state["hosp", , , ], c(1, 2), max),
    peak_hosp_date = apply(summary_state["hosp", , , ], c(1, 2), which.max))
  peak_hosp$peak_hosp_date[] <- dates[peak_hosp$peak_hosp_date]
  peak_hosp <- aperm(abind_quiet(peak_hosp, along = 3), c(3, 1, 2))


  ## reset trajectories to zero
  summary_state <- summary_state[setdiff(keep, c("hosp", "icu")), , , ]
  summary_state <- summary_state[, , , dim(state)[4]] - summary_state[, , , 1]

  # add diagnoses admitted
  diagnoses_admitted <- summary_state["diagnoses", , , drop = FALSE] +
    summary_state["admitted", , , drop = FALSE]
  rownames(diagnoses_admitted) <- "diagnoses_admitted"

  # bind results together
  summary_state <- abind_quiet(summary_state, diagnoses_admitted, peak_hosp,
                               along = 1)

  abind_quiet(summary_state, along = 4)
}


simulate_summary_vaccination <- function(state, index, vaccine_efficacy,
                                         n_strain) {
  n_groups <- sircovid:::carehomes_n_groups()
  regions <- dimnames(state)[[3]]

  ## output the cumulative transitions between vaccine strata
  ## by age / vaccine stratum / region / over time
  n_vaccinated <- apply(state[names(index$n_vaccinated), , , , drop = FALSE],
                        c(1, 3, 4), mean)
  n_strata <- nrow(n_vaccinated) / n_groups
  n_vaccinated <-
    mcstate::array_reshape(n_vaccinated, 1L, c(n_groups, n_strata))

  ## output the number recovered in each vaccine stratum / region / over time
  R_raw <- mcstate::array_reshape(
    state[names(index$R), , , , drop = FALSE], 1L,
    c(n_groups, n_strain, n_strata))
  ## the mean is taken here over age and particle but we want the
  ## number for each age group so you multiply back by n_groups = 19
  ## only extract strain 1
  R_strain_1 <- R_raw[, 1, , , , ]
  R <- apply(R_strain_1, c(2, 4, 5), mean) * n_groups

  # calculate the proportion protected given vaccine efficacy
  n_protected <- simulate_summary_protected(n_vaccinated, R, vaccine_efficacy)
  dimnames(n_protected)[[2]] <- regions

  # Output number of first and second doses
  doses <- n_vaccinated[, c(1, 3), , ]
  dimnames(doses)[2:3] <- list(c("first_dose", "second_dose"), regions)
  doses_inc <- aperm(apply(doses, c(1, 2, 3), diff), c(2, 3, 4, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-4], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")

  n_doses <- abind_quiet(doses, doses_inc, along = 2)

  list(n_vaccinated = n_vaccinated,
       n_protected = n_protected,
       n_doses = n_doses)
}


## multistrain controls if VOC is seeded, if not turn off to remove
##  NA prob_strain errors
simulate_rt <- function(steps, S, pars, critical_dates, R = NULL,
                        prob_strain = NULL, no_seeding = FALSE,
                        prop_voc = NULL) {
  dim_S <- dim(S)
  region_names <- dimnames(S)[[3]]
  ## TODO: add mcstate::array_combine(S, 2:3)
  dim(S) <- c(dim_S[[1]], dim_S[[2]] * dim_S[[3]], dim_S[[4]])

  if (!is.null(R)) {
    dim_R <- dim(R)
    dim(R) <- c(dim_R[[1]], dim_R[[2]] * dim_R[[3]], dim_R[[4]])
  }

  if (!is.null(prob_strain)) {
    if (no_seeding && is.null(prop_voc)) {
      ## if not seeding second VOC then just manually set prob_strain
      ## to avoid NAs
      prob_strain[1, , , ] <- 1
      prob_strain[2, , , ] <- 0
    }
    dim_prob_strain <- dim(prob_strain)
    dim(prob_strain) <- c(dim_prob_strain[[1]],
                          dim_prob_strain[[2]] * dim_prob_strain[[3]],
                          dim_prob_strain[[4]])
  }

  rt_type <- c("Rt_general", "eff_Rt_general")

  rt <- sircovid::carehomes_Rt_trajectories(
    steps, S, pars, type = rt_type,
    initial_step_from_parameters = FALSE,
    interpolate_critical_dates = critical_dates,
    interpolate_every = 7,
    interpolate_min = 3, R = R, prob_strain = prob_strain, weight_Rt = TRUE)

  ## Better format for the Rt values, will move into sircovid.
  ret_rt <- rt[rt_type]
  for (i in rt_type) {
    ret_rt[[i]] <- mcstate::array_reshape(t(rt[[i]]), 1, dim_S[2:3])
    dimnames(ret_rt[[i]]) <- list(NULL, region_names, NULL)
  }

  ret_rt
}


## TODO: almost the same as fit_process.R:calculate_vaccination and
## needs merge; differences are input args slightly different and
## output format appears different
simulate_summary_protected <- function(n_vaccinated, R, vaccine_efficacy) {
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
  #
  # - number vaccinated sum_sr(V[s, r, t])
  # - number protected after vaccination against severe disease:
  #   i.e., sum_asr(V[a, s, r, t] * eff_severe[a, s])
  # - number protected after vaccination against infection
  #   i.e., sum_asr(V[a, s, r, t] * eff_inf[a, s])
  # - number protected after infection (this is added to all of the above):
  #   i.e., sum_sr(R[s, r, t])
  # - number protected after infection only:
  #   i.e., sum_r(R[1, r, t]


  ## create array of number in each vaccine stratum
  n_strata <- ncol(n_vaccinated)
  V <- array(0, dim = dim(n_vaccinated))
  V[, n_strata, , ] <- n_vaccinated[, n_strata - 1L, , ]
  for (i in seq(2, n_strata - 1)) {
    V[, i, , ] <- n_vaccinated[, i - 1, , ] - n_vaccinated[, i, , ]
  }

  sum_sr <- function(x) apply(x, c(2, 3), sum)
  sum_asr <- function(x) apply(x, c(3, 4), sum)

  ret <- list(
    ever_vaccinated = sum_sr(n_vaccinated[, 1, , ]),
    protected_against_infection = sum_asr(c(vp$infection) * V),
    protected_against_severe_disease = sum_asr(c(vp$severe_disease) * V),
    ever_infected = sum_sr(R),
    ever_infected_unvaccinated = R[1, , ]
  )

  aperm(abind_quiet(ret, along = 3), c(3, 1, 2))
}


prepare_inflate_strain <- function(pars, state, info) {
  ## First inflate the parameters, because we need these in order to
  ## inflate the state.

  ## Easy updates:
  update <- list(cross_immunity = c(1, 1),
                 n_strains = 4,
                 strain_transmission = rep(1, 4))

  inflate_pars <- function(p_i) {
    p_i[names(update)] <- update
    for (rel in grep("^rel_", names(p_i), value = TRUE)) {
      rel_old <- p_i[[rel]]
      ## rel_gamma_X
      if (is.null(dim(rel_old))) {
        p_i[[rel]] <- rep(rel_old, 4)
      } else {
        rel_new <- array(0, c(nrow(rel_old), 4, nlayers(rel_old)))
        rel_new[, , ] <- rel_old[, 1, , drop = FALSE]
        p_i[[rel]] <- rel_new
      }
    }
    p_i
  }

  pars_new <- lapply(pars, lapply, inflate_pars)

  ## Then update the states given that:
  info_old <- info[[1]]$info
  info_new <- sircovid::carehomes$new(pars_new[[1]][[1]], 0, 1)$info()
  state_new <- lapply(state, sircovid::inflate_state_strains,
                      info_old, info_new)

  info_new <- lapply(info, function(x) {
    x$info <- info_new
    x
  })

  list(pars = pars_new, state = state_new, info = info_new)
}


move_strain_compartments <- function(state, info, compartment,
                                     original_strain_idx, new_strain_idx,
                                     prop, regions) {
  for (i in seq_along(compartment)) {
    dim <- info$dim[[compartment[[i]]]]
    for (r in seq_along(prop)) {
      new_state <- state[info$index[[compartment[[i]]]], , r]
      tmp_dim <- dim(new_state)
      new_state <- array(new_state, c(dim, ncol(new_state)))

      if (length(dim) == 4) {
        new_state[, 2, , , ] <-
          round(prop[[regions[[r]]]] * new_state[, 1, , , ])
        new_state[, 1, , , ] <- new_state[, 1, , , ] - new_state[, 2, , , ]
      } else {
        stop(sprintf("Unexpected dimensions (%d) in move_strain_compartment",
                     length(dim)))
      }
      state[info$index[[compartment[[i]]]], , r] <- array(new_state, tmp_dim)
    }
  }
  state
}
