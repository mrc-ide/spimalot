##' Prepare combined output for simulations. This carries out
##' maintenance for the state - upgrading state if sircovid has
##' updated (within reason), sampling starting parameters, reducing to
##' selected regions and (potentially) adding empty strain or booster
##' compartments.
##'
##' @title Prepare for simulation
##' @param combined A "combined" object
##'
##' @param simulate_parameters A nested list of simulation parameters
##'
##' @param n_par The number of parameters to sample
##'
##' @param regions Character vector of regions to use
##'
##' @param seed_voc Logical, indicating if seeding a new VOC
##'
##' @export
spim_simulate_prepare <- function(combined, simulate_parameters, n_par,
                                  regions = NULL,
                                  seed_voc = FALSE) {
  if (is.null(regions)) {
    regions <- sircovid::regions("all")
  }

  combined <- simulate_prepare_drop_regions(combined, regions)
  combined <- simulate_prepare_upgrade(combined)

  if (seed_voc) {
    message("Implementing date to seed new VOC")
    combined <- simulate_seed_parameters(combined, regions, simulate_parameters)
  }

  info <- combined$info

  ## Take a random sample of our parameters without replacement.
  n_regions <- length(regions)
  n_par_combined <- nrow(combined$pars[[1]])
  if (n_par > n_par_combined) {
    message(sprintf(
      "Reducing n_par from %d to %d as too few available in combined",
      n_par, n_par_combined))
    n_par <- n_par_combined
  }
  i <- sort(sample(n_par_combined, n_par, replace = FALSE))
  pars_mcmc <- lapply(combined$pars, function(x) x[i, , drop = FALSE])
  state <- lapply(combined$state, function(x) x[, i, drop = FALSE])

  ## TODO: this is quite slow
  message("Creating odin parameters from sampled parameters")
  pars <- Map(function(pars_region, transform)
    apply(pars_region, 1, transform),
    pars_mcmc, combined$transform)

  ## For the final object we will use a list-matrix of parameters and
  ## a 3d array of state as these will feed more easily into dust.
  pars <- unlist(unname(pars), FALSE)
  dim(pars) <- c(n_par, n_regions)
  colnames(pars) <- regions

  state <- array(unlist(unname(state)),
                 c(nrow(state[[1]]), n_par, n_regions))

  ## Our final object that we will use in the simulations
  ret <- combined[c("step", "date", "dt", "steps_per_day", "base")]
  ret$pars <- pars
  ret$state <- state
  ret$info <- info
  ret
}


##' Create a list of simulation parameters
##'
##' @title Create simulation parameters
##'
##' @param grid A scenario grid; this is a data.frame with names being
##'   the type of thing being changed and the values representing
##'   quantities in vars
##'
##' @param vars A list of parameters organised by scenario
##'
##' @param base A list of base parameters
##'
##' @param ignore Character vector of additional values in `grid` that
##'   are not directly model parameters
##'
##' @param regions Character vector of regions
##'
##' @param multistrain Logical, indicating if the simulation will use
##'   multiple strains
##'
##' @return A list of parameters for the model
##' @export
spim_simulate_args <- function(grid, vars, base, ignore, regions, multistrain) {

  simulate_args_validate(grid, vars, base, ignore, multistrain)

  f <- function(i) {
    el <- grid[i, ]
    for (nm in setdiff(names(el), ignore)) {
      level <- el[[nm]]
      if (!(level %in% names(vars[[nm]]))) {
        stop(sprintf("'%s' not found in vars$%s", level, nm))
      }
      ## Special treatment in case the value really is NULL; make sure
      ## we get a named NULL in the result rather than deleting the
      ## entry
      value <- vars[[nm]][[level]]
      if (is.null(value)) {
        base[nm] <- list(NULL)
      } else {
        base[[nm]] <- value
      }
    }
    base[ignore] <- el[ignore]
    base
  }

  ret <- lapply(seq_rows(grid), f)

  message("Validating generated parameters")
  for (i in seq_along(ret)) {
    tryCatch(
      simulate_validate_args1(ret[[i]], regions, multistrain),
      error = function(e)
        stop(sprintf("While checking args[[%d]]: %s", i, e$message)))
  }

  ret
}

##' Run simulations locally
##'
##' @title Run simulations locally
##'
##' @param args Arguments returned by [spim_simulate_args]
##'
##' @param combined Processed combined output returned by
##'   [spim_simulate_prepare]
##'
##' @export
spim_simulate_local <- function(args, combined) {
  lapply(seq_along(args), spim_simulate, args, combined)
}

##' Run simulations with rrq workers
##'
##' @title Run simulations with rrq
##'
##' @param args Arguments returned by [spim_simulate_args]
##'
##' @param combined Processed combined output returned by
##'   [spim_simulate_prepare]
##'
##' @param rrq rrq object
##'
##' @export
spim_simulate_rrq <- function(args, combined, rrq) {
  rrq$lapply(seq_along(args), spim_simulate, args, combined)
}


simulate_args_names <- function(multistrain = TRUE) {
  args <-
    c(## Core simulation parameters
      "end_date", "seed", "n_threads",
      ## Output control
      "output_keep", "output_rt", "output_time_series", "output_vaccination",
      "output_state_by_age", "output_weight_rt",
      ## Rt control
      "rt_type",
      "rt_future",
      ## Seasonality
      "seasonality",
      ## Vaccination
      "vaccine_daily_doses", "vaccine_booster_daily_doses",
      "vaccine_efficacy", "vaccine_booster_efficacy", "vaccine_eligibility",
      "vaccine_uptake", "vaccine_lag_groups", "vaccine_lag_days",
      "vaccine_delay_multiplier",
      ## waning
      "waning_rate"
    )

  if (multistrain) {
    args <-
      c(args,
        "strain_seed_date", "strain_transmission", "strain_seed_rate",
        "strain_vaccine_efficacy", "strain_initial_proportion",
        "strain_vaccine_booster_efficacy", "strain_cross_immunity",
        "strain_severity_modifier")
  }

  args
}


spim_simulate_one <- function(args, combined, move_between_strains = FALSE) {
  ## TODO: run validate here again, requires moving ignore into the
  ## object though.
  multistrain <- combined$info[[1]]$multistrain
  if (multistrain) {
    n_strain <- 4
  } else {
    n_strain <- 1
  }

  regions <- names(combined$info)

  ## Lots of updates to parameters to account for how vaccination
  ## changes over the future.

  step_start <- combined$step
  date_start <- sircovid::sircovid_date(combined$date)
  end_date <- sircovid::sircovid_date(args$end_date)
  dates <- seq(date_start, end_date)
  steps <- dates * combined$steps_per_day
  step_end <- last(steps)

  info <- combined$info[[1]]$info
  index <- simulate_index(info, args$output_keep,
                          args$output_vaccination,
                          multistrain)

  state_start <- combined$state

  S <- mcstate::array_flatten(state_start[index$S, , , drop = FALSE], 2:3)

  if (multistrain) {
    R <- mcstate::array_flatten(state_start[index$R, , , drop = FALSE], 2:3)
    prob_strain <- mcstate::array_flatten(
      state_start[index$prob_strain, , , drop = FALSE], 2:3)
  } else {
    R <- NULL
    prob_strain <- NULL
  }

  pars <- lapply(regions, simulate_one_pars_vaccination, args, combined,
                 n_strain)
  pars <- unlist(pars, FALSE, FALSE)
  attributes(pars) <- attributes(combined$pars)

  ## TODO: if we could reuse the rt that we had it and avoid quite a
  ## bit of time here.
  pars <- setup_future_betas(pars, args$rt_future, S, args$rt_type, step_start,
                             step_end, combined$dt, args$seasonality, R,
                             prob_strain)

  if (!is.null(args$strain_initial_proportion) && move_between_strains) {
    state_start <- move_strain_compartments(
      state_start, info, c("E", "I_A", "I_P", "I_C_1"),
      1, 2, args$strain_initial_proportion, regions)
  }

  for (i in seq_along(pars)) {
    pars[[i]]$steps_per_day <- as.integer(pars[[i]]$steps_per_day)
    pars[[i]]$index_dose <- as.integer(pars[[i]]$index_dose)
  }
  message("Creating dust object")
  obj <- sircovid::lancelot$new(pars, step_start, NULL, pars_multi = TRUE,
                                n_threads = args$n_threads, seed = args$seed)
  obj$update_state(state = state_start)
  obj$set_index(index$run)
  message("Simulating!")
  state <- obj$simulate(steps)
  dimnames(state)[[3]] <- regions

  message("Adding summary statistics")
  ret <- list(
    date = dates,
    summary_state = create_summary_state(state, args$output_keep, dates))

  if (args$output_time_series) {
    ret$state <- state[args$output_keep, , , ]
  }

  if (args$output_state_by_age) {

    ret$state_by_age <- simulate_extract_age_class_state(state, index)
  }

  if (args$output_rt) {
    critical_dates <- unique(sircovid::sircovid_date(args$rt_future$date))
    critical_dates <- critical_dates[critical_dates > date_start]
    message("Calculating Rt")
    rt <- simulate_rt(
      steps,
      state[names(index$S), , , ],
      pars,
      sort(critical_dates),
      args$voc_seeded,
      state[names(index$R), , , ],
      state[names(index$prob_strain), , , ],
      no_seeding = identical(args$strain_seed_rate[[1]], numeric(1)),
      prop_voc = args$strain_initial_proportion,
      weight_Rt = FALSE)
    if (args$output_weight_rt) {
      weighted_rt <- simulate_rt(
        steps,
        state[names(index$S), , , ],
        pars,
        sort(critical_dates),
        args$voc_seeded,
        state[names(index$R), , , ],
        state[names(index$prob_strain), , , ],
        no_seeding = identical(args$strain_seed_rate[[1]], numeric(1)),
        prop_voc = args$strain_initial_proportion,
        weight_Rt = TRUE)

      # combined weighted and strain specific outputs
      rt <- Map(function(rt, weighted_rt) {
        x <- abind_quiet(rt, weighted_rt, along = 4)
        dimnames(x)[[4]] <- c("strain_1", "strain_2", "both")
        x}, rt, weighted_rt)
    }


    ret <- c(ret, rt)
  }

  if (args$output_vaccination) {

    if (is.null(args$strain_vaccine_booster_efficacy)) {

      ret <-
        c(ret,
          simulate_calculate_vaccination(state, index,
                                         args$vaccine_efficacy,
                                         args$vaccine_booster_efficacy,
                                         n_strain,
                                         args$strain_vaccine_efficacy,
                                         args$strain_vaccine_booster_efficacy,
                                         args$strain_cross_immunity))

    } else {

      rel_list <- pars[[1]][names(args$vaccine_efficacy)]
      vaccine_efficacy_strain_1 <- lapply(rel_list, "[", , 1, -5)
      vaccine_efficacy_strain_2 <- lapply(rel_list, "[", , 2, -5)
      booster_efficacy_strain_1 <- lapply(rel_list, "[", , 1, 5)
      booster_efficacy_strain_2 <- lapply(rel_list, "[", , 2, 5)

      ret <-
        c(ret,
          simulate_calculate_vaccination(
            state, index,
            vaccine_efficacy = vaccine_efficacy_strain_1,
            booster_efficacy = booster_efficacy_strain_1,
            n_strain,
            strain_vaccine_efficacy = vaccine_efficacy_strain_2,
            strain_vaccine_booster_efficacy = booster_efficacy_strain_2,
            args$strain_cross_immunity))
    }
  }

  ret
}


### prep
simulate_prepare_upgrade <- function(combined) {
  ours <- packageVersion("sircovid")
  for (i in seq_along(combined$info)) {
    if (combined$info[[i]]$version != ours) {
      message(sprintf("Upgrading state for '%s' (%s => %s)",
                      names(combined$state)[[i]],
                      combined$info[[i]]$version, ours))
      ## NOTE: this won't work well if we have to add new parameters
      ## because we then break the transform function.
      p <- combined$transform[[i]](combined$pars[[i]][1, ])
      p_last_epoch <- p[[length(p)]]$pars
      info_new <- sircovid::lancelot$new(p_last_epoch, 0, 1)$info()
      cmp <- combined$state[[i]]
      combined$state[[i]] <- sircovid::upgrade_state(
        combined$state[[i]],
        combined$info[[i]]$info[[length(p)]],
        info_new)
      combined$info[[i]]$info <- info_new
    }
  }

  combined
}


simulate_seed_parameters <- function(combined, regions, simulate_parameters) {

  stopifnot("strain_seed_date" %in% names(simulate_parameters))
  strain_seed_date <- simulate_parameters$strain_seed_date
  for (i in regions) {
    combined$pars[[i]][, "strain_seed_date"] <- strain_seed_date
  }
  combined

}



simulate_prepare_drop_regions <- function(combined, regions) {
  nms <- setdiff(names(combined),
                 c("date", "step", "dt", "steps_per_day", "simulate"))
  for (i in nms) {
    msg <- setdiff(regions, names(combined[[i]]))
    if (length(msg) > 0) {
      stop(sprintf("Missing regions from %s: %s", i,
                   paste(squote(msg), collapse = ", ")))
    }
    combined[[i]] <- combined[[i]][regions]
  }
  combined
}


simulate_prepare_inflate_strain <- function(pars, state, info) {
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
  info_new <- sircovid::lancelot$new(pars_new[[1]][[1]], 0, 1)$info()
  state_new <- lapply(state, sircovid::inflate_state_strains,
                      info_old, info_new)

  info_new <- lapply(info, function(x) {
    x$info <- info_new
    x$multistrain <- TRUE
    x
  })

  list(pars = pars_new, state = state_new, info = info_new)
}


## Sometimes if a run is done with no booster we still want
## to use it with boosters; this will update the parameters and the
## state in order to allow this. This takes the (almost) final output
## of simulate_prepare_
simulate_prepare_inflate_vacc_classes <- function(pars, state, info) {
  ## First inflate the parameters, because we need these in order to
  ## inflate the state.

  ## How many strata are we adding
  n_new_vacc_classes <- 1L # hardcoded for now

  ## of the new strata, are any corresponding to a dose (with a schedule)
  ## and where are they in the additional strata.
  ## hardcoded for now # first of the additional vaccine classes
  idx_new_dose <- 1L
  n_new_doses <- length(idx_new_dose)

  old_n_vacc_classes <- pars[[1]][[1]]$n_vacc_classes
  new_n_vacc_classes <- old_n_vacc_classes + n_new_vacc_classes
  old_idx_dose <- pars[[1]][[1]]$index_dose
  idx_booster <- old_idx_dose[length(old_idx_dose)] + idx_new_dose
  new_idx_dose <- c(old_idx_dose, idx_booster)
  new_n_vacc_classes <- pars[[1]][[1]]$n_vacc_classes + n_new_vacc_classes

  ## Easy updates: ## TOO: check if need anything else here
  update <- list(n_vacc_classes = new_n_vacc_classes,
                 n_doses = pars[[1]][[1]]$n_doses + n_new_doses,
                 vaccine_index_booster = idx_booster,
                 index_dose = new_idx_dose,
                 index_dose_inverse =
                   sircovid:::create_index_dose_inverse(new_n_vacc_classes,
                                                        new_idx_dose))

  inflate_pars <- function(p_i) {
    p_i[names(update)] <- update
    for (rel in grep("^rel_", names(p_i), value = TRUE)) {
      rel_old <- p_i[[rel]]
      if (!is.null(dim(rel_old))) {
        rel_new <- array(0, c(nrow(rel_old), ncol(rel_old),
                              nlayers(rel_old) + n_new_vacc_classes))
        rel_new[, , seq_len(nlayers(rel_old))] <- rel_old
        p_i[[rel]] <- rel_new
      }
    }
    rel <- "vaccine_progression_rate_base"
    rel_old <- p_i[[rel]]
    rel_new <- matrix(0, nrow(rel_old), ncol(rel_old) + n_new_vacc_classes)
    rel_new[, seq_len(ncol(rel_old))] <- rel_old
    p_i[[rel]] <- rel_new
    p_i
  }

  pars_new <- lapply(pars, lapply, inflate_pars)

  ## Then update the states given that:
  info_old <- info[[1]]$info
  info_new <- sircovid::lancelot$new(pars_new[[1]][[1]], 0, 1)$info()
  state_new <- lapply(state, sircovid::inflate_state_vacc_classes,
                      info_old, info_new)

  info_new <- lapply(info, function(x) {
    x$info <- info_new
    x
  })

  list(pars = pars_new, state = state_new, info = info_new)
}


simulate_one_pars_vaccination <- function(region, args, combined, n_strain) {
  priority_population <- sircovid::vaccine_priority_population(
    region,
    uptake = args$vaccine_uptake * args$vaccine_eligibility)

  pars <- combined$pars[, region]
  vaccine <- combined$vaccine[[region]]

  vaccine_index_dose2 <- pars[[1]]$index_dose[[2]]
  vaccine_progression_rate <- pars[[1]]$vaccine_progression_rate_base
  N_tot <- pars[[1]]$N_tot

  ## TODO: in the validation, if booster doses is non-empty, we should
  ## check that we have a model with boosters
  if (!is.null(args$vaccine_booster_daily_doses)) {
    args$vaccine_efficacy <-
      Map(cbind, args$vaccine_efficacy, args$vaccine_booster_efficacy)
    args$strain_vaccine_efficacy <-
      Map(cbind, args$strain_vaccine_efficacy,
          args$strain_vaccine_booster_efficacy)
    ## TODO: index_dose for boosters is assumed to be fourth column, when in
    ## fact it's now the fifth!! This needs fixing ASAP
    vaccine_index_booster <- pars[[1]]$index_dose[[3]] + 1
  } else {
    vaccine_index_booster <- NULL
  }

  mean_days_between_doses <- round(vaccine$mean_days_between_doses *
                                     args$vaccine_delay_multiplier)

  ## TODO: potentially a big problem here! We cannot fully disentangle booster
  ## eligibility from general vaccine eligibility, cna we? Especially, it seems
  ## the wya the latter is inputed (i.e. eligibility * uptake = calculation of
  ## priority population) seems to heavily affect actual boosters uptake! So
  ## we are potentially overestimating by quite a lot the effect of a
  ## boster vaccination programme!
  vaccine_schedule <- sircovid::vaccine_schedule_scenario(
    ## schedule_past got renamed to schedule_real in new params task
    schedule_past = vaccine$schedule_real,
    doses_future = args$vaccine_daily_doses[[region]],
    end_date = args$end_date,
    mean_days_between_doses = mean_days_between_doses,
    priority_population = priority_population,
    lag_groups = args$vaccine_lag_groups,
    lag_days = args$vaccine_lag_days,
    boosters_future = args$vaccine_booster_daily_doses[[region]],
    boosters_prepend_zero = TRUE,
    booster_proportion = args$vaccine_booster_eligibility)

  ## check boosters
  # > par(mfrow = c(3, 1))
  # > image(t(vaccine_schedule$doses[, 1, ]))
  # > image(t(vaccine_schedule$doses[, 2, ]))
  # > image(t(vaccine_schedule$doses[, 3, ]))

  ## TODO: potential placeholder for manually setting
  ## vaccine_schedule$doses[, 3, ] to zero if we don't want to give boosters
  ## to all age groups
  ## but probably better done inside sircovid::vaccine_schedule_scenario
  # > vaccine_schedule$doses[1:17, 3, 0] <- 0

  rel_list <- sircovid::modify_severity(
    args$vaccine_efficacy,
    args$strain_vaccine_efficacy,
    args$strain_severity_modifier)

  extra <- sircovid:::lancelot_parameters_vaccination(
    N_tot,
    pars[[1]]$dt,
    rel_susceptibility = rel_list$rel_susceptibility,
    rel_p_sympt = rel_list$rel_p_sympt,
    rel_p_hosp_if_sympt = rel_list$rel_p_hosp_if_sympt,
    rel_p_death = rel_list$rel_p_death,
    rel_infectivity = rel_list$rel_infectivity,
    vaccine_schedule = vaccine_schedule,
    vaccine_index_dose2 = vaccine_index_dose2,
    vaccine_index_booster = vaccine_index_booster,
    vaccine_progression_rate = vaccine_progression_rate,
    n_strains = n_strain,
    n_doses = pars[[1]]$n_doses)

  n_group <- nrow(rel_list$rel_p_sympt)
  n_vacc_class <- extra$n_vacc_classes

  rel_severity <- sircovid:::build_rel_param(extra$rel_p_death, n_strain,
                                             n_vacc_class, "rel_p_death")

  extra$rel_p_ICU <- extra$rel_p_R <-
    array(1, c(n_group, n_strain, n_vacc_class))

  extra$rel_p_ICU_D <- extra$rel_p_H_D <- extra$rel_p_W_D <-
    extra$rel_p_G_D <- rel_severity

  if (!is.null(args$strain_transmission)) {
    strain_params <- sircovid:::lancelot_parameters_strain(
      args$strain_transmission,
      sircovid::sircovid_date(args$strain_seed_date),
      args$strain_seed_rate[[region]],
      pars[[1]]$dt)
    strain_params$cross_immunity <- args$strain_cross_immunity
    strain_params$waning_rate <- rep(args$waning_rate, 19)
    extra <- c(extra, strain_params)
  }

  for (i in seq_along(pars)) {
    pars[[i]][names(extra)] <- extra
  }

  ## TODO: someone could explain why this is the check function wanted
  ## here, as it is not obvious.
  lapply(pars, sircovid::lancelot_check_severity)
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
    names(index_prob_strain) <- paste0(
      "prob_strain_", seq_along(index_prob_strain))
  } else {
    index_prob_strain <- NULL
  }

  index <- c(sircovid::lancelot_index(info)$state[keep],
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


## TODO: this would be heaps nicer if we saved the final rt into the object
setup_future_betas <- function(pars, rt_future, S, rt_type,
                               step_current, step_end, dt, seasonality,
                               R, prob_strain) {
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
    sircovid::lancelot_Rt(step_current, S[, i, drop = FALSE], pars[[i]],
                          type = rt_type, R = R[, i, drop = FALSE],
                          prob_strain = prob_strain[, i, drop = FALSE],
                          weight_Rt = FALSE)[[rt_type]][1])
  beta_rt_ratio <- beta[, , step_current] / rt

  for (region_index in seq_len(ncol(pars))) {
    r <- colnames(pars)[[region_index]]
    rt_future_r <- rt_future[rt_future$region == r, ]
    rt_future_r$step_start <- sircovid::sircovid_date(rt_future_r$date) / dt
    rt_future_r$step_end <- c(rt_future_r$step_start[-1L] - 1L, step_end)

    for (i in seq_rows(rt_future_r)) {
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
  beta_mult <- spim_calc_seasonality(date, seasonality_date_peak, seasonality)

  i <- seq(to = ncol(beta), length.out = n_beta_add)
  beta[, i] <- beta[, i] * rep(beta_mult, each = nrow(beta))

  for (i in seq_along(pars)) {
    pars[[i]]$beta_step <- beta[i, ]
  }

  pars
}


##' Seasonality
##'
##'@title Seasonality
##'
#' @param date date in sircovid date format
#' @param seasonality_date_peak date of peak in sircovid date format
#' @param seasonality peak to annual average difference
#' @return a value or vector of values with seasonal effect for each date
#'
#' @export
spim_calc_seasonality <- function(date, seasonality_date_peak, seasonality) {
  delta <- ((date - seasonality_date_peak) %% 365) / 365
  1 + cos(2 * pi * delta) * seasonality
}


## TODO: See spim_restart_initial_inflate_strain in ncov-outputs
## (rtm_inference_pmcmc_spim_fits2_restart_variant/support.R) for
## another shot at this which is stochastic.
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
        tmp_array <- round(prop[[regions[[r]]]] * new_state[, 1, , , ])
        new_state[, 2, , , ] <- new_state[, 2, , , ] + tmp_array
        new_state[, 1, , , ] <- new_state[, 1, , , ] - tmp_array
      } else {
        stop(sprintf("Unexpected dimensions (%d) in move_strain_compartment",
                     length(dim)))
      }
      state[info$index[[compartment[[i]]]], , r] <- array(new_state, tmp_dim)
    }
  }
  state
}


create_summary_state <- function(state, keep, dates) {
  summary_state <- state[keep, , , ]
  # aggregate over regions to get England peak hosp
  summary_state_england <-
    apply(summary_state[, , sircovid::regions("england"), ], c(1, 2, 4), sum)
  summary_state <- abind_quiet(list(summary_state,
                                    england = summary_state_england),
                               along = 3)

  # calculate number and date of peak hospital bed occupancy
  peak_hosp <-
    list(peak_hosp = apply(summary_state["hosp", , , ], c(1, 2), max),
         peak_hosp_date =
           apply(summary_state["hosp", , , ], c(1, 2), which.max))
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

  # reshape to add empty time dimension - this will make summarising easier
  abind_quiet(summary_state, along = 4)
}


simulate_rt <- function(steps, S, pars, critical_dates, voc_seeded,
                        R = NULL, prob_strain = NULL, no_seeding = FALSE,
                        prop_voc = NULL, weight_Rt = TRUE,
                        rt_type = c("Rt_general", "eff_Rt_general")) {
  dim_S <- dim(S)
  region_names <- dimnames(S)[[3]]
  ## TODO: add mcstate::array_combine(S, 2:3)
  dim(S) <- c(dim_S[[1]], dim_S[[2]] * dim_S[[3]], dim_S[[4]])

  if (!is.null(R)) {
    dim_R <- dim(R)
    dim(R) <- c(dim_R[[1]], dim_R[[2]] * dim_R[[3]], dim_R[[4]])
  }

  if (!is.null(prob_strain)) {
    ## Manually set prob strain to remove VOC if:
    ##  1. No VOC seeded in the combined; and
    ##  2. Not seeding in the simulation; and
    ##  3. Not moving any proportion in the simulation
    if (!voc_seeded && no_seeding && is.null(prop_voc)) {
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

  rt <- sircovid::lancelot_Rt_trajectories(
    steps, S, pars, type = rt_type,
    initial_step_from_parameters = FALSE,
    interpolate_critical_dates = critical_dates,
    interpolate_every = 7,
    interpolate_min = 3, R = R, prob_strain = prob_strain,
    weight_Rt = weight_Rt)

  ## Better format for the Rt values, will move into sircovid.
  ret_rt <- rt[rt_type]
  for (i in rt_type) {
    if (weight_Rt) {
      ret_rt[[i]] <- mcstate::array_reshape(t(rt[[i]]), 1, dim_S[2:3])
      dimnames(ret_rt[[i]]) <- list(NULL, region_names, NULL)
    } else {
      ret_rt[[i]] <- vapply(
        seq_len(ncol(ret_rt[[i]])),
        function(x) mcstate::array_reshape(t(rt[[i]][, x, ]), 1, dim_S[2:3]),
        array(0, dim_S[2:4])
      )
      dimnames(ret_rt[[i]]) <- list(NULL, region_names, NULL)
    }
  }

  ret_rt
}


simulate_calculate_vaccination <- function(state, index, vaccine_efficacy,
                                           booster_efficacy, n_strain,
                                           strain_vaccine_efficacy,
                                           strain_vaccine_booster_efficacy,
                                           strain_cross_immunity) {
  n_groups <- sircovid:::lancelot_n_groups()
  regions <- dimnames(state)[[3]]

  ## output the cumulative transitions between vaccine strata
  ## by age / vaccine stratum / region / over time
  n_vaccinated <- apply(state[names(index$n_vaccinated), , , , drop = FALSE],
                        c(1, 3, 4), mean)
  n_strata <- nrow(n_vaccinated) / n_groups
  n_vaccinated <-
    mcstate::array_reshape(n_vaccinated, 1L, c(n_groups, n_strata))

  ## output the number recovered in each vaccine stratum / region / over time
  ## R_raw: [age, strain, vaccine, particle, region, time]
  R_raw <- mcstate::array_reshape(
    state[names(index$R), , , , drop = FALSE],
    1L, c(n_groups, n_strain, n_strata))

  ## Take the sum over age
  R <- apply(R_raw, seq(2, 6), sum)
  ## take the mean over the particles
  R <- apply(R, c(1, 2, 4, 5), mean)
  ## R: [strain, vaccine, region, time]

  ## need to allow for imperfect cross-strain immunity
  ## R_strain_1: 100% of strain-level 1, 3, 4 +
  ##   args$cross_protection[2] * strain-level 2
  ## R_strain_2: 100% of strain-level 2, 3, 4 +
  ##   args$cross_protection[1] * strain-level 1
  calc_strain_immunity <- function(strain, R, strain_cross_immunity) {
    apply(R[-strain, , , ], seq(2, 4), sum) +
      R[strain, , , ] * strain_cross_immunity[strain]
  }

  R_strain <- list(strain_1 = calc_strain_immunity(2, R, strain_cross_immunity),
                   strain_2 = calc_strain_immunity(1, R, strain_cross_immunity))


  ## R_strain: [vaccine, region, time]

  # calculate the proportion protected given strain-specific vaccine efficacy
  # and cross-strain immunity

  n_protected_strain_1 <- calculate_n_protected(
    n_vaccinated, R_strain$strain_1, vaccine_efficacy, booster_efficacy)
  dimnames(n_protected_strain_1)[[2]] <- regions

  n_protected_strain_2 <- calculate_n_protected(
    n_vaccinated, R_strain$strain_2, strain_vaccine_efficacy,
    strain_vaccine_booster_efficacy)
  dimnames(n_protected_strain_1)[[2]] <- regions

  # Output number of first, second and booster doses

  idx_doses <- c("first_dose" = 1, "second_dose" = 3, "booster_dose" = 5)
  doses <- n_vaccinated[, idx_doses, , ]
  dimnames(doses)[2:3] <- list(names(idx_doses), regions)
  doses_inc <- aperm(apply(doses, c(1, 2, 3), diff), c(2, 3, 4, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-4], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")

  n_doses <- abind_quiet(doses, doses_inc, along = 2)

  list(n_vaccinated = n_vaccinated,
       n_protected = list(strain_1 = n_protected_strain_1,
                          strain_2 = n_protected_strain_2),
       n_doses = n_doses)
}


## TODO: overlap considerably with calculate_n_protected
## make R strain specific
calculate_n_protected <- function(n_vaccinated, R, vaccine_efficacy,
                                  booster_efficacy) {
  vp <- get_vaccine_protection(vaccine_efficacy, booster_efficacy)

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
  #   > sum_asr(V[a, s, r, t] * eff_severe[a, s])
  # - number protected after vaccination against infection
  #   > sum_asr(V[a, s, r, t] * eff_inf[a, s])
  # - number protected after infection (this is added to all of the above):
  #   > sum_sr(R[s, r, t])
  # - number protected after infection only:
  #   > sum_r(R[1, r, t]


  ## create array of number in each vaccine stratum
  n_strata <- ncol(n_vaccinated)

  # check arrays are conformable
  stopifnot(ncol(vp[[1]]) == n_strata)

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
    protected_against_death = sum_asr(c(vp$death) * V),
    ever_infected = sum_sr(R),
    ever_infected_unvaccinated = R[1, , ]
  )

  aperm(abind_quiet(ret, along = 3), c(3, 1, 2))
}




simulate_args_validate <- function(grid, vars, base, ignore, multistrain) {
  err <- intersect(ignore, names(vars))
  if (length(err) > 0) {
    stop("Names in ignore must not occur in vars: ",
         paste(squote(err), collapse = ", "))
  }

  err <- intersect(names(vars), names(base))
  if (length(err) > 0) {
    stop("Names in base must not occur in vars: ",
         paste(squote(err), collapse = ", "))
  }

  msg <- setdiff(names(grid), union(names(vars), ignore))
  if (length(msg) > 0) {
    stop("grid elements not found in vars: ",
         paste(squote(msg), collapse = ", "))
  }

  for (v in setdiff(names(grid), ignore)) {
    msg <- setdiff(grid[[v]], names(vars[[v]]))
    if (length(msg) > 0) {
      stop(sprintf("Missing elements in %s: %s (found %s)",
                   v, paste(squote(msg), collapse = ", "),
                   paste(names(vars[[v]]), collapse = ", ")))

    }
  }

  msg <- setdiff(simulate_args_names(multistrain),
                 union(names(vars), names(base)))
  if (length(msg) > 0) {
    stop("Required elements not found in vars or base: ",
         paste(squote(msg), collapse = ", "))
  }
}


## this is where we check that any dependencies among parameters will
## be satisfied. For example if we have booster daily doses then we
## must also have efficacy information. Other parameters depend on
## strains.
##
## if multistrain is FALSE, then we need all strain_ parameters to be NULL
simulate_validate_args1 <- function(args, regions, multistrain) {
  has_boosters <- !is.null(args$vaccine_booster_daily_doses)
  n_vacc_strata <- ncol(args$vaccine_efficacy[[1]])
  n_groups <- nrow(args$vaccine_efficacy[[1]])

  expected <- simulate_args_names(multistrain)

  msg <- setdiff(expected, names(args))
  if (length(msg) > 0) {
    stop("Missing expected values from args: ",
         paste(squote(msg), collapse = ", "))
  }

  assert_is(args$end_date, "Date")
  ## seed: null or raw
  assert_scalar_positive_integer(args$n_threads)
  assert_character(args$output_keep)
  assert_scalar_logical(args$output_rt)
  assert_scalar_logical(args$output_time_series)
  assert_scalar_logical(args$output_vaccination)
  assert_scalar_logical(args$output_state_by_age)
  assert_scalar_logical(args$output_weight_rt)
  match_value(args$rt_type, c("Rt_general", "eff_Rt_general"))

  assert_scalar_numeric(args$seasonality)
  assert_length(args$vaccine_uptake, n_groups)
  assert_length(args$vaccine_eligibility, n_groups)

  validate_rt_future(args$rt_future, regions)

  validate_vaccine_efficacy(args$vaccine_efficacy, n_groups, n_vacc_strata)
  validate_vaccine_doses(args$vaccine_daily_doses, regions,
                         "vaccine_daily_doses")

  assert_numeric(args$waning_rate)
  assert_length(args$waning_rate, 1)


  if (multistrain) {
    assert_is(args$strain_seed_date, "Date")
    validate_strain_seed_rate(args$strain_seed_rate, regions)
    validate_vaccine_efficacy(args$strain_vaccine_efficacy,
                              n_groups, n_vacc_strata)
    assert_length(args$strain_cross_immunity, 2)
    assert_numeric(args$strain_cross_immunity)
    validate_strain_severity_modifier(
      args$strain_severity_modifier)
  } else {
    for (i in grep("^strain_", names(args), value = TRUE)) {
      assert_is(args[[i]], "NULL", i)
    }
  }

  if (has_boosters) {
    validate_vaccine_doses(args$vaccine_booster_daily_doses, regions,
                           "vaccine_booster_daily_doses")
    validate_vaccine_efficacy(args$vaccine_booster_efficacy,
                              n_groups, 1)
  } else {
    for (v in c("vaccine_booster_daily_doses", "vaccine_booster_efficacy")) {
      assert_is(args[[v]], "NULL", v)
    }
  }

  if (multistrain && has_boosters) {
    validate_vaccine_efficacy(args$strain_vaccine_booster_efficacy,
                              n_groups, 1)
  } else {
    assert_is(args$strain_vaccine_booster_efficacy, "NULL")
  }
}


validate_vaccine_doses <- function(x, regions, name = deparse(substitute(x))) {
  assert_is(x, "list", name)
  msg <- setdiff(regions, names(x))
  if (length(msg) > 0) {
    stop(sprintf("Missing regions from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  for (r in regions) {
    el <- x[[r]]
    if (!is.numeric(el) || anyNA(el)) {
      stop(sprintf("%s$%s must be numeric and non-NA", name, r))
    }
    if (is.null(names(el))) {
      stop(sprintf("%s$%s must be named", name, r))
    }
    if (anyNA(suppressWarnings(as.Date(names(el))))) {
      stop(sprintf("names of %s$%s must be dates (YYYY-MM-DD)", name, r))
    }
  }
}


validate_strain_seed_rate <- function(x, regions,
                                      name = deparse(substitute(x))) {
  assert_is(x, "list", name)
  msg <- setdiff(regions, names(x))
  if (length(msg) > 0) {
    stop(sprintf("Missing regions from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  for (r in regions) {
    assert_scalar_numeric(x[[r]], sprintf("%s:%s", name, r))
  }
}


validate_vaccine_efficacy <- function(x, n_groups, n_vacc_strata,
                                      name = deparse(substitute(x))) {
  expected <- c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt",
                "rel_infectivity", "rel_p_death")
  if (!setequal(names(x), expected)) {
    stop(sprintf("Invalid names for %s, expected %s",
                 name, paste(squote(expected), collapse = ", ")))
  }
  ok <- vlapply(x, function(e)
    identical(dim(e), as.integer(c(n_groups, n_vacc_strata))))
  if (!all(ok)) {
    stop(sprintf("All elements of %s must have size %d x %d",
                 name, n_groups, n_vacc_strata))
  }
}


## TODO: this data structure can be replaced by a single number
## ([[3]]$rep_p_hosp_if_sympt) and is generally horrific. Can we
## simplify this please?
validate_strain_severity_modifier <- function(x,
                                              name = deparse(substitute(x))) {
  assert_length(x, 4)
  for (i in seq_along(x)) {
    el <- x[[i]]
    expected <- c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt",
                  "rel_infectivity", "rel_p_death")
    if (!setequal(names(el), expected)) {
      stop(sprintf("Invalid names for %s[[%d]] expected %s",
                   name, i, paste(squote(expected), collapse = ", ")))
    }
    for (v in names(el)) {
      assert_scalar_numeric(el[[v]], sprintf("%s[[%d]]$%s", name, i, v))
    }
  }
}

## TODO: We might filter off dates in the past, and/or require that
## they are gone here. I (Rich) remember that being nasty in getting
## future betas done, so be careful.
validate_rt_future <- function(x, regions, name = deparse(substitute(x))) {
  assert_is(x$date, "Date", sprintf("%s$date", name))
  msg <- setdiff(regions, x$region)
  if (length(msg) > 0) {
    stop("No rt values found for regions: ",
         paste(squote(msg), collapse = ", "))
  }
  assert_numeric(x$Rt, sprintf("%s$Rt", name))
  assert_numeric(x$Rt_sd, sprintf("%s$Rt_sd", name))
}


simulate_extract_age_class_state <- function(state, index) {
  n_groups <- sircovid:::lancelot_n_groups()

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

    if (ncol(x) == 3) {
      colnames(x) <- c("unvaccinated", "partial_protection", "full_protection")
    } else  if (ncol(x) == 5) {
      colnames(x) <- c("unvaccinated", "partial_protection", "full_protection",
                       "waned_protection", "booster")
    }

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
    ret <- apply(abind_quiet(res, along = 5), c(1, 3, 4, 5), mean)

    # [age, vaccine status, region, time]
    round(aperm(ret, c(4, 1, 2, 3)))
  }

  lapply(arrays, f)
}


##' Create expanded run grid for simulation
##'
##' @title Create expanded run grid
##'
##' @param ... named variables to expand over, should be in
##'  `simulation_central_analysis`, omitted variables will take central value
##' @param full_run If `TRUE` saves trajectories for expanded scenarios, this
##'   should very rarely be `TRUE`, change default with care as this will lead
##'   to massive objects
##' @param prefix prefix for analysis name, prefixes row number
##'
##' @return A grid of scenarios to run
##' @export
spim_expand_grid <- function(..., full_run = FALSE, prefix = "Grid_") {

  central <- simulation_central_analysis(full_run)

  actual <- names(list(...))
  expected <- setdiff(colnames(central), c("RUN", "full_run"))
  mtc <- is.na(match(actual, expected))
  if (any(mtc)) {
    stop(sprintf("Unexpected variables(s) %s", str_collapse(actual[mtc])))
  }

  tidyr::expand_grid(
    RUN = TRUE,
    full_run,
    ...
  ) %>%
    ## adds central for missing variables
    tidyr::expand_grid(dplyr::select(central, -names(.)), .) %>%
    ## add analysis name
    dplyr::mutate(analysis = paste0(prefix, seq_rows(.)))
}


##' Create grid of scenarios to run for simulation
##'
##' @title Create scenario run grid
##'
##' @param scenarios Scenarios to run simulation over
##' @param csv Path of csv to load run grid from
##' @param expand_grid Optional large grid of scenarios such as from
##'   [spim_expand_grid]
##' @param force_central If `TRUE` (default) then central analysis is always
##'  included as specified in `simulation_central_analysis`. This should rarely
##'  be `FALSE` as often required for basic checking plots.
##' @param set_strain_params If `TRUE` automatically sets strain parameters
##'   `strain_cross_immunity`, `strain_severity_modifier` and
##'   `waning_rate`, as equivalent to `strain_vaccine_efficacy`
##' @param multistrain If `FALSE` then removes all columns related to a second
##'   strain
##' @param analyses If not `NULL` then filters grid by given analyses
##'
##' @return A grid of scenarios to run
##' @export
spim_run_grid <- function(scenarios, csv = NULL, expand_grid = NULL,
                          force_central = TRUE, set_strain_params = TRUE,
                          multistrain = TRUE, analyses = NULL) {

  if (is.null(csv) && is.null(expand_grid) && !force_central) {
    stop("At least one of 'csv', 'expand_grid', 'force_central' must be
    non-NULL/TRUE")
  }

  run_grid <- simulation_central_analysis(TRUE, multistrain)
  if (!force_central) {
    run_grid <- run_grid[-1, ]
  }

  if (!is.null(csv)) {
    csv_grid <- read.csv(csv)
    if (!multistrain) {
      csv_grid <- csv_grid %>% select(-starts_with("strain_"))
    }
    run_grid <- rbind(run_grid, csv_grid)
  }

  if (!is.null(expand_grid)) {
    run_grid <- rbind(run_grid, expand_grid)
  }

  if (multistrain && set_strain_params) {
    run_grid <- run_grid %>%
      dplyr::mutate(strain_cross_immunity = strain_vaccine_efficacy,
                    strain_severity_modifier = strain_vaccine_efficacy,
                    waning_rate = strain_vaccine_efficacy)
  }

  if (!is.null(analyses)) {
    run_grid <- run_grid %>%
      dplyr::filter(analysis %in% analyses)
  }

  run_grid %>%
    dplyr::filter(RUN) %>%
    ## expand over scenarios
    tidyr::expand_grid(scenario = scenarios) %>%
    ## de-duplicate
    dplyr::distinct(across(!any_of("analysis")), .keep_all = TRUE) %>%
    ## set cross immunity and modifier
    dplyr::mutate(rt_future =
                    paste(scenario, adherence_to_baseline_npis, sep = ": "))
}

##' Calculate SHAPs from a tidy summary of predictions over various features.
##'  SHAPs calculated as the expected difference in predicted states with and
##'  without the given feature of interest (over all feature levels).
##'
##' @title Calculate SHAPs over predicted states
##'
##' @param summary A tidy summary object such as that returned by
##'   `create_summary`
##' @param feats Features to calculate SHAPS for. If NULL then uses default
##'   selection returned by `spim_simulation_predictors`
##'
##' @export
spim_simulation_shaps <- function(summary, feats = NULL) {

  if (is.null(feats)) {
    feats <- spim_simulation_predictors(summary)
  }

  states <- unique(summary$state)
  out <- set_names(vector("list", length(states)), states)

  for (State in states) {
    state_df <- summary %>%
      dplyr::filter(state == State) %>%
      dplyr::select(`50%`, feats)

    shaps <- set_names(vector("list", length(feats)), feats)
    for (i in seq_along(shaps)) {
      lvls <- unique(state_df[[feats[[i]]]])
      lvls <- set_names(numeric(length(lvls)), lvls)

      for (j in names(lvls)) {
        for (k in names(lvls)) {
          if (!identical(j, k)) {
            with <- state_df %>%
              dplyr::filter(!!as.symbol(feats[[i]]) == j) %>%
              dplyr::select(`50%`) %>%
              unlist()
            without <- state_df %>%
              dplyr::filter(!!as.symbol(feats[[i]]) == k) %>%
              dplyr::select(`50%`) %>%
              unlist()

            ## only compare possible scenarios
            which <- intersect(names(with), names(without))

            lvls[[j]] <- mean(c(lvls[[j]], mean(with[which] - without[which])))
          }
        }
      }

      shaps[[i]] <- lvls
    }

    mshaps <- reshape2::melt(shaps)
    mshaps$lvl <- unlist(lapply(shaps, names))
    out[[State]] <- mshaps
  }

  mout <- do.call(rbind, out)
  mout$state <- vapply(strsplit(rownames(mout), ".", TRUE), "[[",
                       character(1), 1)
  rownames(mout) <- NULL
  colnames(mout)[[2]] <- "Var"
  mout
}

##' Names of 'predictive' variables from a tidy simulation summary object, i.e.
##'  those variables which: i) are not purely informative;
##'  ii) impact upon predictions; iii) are not outcome variables.
##'
##' @title Return predictive simulation variables
##'
##' @param summary A tidy summary object such as that returned by
##'   `create_summary`
##'
##' @export
spim_simulation_predictors <- function(summary) {
  vars <- names(which(vapply(summary, function(i)
    length(unique(i)) > 1, logical(1))))
  vars <- setdiff(vars, c(## not needed
    "analysis", "full_run", "state", "rt_future",
    ## taken into account in scenario
    "adherence_to_baseline_npis",
    ## outcomes
    "2.5%", "50%", "97.5%",
    ## these two are identical to
    ## vaccine_efficacy_strain_2
    "strain_severity_modifier", "rel_strain_modifier"))

  vars
}


##' Prepare csv with NPI keys for simulation task
##'
##' @title Prepare NPI key csv for simulation
##'
##' @param path File path to a csv containing columns: nation
##'  (at least one of england, scotland, wales, northern ireland), npi
##'  (name of NPI step for associated schedule), Rt (mean value for Rt),
##'  Rt_sd (standard deviation for Rt), adherence (low, central, or high)
##'
##' @param schools Should be one of 'open' or 'closed' and specifies which
##'  state values in the csv correspond to
##'
##' @param schools_modifier Amount to modify first strain Rt (excl. immunity).
##'  If `schools = 'open'` then the modifier will be subtracted
##'  from the given csv to create `schools = "closed"` scenarios, otherwise
##'  added - `abs(schools_modifier)` is applied internally
##'
##' @param country If `"england"` then all other nations filtered, otherwise
##'  no filtering.
##'
##' @param gradual_start,gradual_end,gradual_steps If all three of these are
##'  non-NULL then creates a scenario of gradual transition from Rt in scenario
##'  `gradual_start` to Rt in `gradual_end` in `gradual_steps` number of steps
##'
##' @param overwrite_central_adherence
##'  If non-NULL then a list with names corresponding to `npi` and values to
##'  overwrite the csv in the central/low/high adherence scenario. This may be
##'  be useful if automating some, but not all, of the Rt values
##'
##' @param overwrite_low_adherence,overwrite_high_adherence As for
##'   `overwrite_central_adherence`
##'
##' @param overwrite_mtp_adherence As for `overwrite_central_adherence`
##'
##' @param npi_key Optional data.frame, instead of reading from `path`
##'
##' @param modify_gradual Logical, stating whether to adjust gradual `npi`
##'   values accounting for schools open/close effect
##'
##' @return tibble for passing to [spim_prepare_rt_future]
##'
##' @export
spim_prepare_npi_key <- function(schools, schools_modifier, country,
                                 path = NULL,
                                 gradual_start = NULL, gradual_end = NULL,
                                 gradual_steps = NULL,
                                 overwrite_central_adherence = NULL,
                                 overwrite_low_adherence = NULL,
                                 overwrite_high_adherence = NULL,
                                 overwrite_mtp_adherence = NULL,
                                 npi_key = NULL, modify_gradual = FALSE) {

  schools_modifier <- abs(schools_modifier)

  all_schools <- c("open", "closed")
  stopifnot(schools %in% all_schools)
  stopifnot(length(c(gradual_start, gradual_end, gradual_steps)) %in% c(0, 3))
  stopifnot(sum(is.null(path), is.null(npi_key)) == 1)

  if (is.null(npi_key)) {
    npi_key <- read_csv(path)
  } else {
    stopifnot(identical(colnames(npi_key),
                        c("nation", "npi", "Rt", "Rt_sd", "scenario")))
  }

  npi_key <- npi_key %>%
    dplyr::filter(nation == case_when(country == "england" ~ "england",
                                      TRUE ~ nation)) %>%
    dplyr::mutate(npi = sprintf("%s_schools_%s", npi, schools))

  nations <- unique(npi_key$nation)

  for (ad in unique(npi_key$scenario)) {
    for (n in unique(npi_key$npi)) {
      if (sum(npi_key$npi == n & npi_key$scenario == ad) == 0) {
        tmp <- npi_key[npi_key$npi == n & npi_key$scenario == "central", ]
        tmp$scenario <- ad
        npi_key <- rbind(npi_key, tmp)
      }
    }
  }

  npi_key <- dplyr::bind_rows(
    npi_key,
    npi_key %>%
      dplyr::mutate(npi = gsub(schools, setdiff(all_schools, schools), npi),
                    Rt = case_when(schools == "open" ~ Rt - schools_modifier,
                                   schools == "closed" ~ Rt + schools_modifier))
  )

  if (!is.null(gradual_start)) {

    gradualise <- function(start, end, steps, ad) {
      lapply(unique(npi_key$nation), function(nat) {
        to <- npi_key %>%
          dplyr::filter(npi == end, scenario == ad, nation == nat) %>%
          dplyr::select(Rt) %>%
          as.numeric()
        from <- npi_key %>%
          dplyr::filter(npi == start, scenario == ad, nation == nat) %>%
          dplyr::select(Rt) %>%
          as.numeric()
        sd <- npi_key %>%
          dplyr::filter(npi == end, scenario == ad, nation == nat) %>%
          dplyr::select(Rt_sd) %>%
          as.numeric()
        steps <- round(seq.int(from, to, length.out = steps + 1)[2:steps], 3)

        npi <- sprintf("p%d_%s", seq_along(steps - 1), end)

        ## check if we've accounted for open or closed
        end_open <- grepl("open", end)

        if (end_open) {
          npi <- c(npi, gsub("open", "closed", npi))
        } else {
          npi <- c(npi, gsub("closed", "open", npi))
        }

        ## if modifying gradual then adjust with schools_modifier else
        ##  use same values
        if (modify_gradual) {
          if (end_open) {
            steps <- c(steps, steps - schools_modifier)
          } else {
            steps <- c(steps, steps + schools_modifier)
          }
        }

        data.frame(nation = nat, npi = npi, Rt = steps, Rt_sd = sd,
                   scenario = ad)
      }) %>%
        dplyr::bind_rows()
    }

    f <- function(ad) {
      gradualise(gradual_start, gradual_end, gradual_steps, ad)
    }
    npi_key <- npi_key %>%
      dplyr::bind_rows(
        dplyr::bind_rows(
          lapply(npi_key[npi_key$npi == gradual_end, "scenario"], f)))
  }

  npi_key %>%
    dplyr::arrange(scenario, nation, npi) %>%
    dplyr::mutate(Rt = round(Rt, 3)) %>%
    `rownames<-`(NULL)

}


##' Prepare Rt future grid for simulation task
##'
##' @title Prepare Rt future for simulation
##'
##' @param npi_key data.frame, rather than reading from `path`
##'
##' @param start_date Start date of simulation, all changes in schedule before
##'  this date are removed.
##'
##' @param end_date End date of simulation, all changes in schedule after
##'  this date are removed.
##'
##' @param path File path to a csv containing columns: nation
##'  (at least one of england, scotland, wales, northern ireland), scenario
##'  (scenario ID), year (YYYY), month (MM), (DD), npi
##'  (NPI ID to match npi_key$npi)
##'
##' @param schedule data.frame alternative to csv in `path` except with
##'  `date` column already create (and not year, month, day columns)
##'
##' @return tibble which is essentially a column binds of npi_key and
##'  Rt_schedule datasets after cleaning
##'
##' @export
spim_prepare_rt_future <- function(npi_key, start_date, end_date,
                                   path = NULL, schedule = NULL) {

  stopifnot(sum(is.null(path), is.null(schedule)) == 1)
  if (!is.null(path)) {
    schedule <-
      read_csv(path) %>%
      dplyr::mutate(date = as.Date(sprintf("%s-%s-%s", year, month, day))) %>%
      dplyr::select(-day, -month, -year)
  }

  schedule <- schedule %>%
    dplyr::filter(
      nation %in% unique(npi_key$nation),
      ## remove all dates after the end date and before the start date
      date >= as.Date(start_date),
      date <= as.Date(end_date)
    ) %>%
    dplyr::left_join(npi_key, by = c("nation", "npi")) %>%
    dplyr::mutate(key = paste(scenario_key, scenario, sep = ": "))

  schedule_celtic <- schedule %>%
    dplyr::filter(nation != "england") %>%
    dplyr::mutate(region = nation)

  schedule %>%
    dplyr::filter(nation == "england") %>%
    tidyr::expand_grid(region = sircovid::regions("england")) %>%
    dplyr::bind_rows(schedule_celtic)
}


##' Find little r from big R
##'
##' @title Find daily little r values from big R
##'
##' @param summary Simulation summary object
##' @param dates Dates for which to compute little r
##' @param scenarios Scenarios for which to compute little r . If `NULL`
##'  computed over all scenarios
##' @param analyses Analyses for which to compute little r . If `NULL`
##'  computed over all analyses
##' @param reg Region to filter summary object by
##' @param wide If `FALSE` (default) returns results in long format as numeric,
##'  otherwise returns wide format with dates as columns and scenarios/analyses
##'  as rows and entries are given as `central (low, high)`
##'
##' @return data.frame with daily little r at given dates, scenarios, analyses
##'
##' @export
spim_rejuvenatoR <- function(summary, dates, scenarios = NULL, analyses = NULL,
                             reg = "england", wide = FALSE) {
  ## generation time distribution from STM paper
  ## https://dx.doi.org/10.1126/scitranslmed.abg4262
  if (is.null(scenarios)) {
    scenarios <- unique(summary$state$scenario)
  }
  if (is.null(analyses)) {
    analyses <- unique(summary$state$analysis)
  }

  t <- seq(0, 30, 1)
  mean_gt <- 6.7
  sd_gt <- 3.5
  w <- EpiEstim::discr_si(t, mean_gt, sd_gt)

  ### grid of r values
  r_grid <- seq(-0.5, 0.5, 0.0001)
  R_grid <- epitrix::r2R0(r = r_grid, w = w)
  r_R_corresp <- data.frame(r = r_grid, R = R_grid)

  find_r_from_R <- function(R) {
    if (R < min(r_R_corresp$R) || R > max(r_R_corresp$R)) {
      stop("R value outside of grid, expand grid further")
    }
    r_R_corresp$r[which.min(abs(r_R_corresp$R - R))]
  }

  obj <- summary$state %>%
    dplyr::filter(state == "eff_Rt_general",
                  region == reg,
                  date %in% dates,
                  analysis %in% analyses,
                  scenario %in% scenarios)

  out <- apply(obj, 1, function(x) {
    x <- vapply(as.numeric(c(x[["50%"]], x[["2.5%"]], x[["97.5%"]])),
                find_r_from_R, numeric(1))
  }) %>%
    `rownames<-`(c("50%", "2.5%", "97.5%")) %>%
    t() %>%
    data.frame(obj$date, obj$analysis, obj$scenario) %>%
    `colnames<-`(c("50%", "2.5%", "97.5%", "date", "analysis", "scenario")) %>%
    dplyr::arrange(date, scenario, analysis)

  if (wide) {
    out <- out %>%
      dplyr::mutate(r = sprintf("%#.4f (%#.4f, %#.4f)",
                                `50%`, `2.5%`, `97.5%`)) %>%
      dplyr::select(date, scenario, analysis, r) %>%
      tidyr::pivot_wider(names_from = date, values_from = r)
  }

  out
}


#' Calculate when all second doses given out
#'
#' @title Calculate final doses date
#' @param summary Simulation summary object
#'
#' @export
spim_simulate_complete_doses <- function(summary) {
  summary$n_doses %>%
    dplyr::filter(state == "second_dose_inc",
                  region == "england") %>%
    dplyr::group_by(across(-c(value, group))) %>%
    dplyr::summarise(mean = sum(value)) %>%
    dplyr::filter(mean < 5e3) %>%
    dplyr::filter(date == min(date)) %>%
    dplyr::ungroup() %>%
    dplyr::select(analysis, scenario, date, mean)
}


#' Save NPI key with quantiles, mean, and standard deviation
#' @title Save NPI key
#' @param npi_key Output from [spim_prepare_npi_key]
#' @param filename File to save NPI key to
#' @export
spim_write_npi_key <- function(npi_key, filename) {
  unique(npi_key$nation) %>%
    lapply(function(i) {
      key <- npi_key %>%
        dplyr::filter(nation == i) %>%
        dplyr::select(-nation)

      npi_pars <- mapply(function(mean, sd) {
        dist <- distr6::dstr("Lognormal", mean = mean, sd = sd)
        unlist(c(q2.5 = dist$quantile(0.025),
                 q97.5 = dist$quantile(0.975),
                 meanlog = dist$parameters("meanlog")$value,
                 sdlog = dist$parameters("sdlog")$value))
      }, mean = key$Rt, sd = key$Rt_sd)
      npi_pars <- cbind(key, signif(t(npi_pars), 3))
      npi_pars$region <- i
      npi_pars
    }) %>%
    dplyr::bind_rows() %>%
    write.csv(filename, row.names = FALSE)
}


#' Write first and second dose uptake at a given date
#' @title Save dose uptake at a given date
#' @param report_date Date to calculate uptake
#' @param doses Output from [spim_calculate_doses]
#' @param filename File to save uptake to
#' @export
spim_write_uptake <- function(report_date, doses, filename) {
  doses %>%
    dplyr::filter(state %in% c("state_first_dose", "state_second_dose"),
                  date == report_date) %>%
    dplyr::select(group, state, analysis, prop) %>%
    arrange(state, analysis, group, prop) %>%
    write.csv(filename, row.names = FALSE)
}


#' Get variables used in simulation tasks
#' @title Get simulation variables
#' @export
spim_simulation_vars <- function() {
  c(
    "seasonality",
    "rt_future",
    "vaccine_daily_doses",
    "vaccine_booster_daily_doses",
    "vaccine_efficacy",
    "vaccine_booster_efficacy",
    "vaccine_booster_eligibility",
    "vaccine_eligibility",
    "vaccine_uptake",
    "vaccine_lag_groups",
    "vaccine_lag_days",
    "strain_transmission",
    "strain_seed_rate",
    "strain_vaccine_efficacy",
    "strain_initial_proportion",
    "strain_vaccine_booster_efficacy",
    "strain_cross_immunity",
    "strain_severity_modifier",
    "waning_rate"
  )
}
