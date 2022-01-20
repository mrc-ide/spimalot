
##' Process fit data for `lancelot` model
##'
##' @title Process a fit for `lancelot` model
##' @param samples The pmcmc_samples object from [mcstate::pmcmc]
##'
##' @param parameters The parameter specification
##'   ([spimalot::spim_pars_pmcmc_load])
##'
##' @param data Data sets used in fitting
##'
##' @export
spim_fit_process <- function(samples, parameters, data) {
  region <- samples$info$region

  ## This is just the info/prior/proposal + base of the parameter used
  ## to fit the model. Below we will also create a refitted
  ## parameters_new object with new proposal matrix
  parameters_raw <- parameters$raw
  parameters_raw$base <- parameters$base

  message("Computing restart information")
  restart <- fit_process_restart(samples, parameters_raw)
  samples$restart <- NULL

  samples$trajectories$date <-
    samples$trajectories$step / samples$trajectories$rate

  ## The Rt calculation is slow and runs in serial; it's a surprising
  ## fraction of the total time.
  message("Computing Rt")
  rt <- calculate_lancelot_Rt(samples, TRUE)
  variant_rt <- calculate_lancelot_Rt(samples, FALSE)

  ## TODO: someone needs to document what this date is for (appears to
  ## filter trajectories to start at this date) and when we might
  ## change it.
  message("Preparing onward simulation object")
  start_date_sim <- "2021-06-01"
  simulate <- create_simulate_object(samples, start_date_sim, samples$info$date)

  ## Reduce trajectories in samples before saving
  message("Reducing trajectories")
  samples <- reduce_trajectories(samples)

  message("Computing parameter MLE and covariance matrix")
  parameters_new <- spim_fit_parameters(samples, parameters_raw)

  ## This can possibly merge in with the initial restart processing if
  ## we're careful now.  We need to have computed Rt first is the only
  ## trick.
  if (!is.null(restart)) {
    ## When adding the trajectories, we might as well strip them down
    ## to the last date in the restart
    restart_date <- max(restart$state$time)
    i <- samples$trajectories$date <= restart_date

    restart$parent <- list(
      trajectories = trajectories_filter_time(samples$trajectories, i),
      rt = rt_filter_time(rt, i, samples$info$multiregion),
      data = data,
      prior = parameters_raw$prior)
  }

  ## This list returns:
  ##
  ## 1. 'fit' list of;
  ## samples - reduced mcstate trajectory samples
  ## rt - calculated spimalot lancelot Rt value (calculate_lancelot_Rt output)
  ## variant_rt - variant specific Rt
  ## simulate - object containing spimalot objects (used for onward simulation)
  ##            such as; thinned samples, start date of sim, samples date.
  ## parameters - parameter MLE and covariance matrix
  ## data - fitted and full data
  ##
  ## 2. 'restart' list of;
  ## restart information from spimalot parent objects (i.e trajectories, prior,
  ##   rt, etc.)
  list(
    fit = list(samples = samples,
               rt = rt,
               variant_rt = variant_rt,
               simulate = simulate,
               parameters = parameters_new,
               ## NOTE: fit$data$fitted is assumed to exist by the restart
               data = list(fitted = data$fitted, full = data$full)),
    restart = restart)
}


create_simulate_object <- function(samples, start_date_sim, date) {
  start_date_sim <- sircovid::sircovid_date(start_date_sim)
  fit_dates <- samples$trajectories$date
  idx_dates <- (fit_dates >= start_date_sim) &
    (fit_dates <= sircovid::sircovid_date(date))
  date <- fit_dates[idx_dates]

  state_keep <- c("deaths", "deaths_comm", "deaths_hosp", "admitted",
                  "diagnoses", "infections", "hosp", "icu")
  state_full <- samples$trajectories$state

  if (samples$info$multiregion) {
    region <- samples$info$region
    state_by_age <- lapply(region, function(r)
      extract_age_class_state(state_full[, r, , idx_dates]))
    names(state_by_age) <- region

    ## thin trajectories, but here we already have a regional dimension
    state <- state_full[state_keep, , , idx_dates]
    ## However, we have to push it out
  } else {
    state_by_age <- extract_age_class_state(state_full[, , idx_dates])

    ## thin trajectories and reshape to add a regional dimension:
    state <- state_full[state_keep, , idx_dates]
    n_samples <- ncol(state_full)
    ## TODO: this needs to be moved to be (1, n_samples) (swapping
    ## dimensions) to match the multiregion filter.
    state <- mcstate::array_reshape(state, i = 2, c(n_samples, 1))
  }

  list(date = date, state = state, state_by_age = state_by_age)
}


calculate_lancelot_Rt <- function(samples, weight_Rt) {
  step <- samples$trajectories$step
  info <- samples$info$info
  epoch_dates <- samples$info$epoch_dates

  state <- samples$trajectories$state
  pars <- samples$pars
  transform <- samples$predict$transform

  multiregion <- samples$info$multiregion

  if (multiregion) {
    ## pars, state, step, info, epoch_dates
    ret <- lapply(samples$info$region, function(r)
      calculate_lancelot_Rt_region(pars[, , r], state[, r, , ], transform[[r]],
                                   step, info, epoch_dates, weight_Rt))
    names(ret) <- samples$info$region
  } else {
    ret <- calculate_lancelot_Rt_region(pars, state, transform,
                                        step, info, epoch_dates, weight_Rt)
  }
  ret
}


calculate_lancelot_Rt_region <- function(pars, state, transform,
                                         step, info, epoch_dates, weight_Rt) {
  index_S <- grep("^S_", rownames(state))
  index_R <- grep("^R_", rownames(state))
  index_ps <- grep("^prob_strain", rownames(state))

  S <- state[index_S, , , drop = FALSE]
  R <- state[index_R, , , drop = FALSE]
  prob_strain <- state[index_ps, , , drop = FALSE]

  pars_model <- lapply(seq_len(length(info)), function(j)
    lapply(seq_rows(pars), function(i)
      transform(pars[i, ])[[j]]$pars))

  pars_model[[1]]$steps_per_day
  dates <- step / 4

  ## TODO: currently we'll just deal with multistage here, but it would
  ## be good to adapt the sircovid function to deal with multistage
  ## parameters
  rt <- list(step = numeric(0),
             date = numeric(0),
             beta = numeric(0),
             eff_Rt_all = numeric(0),
             eff_Rt_general = numeric(0),
             Rt_all = numeric(0),
             Rt_general = numeric(0))

  for (i in seq_len(length(epoch_dates) + 1L)) {
    if (i == 1) {
      dates1 <- which(dates <= epoch_dates[1])
      initial_step_from_parameters <- TRUE
    } else  if (i <= length(epoch_dates)) {
      dates1 <- which(dates > epoch_dates[i - 1] & dates <= epoch_dates[i])
      initial_step_from_parameters <- FALSE
    } else {
      dates1 <- which(dates > epoch_dates[i - 1])
      initial_step_from_parameters <- FALSE
    }

    if (length(dates1) == 0) {
      next
    }

    n_strains <- pars_model[[i]][[1]]$n_real_strains
    n_vacc_classes <- pars_model[[i]][[1]]$n_vacc_classes

    suffix <- paste0("_", c(sircovid:::sircovid_age_bins()$start, "CHW", "CHR"))
    S_nms <- get_names("S", list(n_vacc_classes), suffix)

    step1 <- step[dates1]
    S1 <- S[S_nms, , dates1, drop = FALSE]

    if (n_strains == 1) {
      R1 <- NULL
      prob_strain1 <- NULL
    } else {
      R_nms <- get_names("R",
                         list(S = n_real_strains, V = n_vacc_classes),
                         suffix)
      R1 <- R[R_nms, , dates1, drop = FALSE]
      prob_strain1 <- prob_strain[, , dates1, drop = FALSE]
    }

    if (!(n_strains == 1 && !weight_Rt)) {
      rt1 <- sircovid::lancelot_Rt_trajectories(
        step1, S1, pars_model[[i]],
        initial_step_from_parameters = initial_step_from_parameters,
        shared_parameters = FALSE, R = R1, prob_strain = prob_strain1,
        weight_Rt = weight_Rt)
      for (nm in names(rt)) {
        if (length(rt[[nm]]) == 0) {
          rt[[nm]] <- rt1[[nm]]
        } else if (length(dim(rt1[[nm]])) == 2) {
          rt[[nm]] <- rbind(rt[[nm]], rt1[[nm]])
        } else {
          rt[[nm]] <- abind1(rt[[nm]], rt1[[nm]])
        }
      }
    }
  }

  class(rt) <- c("Rt_trajectories", "Rt")
  rt
}


## All the functions below here have awful names, but none are
## exported so we can tidy this up later.
extract_age_class_state <- function(state) {
  n_groups <- sircovid:::lancelot_n_groups()
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
  remove_strings <- c("prob_strain", "S_", "R_", "I_weighted_", "D_hosp_",
                      "D_all_", "diagnoses_admitted_", "cum_infections_disag_",
                      "cum_n_vaccinated")
  re <- sprintf("^(%s)", paste(remove_strings, collapse = "|"))

  multiregion <- samples$info$multiregion
  state <- samples$trajectories$state
  keep <- !grepl(re, rownames(state))

  if (multiregion) {
    samples$trajectories$state <- samples$trajectories$state[keep, , , ]
  } else {
    samples$trajectories$state <- samples$trajectories$state[keep, , ]
  }

  ## TODO: this should be moved into a function with a more
  ## appropriate name

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
  if (length(dim(trajectories$state)) == 3) {
    trajectories$state <- trajectories$state[, , i, drop = FALSE]
  } else {
    trajectories$state <- trajectories$state[, , , i, drop = FALSE]
  }
  trajectories
}


rt_filter_time <- function(rt, i, multiregion) {
  if (multiregion) {
    lapply(rt, rt_filter_time, i, FALSE)
  } else {
    ret <- lapply(rt, function(x) x[i, , drop = FALSE])
    class(ret) <- class(rt)
    ret
  }
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
  keep <- function(x, region) {
    is.na(x) | x %in% region
  }

  region <- samples$info$region
  info <- parameters$info[keep(parameters$info$region, region), ]
  rownames(info) <- NULL

  prior <- parameters$prior[keep(parameters$prior$region, region), ]
  rownames(prior) <- NULL

  ## This will certainly want moving elsewhere, possibly mcstate, as
  ## it's very fiddly.  There's no good reason it needs to be here
  ## (see https://github.com/mrc-ide/mcstate/issues/189)
  if (samples$info$multiregion) {
    i <- which.max(rowSums(samples$probabilities[, "log_posterior", ]))
    initial <- samples$pars[i, , ]

    nms_fixed <- samples$info$pars$fixed
    nms_varied <- samples$info$pars$varied

    ## Spreading these back out is pretty nasty; we need to work out
    ## which are nested and which are not; take the fixed pars from
    ## the first region, then work out the mapping for the rest.
    info$key <- paste(info$name, info$region, sep = "\r")
    info$initial[match(nms_fixed, info$name)] <- unname(initial[nms_fixed, 1])

    ## Not the fastest, but the clearest...
    for (r in region) {
      for (nm in nms_varied) {
        info$initial[info$region == r & info$name == nm] <- initial[nm, r]
      }
    }

    ## For the covariance we need to split things up carefully:
    nms_all <- unique(info$name)
    n_all <- length(nms_all)
    n_fixed <- length(nms_fixed)
    n_varied <- length(nms_varied)
    cov_fixed <- matrix(0, n_fixed, n_all)
    colnames(cov_fixed) <- nms_all
    cov_fixed[, nms_fixed] <- cov(samples$pars[, nms_fixed, 1])

    f <- function(r) {
      cov_varied <- matrix(0, n_varied, n_all)
      colnames(cov_varied) <- nms_all
      cov_varied[, nms_varied] <- cov(samples$pars[, nms_varied, r])
      cov_varied
    }
    cov_varied <- lapply(region, f)

    covariance <- abind_quiet(c(list(cov_fixed), cov_varied), along = 1)

    proposal <- data.frame(
      region = c(rep(NA, n_fixed), rep(region, each = n_varied)),
      name = c(nms_fixed, rep(nms_varied, length(region))),
      covariance,
      stringsAsFactors = FALSE)

  } else {
    i <- which.max(samples$probabilities[, "log_posterior"])
    initial <- samples$pars[i, ]
    info$initial[match(names(initial), info$name)] <- unname(initial)

    covariance <- cov(samples$pars)
    rownames(covariance) <- NULL
    proposal <- data_frame(region = samples$info$region,
                           name = colnames(covariance),
                           covariance)
  }

  parameters$info <- info
  parameters$prior <- prior
  parameters$proposal <- proposal

  parameters
}


calculate_positivity <- function(samples) {
  date <- sircovid::sircovid_date_as_date(samples$trajectories$date)
  state <- samples$trajectories$state
  transform <- samples$predict$transform
  multiregion <- samples$info$multiregion

  if (multiregion) {
    region <- samples$info$region
    pars_model <- lapply(region, function(r)
      lapply(seq_len(nrow(samples$pars)), function(i)
        last(transform[[r]](samples$pars[i, , r]))$pars))
    positivity <- lapply(seq_along(region), function(i)
      calculate_positivity_region(state[, i, , ], pars_model[[i]], date))

    ## Add an extra dimension, then bind together:
    positivity <- lapply(positivity, mcstate::array_reshape,
                         2, c(1, dim(state)[[3]]))
    positivity <- abind_quiet(positivity, along = 2)
  } else {
    pars_model <- lapply(seq_len(nrow(samples$pars)), function(i)
      last(transform(samples$pars[i, ]))$pars)
    positivity <- calculate_positivity_region(state, pars_model, date)
  }

  ## Then bind into the main object:
  samples$trajectories$state <- abind_quiet(state, positivity, along = 1)

  samples
}


calculate_positivity_region <- function(state, pars_model, date) {
  p_NC_names <- c("p_NC_under15", "p_NC_15_24", "p_NC_25_49",
                  "p_NC_50_64", "p_NC_65_79", "p_NC_80_plus",
                  "p_NC_weekend_under15", "p_NC_weekend_15_24",
                  "p_NC_weekend_25_49", "p_NC_weekend_50_64",
                  "p_NC_weekend_65_79", "p_NC_weekend_80_plus")

  pos_under15 <- state[paste0("sympt_cases_under15_inc"), , ]
  pos_15_24 <- state[paste0("sympt_cases_15_24_inc"), , ]
  pos_25_49 <- state[paste0("sympt_cases_25_49_inc"), , ]
  pos_50_64 <- state[paste0("sympt_cases_50_64_inc"), , ]
  pos_65_79 <- state[paste0("sympt_cases_65_79_inc"), , ]
  pos_80_plus <- state[paste0("sympt_cases_80_plus_inc"), , ]

  pos_over25 <- pos_25_49 + pos_50_64 + pos_65_79 + pos_80_plus
  pos_all <- pos_under15 + pos_15_24 + pos_over25

  pars_base <- pars_model[[1]]
  pars <- t(vapply(pars_model, function(p) unlist(p[p_NC_names]),
                   numeric(length(p_NC_names))))

  calc_negs <- function(group) {
    neg <- pars_base[[paste0("N_tot_", group)]] -
      state[paste0("sympt_cases_", group, "_inc"), , ]

    neg[, grepl("^S", weekdays(date))] <-
      neg[, grepl("^S", weekdays(date))] *
      pars[, paste0("p_NC_weekend_", group)]
    neg[, !grepl("^S", weekdays(date))] <-
      neg[, !grepl("^S", weekdays(date))] *
      pars[, paste0("p_NC_", group)]

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
      (pos * pars_base$pillar2_sensitivity +
         neg * (1 - pars_base$pillar2_specificity)) / (pos + neg) * 100
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
  positivity
}


calculate_cases <- function(samples) {
  date <- sircovid::sircovid_date_as_date(samples$trajectories$date)
  state <- samples$trajectories$state
  transform <- samples$predict$transform
  multiregion <- samples$info$multiregion

  if (multiregion) {
    region <- samples$info$region
    pars_model <- lapply(region, function(r)
      lapply(seq_len(nrow(samples$pars)), function(i)
        last(transform[[r]](samples$pars[i, , r]))$pars))
    cases <- lapply(seq_along(region), function(i)
      calculate_cases_region(state[, i, , ], pars_model[[i]], date))

    ## Add an extra dimension, then bind together:
    cases <- lapply(cases, mcstate::array_reshape,
                    2, c(1, dim(state)[[3]]))
    cases <- abind_quiet(cases, along = 2)
  } else {
    pars_model <- lapply(seq_len(nrow(samples$pars)), function(i)
      last(transform(samples$pars[i, ]))$pars)
    cases <- calculate_cases_region(state, pars_model, date)
  }

  ## Then bind into the main object:
  samples$trajectories$state <- abind_quiet(state, cases, along = 1)

  samples
}


calculate_cases_region <- function(state, pars_model, date) {
  cases_names <-
    c("phi_pillar2_cases_under15", "phi_pillar2_cases_15_24",
      "phi_pillar2_cases_25_49", "phi_pillar2_cases_50_64",
      "phi_pillar2_cases_65_79", "phi_pillar2_cases_80_plus",
      "phi_pillar2_cases_weekend_under15", "phi_pillar2_cases_weekend_15_24",
      "phi_pillar2_cases_weekend_25_49", "phi_pillar2_cases_weekend_50_64",
      "phi_pillar2_cases_weekend_65_79", "phi_pillar2_cases_weekend_80_plus")


  pars <- t(vapply(pars_model, function(p) unlist(p[cases_names]),
                   numeric(length(cases_names))))

  calc_cases <- function(group) {
    cases1 <- state[paste0("sympt_cases_", group, "_inc"), , ]
    cases1[, grepl("^S", weekdays(date))] <-
      cases1[, grepl("^S", weekdays(date))] *
      pars[, paste0("phi_pillar2_cases_weekend_", group)]
    cases1[, !grepl("^S", weekdays(date))] <-
      cases1[, !grepl("^S", weekdays(date))] *
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

  cases
}

## adapted from sircovid:::calculate_index
get_names <- function(state_name, suffix_list, suffix0 = NULL) {
  if (is.null(suffix0)) {
    suffixes <- list()
  } else {
    suffixes <- list(suffix0)
  }
  for (i in seq_along(suffix_list)) {
    nm <- names(suffix_list)[[i]]
    if (length(nm) == 0) {
      nm <- ""
    }
    suffixes <- c(suffixes,
                  list(c("", sprintf("_%s%s", nm,
                                     seq_len(suffix_list[[i]] - 1L)))))
  }
  suffixes <- expand.grid(suffixes)
  nms <- apply(suffixes, 1,
               function(x) sprintf("%s%s",
                                   state_name, paste0(x, collapse = "")))
  nms
}
