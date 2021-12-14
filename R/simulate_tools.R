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
    vaccine_booster_efficacy = "central",
    waning_rate = "central", strain_cross_immunity = "central"
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

##' Aggregate simulation output by region to produce national estimates
##'
##' @title Combine regional simulated trajectories
##' @param res A list of simulated objects
##'
##' @param name The name of the aggregated region (e.g. 'england)
##'
##' @param regions Character vector of regions to use
##'
##' @param rm.rtUK Logical, indicating if the aggregated UK trajectory should be
##' removed
##'
##' @export
spim_simulate_combine_trajectories <- function(res, name, regions = NULL,
                                               rm.rtUK = FALSE) {

  # return immediately for summary run
  if (!("state" %in% names(res))) {
    return(res)
  }

  state <- res$state
  dim_state <- dim(state)
  dim_region <- 3L

  regions <- regions %||% all_regions
  rt_names <- intersect(grep("Rt", names(res), value = TRUE),
                        grep("multi", names(res), invert = TRUE, value = TRUE))

  dim_combined <- dim_state
  combined_state <- state[, , 1, , drop = TRUE]
  combined_state[] <- 0

  for (r in regions) {
    for (k in seq_len(dim_state[1])) {
      combined_state[k, , ] <- combined_state[k, , ] + state[k, , r, ]
    }
  }


  combined_Rt <- NULL

  if (length(rt_names) > 0) {

    if (rm.rtUK) {
      for (rt_type in rt_names) {
        which <- match("uk", colnames(res[[rt_type]]))
        if (!is.na(which)) {
          res[[rt_type]] <- res[[rt_type]][, -which, , ]
        }
      }
    }

    which <- vapply(
      rt_names, function(x) all(name %in% colnames(res[[x]])),
      logical(1)
    )

    if (!(sum(which) %in% c(0, length(rt_names)))) {
      stop(sprintf("Region %s is in some of 'rt_names' but not all", name))
    } else if (sum(which) == 0) {
      incid <- apply(state["infections", , , ], c(1, 2), diff)
      combined_Rt <- setNames(vector("list", length(rt_names)), rt_names)
      combined_Rt[] <- 0
      sum_w <- combined_Rt

      for (rt_type in rt_names) {

        # add strain dimension to rt if needed
        dim_rt <- dim(res[[rt_type]])
        if (length(dim_rt) == 3) {
          dim(res[[rt_type]]) <- c(dim_rt, 1)
          dimnames(res[[rt_type]])[[4]] <- "both"
        }
        for (r in regions) {

          x <- res[[rt_type]][, r, , , drop = FALSE]
          w <- c(cbind(NA, t(incid[, , r])))
          combined_Rt[[rt_type]] <- combined_Rt[[rt_type]] + x * w
          sum_w[[rt_type]] <- sum_w[[rt_type]] + w

        }
      }

      combined_Rt <- Map("/", combined_Rt, sum_w)

      if (sum(incid[, , regions]) != sum(sum_w[[1]], na.rm = TRUE)) {
        stop("error in weighting calc")
      }

      if (sum(apply(incid[, , regions], c(1, 2), sum)) !=
          sum(sum_w[[1]], na.rm = TRUE)) {
        stop("error in weighting calc")
      }
    }
  }

  bind_arrays <- function(main, combined, dim_region, name) {
    all_regions <- dimnames(main)[[dim_region]]
    x <- abind_quiet(main, combined, along = dim_region)
    dimnames(x)[[dim_region]] <- c(all_regions, name)
    x
  }

  combine_n_protected <- function(x, regions, name) {
    combined_n_protected <- apply(x[, regions, ], c(1, 3), sum)
    bind_arrays(x, combined_n_protected, 2L, name)
  }

  if (is.list(res$n_protected)) {
    n_protected <- lapply(res$n_protected, combine_n_protected, regions, name)
  } else {
    n_protected <-
      list(strain_1 = combine_n_protected(res$n_protected, regions, name))
  }

  agg_regions <- function(x) {
    i_region <- 3
    dims <- seq_along(dim(x))
    sum_x <- apply(x[, , regions, ], dims[-i_region], sum)
    bind_arrays(x, sum_x, i_region, name)
  }

  ret <- res
  ret$state <- bind_arrays(state, combined_state, dim_region, name)
  ret$n_protected <- n_protected
  ret$n_doses <- agg_regions(res$n_doses)
  if (!is.null(combined_Rt)) {
    ret[rt_names] <- Map(bind_arrays, res[rt_names], combined_Rt, 2, name)
  }

  if ("state_by_age" %in% names(res)) {
    ret$state_by_age <- lapply(res$state_by_age, agg_regions)
  }

  if ("n_vaccinated" %in% names(res)) {
    ret$n_vaccinated <- agg_regions(res$n_vaccinated)
  }

  ret
}


##' Simplify simulation object by moving Rt in with trajectories
##' @title  Simplify simulation object by moving Rt in with trajectories
##' @param x A simulated object

##' @export
spim_simulate_simplify_rt <- function(x) {
  # return immediately for summary run
  if (!("state" %in% names(x))) return(x)

  names_rt <- c("Rt_general", "eff_Rt_general")
  for (v in names_rt) {
    if (length(dim(x[[v]])) == 3) {
      x[[v]] <- array(x[[v]], c(1, dim(x[[v]])))
    } else if (length(dim(x[[v]])) == 4) {
      ## move strain to front
      x[[v]] <- aperm(x[[v]], c(4, 1, 2, 3))
      rownames(x[[v]]) <- paste(v, rownames(x[[v]]), sep = "_")
    }
  }

  rt <- abind_quiet(x[names_rt], along = 1L)
  x$state <- abind_quiet(x$state, rt, along = 1L)


  x[setdiff(names(x), names_rt)]
}




##' Combine diagnoses and admitted trajectories in simulated object
##' @title  Combine diagnoses and admitted trajectories
##' @param obj A simulated object
##'
##' @param incidence Logical, whether you are combining incidence trajectories
##'    or not
##'
##' @export
spim_simulate_add_diagnoses_admitted <- function(obj, incidence = FALSE) {

  # return immediately for summary run
  if (!("state" %in% names(obj))) return(obj)

  nms <- c("diagnoses", "admitted")
  if (incidence) {
    nms <- paste0(nms, "_inc")
  }

  if (all(nms %in% rownames(obj$state))) {
    state <- apply(obj$state[nms, , , , drop = FALSE], c(2, 3, 4), sum)
    dim(state) <- c(1, dim(state))
    if (incidence) {
      rownames(state) <- "diagnoses_admitted_inc"
    } else {
      rownames(state) <- "diagnoses_admitted"
    }

    obj$state <- abind_quiet(obj$state, state, along = 1L)
  }

  obj
}


##' Combine all death trajectories in simulated object
##' @title  Combine all death trajectories
##' @param obj A simulated object
##'
##' @param incidence Logical, whether you are combining incidence trajectories
##'    or not
##'
##' @export
spim_simulate_add_all_deaths <- function(obj, incidence = FALSE) {

  # return immediately for summary run
  if (!("state" %in% names(obj))) return(obj)

  nms <- c("deaths_hosp", "deaths_carehomes", "deaths_comm")
  if (incidence) {
    nms <- paste0(nms, "_inc")
  }

  if (all(nms %in% rownames(obj$state))) {
    state <- apply(obj$state[nms, , , , drop = FALSE], c(2, 3, 4), sum)
    dim(state) <- c(1, dim(state))
    if (incidence) {
      rownames(state) <- "deaths_inc"
    } else {
      rownames(state) <- "deaths"
    }

    obj$state <- abind_quiet(obj$state, state, along = 1L)
  }

  obj
}


##' Removes all results from a simulation object up to and including a
##'  given date
##' @title  Remove simulations up to given date
##' @param obj A simulated object
##' @param date If not `NULL` then the date to remove all results up
##'  to
spim_simulate_remove_dates_to <- function(obj, date) {
  if (is.null(date)) {
    return(obj)
  }

  if (min(obj$date) == sircovid::sircovid_date(date)) {
    # For most simulations, we include state variables from the date of the
    # simulation start onward
    remove_to_date <- date

  } else if (min(obj$date) < sircovid::sircovid_date(date)) {
    # e.g. for an MTP that is run off fits from Friday, we would remove
    # state variables prior to the Monday when submission is due
    remove_to_date <- date - 1

  } else {
    stop(message("Error, obj$date cannot be greater than date"))
  }

  id0 <- seq(which(sircovid::sircovid_date_as_date(obj$date) == remove_to_date))
  obj$date <- obj$date[-id0]
  obj$state <- obj$state[, , , -id0, drop = FALSE]
  obj$state_by_age <- lapply(obj$state_by_age, function(x)
    x[, , , -id0, drop = FALSE])
  obj$n_vaccinated <- obj$n_vaccinated[, , , -id0, drop = FALSE]
  obj$n_protected <- lapply(obj$n_protected, function(x)
    x[, , -id0, drop = FALSE])
  obj$n_doses <- obj$n_doses[, , , -id0, drop = FALSE]
  obj
}


##' Calculate incidence trajectories for a simulated object
##' @title Calculate incidence trajectories for a simulated object
##' @param obj A simulated object
##' @param states Character vector of states to calculate incidence for
##' @param suffix Character giving suffix to add to names of new states, default
##' "_inc"
##' @export
spim_simulate_add_trajectory_incidence <- function(obj, states,
                                                   suffix = "_inc") {

  # return immediately for summary run
  if (!("state" %in% names(obj))) return(obj)

  # TODO:: move into sircovid (this is an updated version of
  #  sircovid::add_trajecctory_incidence)
  if (length(states) == 0) {
    return(obj)
  }

  calc_incidence <- function(x) {
    res <- apply(x, c(1, 2, 3), function(x) c(NA, NA, diff(x[-1L])))
    aperm(res, c(2, 3, 4, 1))
  }

  traj_inc <- calc_incidence(obj$state[states, , , , drop = FALSE])
  rownames(traj_inc) <- paste0(states, suffix)
  obj$state <- abind_quiet(obj$state, traj_inc, along = 1L)

  if ("state_by_age" %in% names(obj)) {

    state_by_age_nms <- intersect(names(obj$state_by_age), states)

    state_by_age_inc <- lapply(obj$state_by_age[state_by_age_nms],
                               calc_incidence)
    names(state_by_age_inc) <- paste0(state_by_age_nms, suffix)

    obj$state_by_age <- c(obj$state_by_age, state_by_age_inc)

  }


  obj
}


##' @title Reset cumulative states to zero at start of simulation
##' @param res A simulated object
##' @param state_names Character vector of states to reset
##' @export
spim_simulate_reset_cumulative_states <- function(res, state_names) {
  # return immediately for summary run
  if (!("state" %in% names(res))) return(res)

  f <- function(x) x - c(x[, , , 1])

  res$state[state_names, , , ] <- f(res$state[state_names, , , , drop = FALSE])

  check <- all(res$state[state_names, , , 3] - res$state[state_names, , , 2] ==
                 res$state[paste0(state_names, "_inc"), , , 3])
  stopifnot(check)

  if ("state_by_age" %in% names(res)) {
    state_by_age_nms <- intersect(names(res$state_by_age), state_names)
    res$state_by_age[state_by_age_nms] <-
      lapply(res$state_by_age[state_by_age_nms], f)
  }

  res
}


##' @title Process simulation output
##' @param obj A simulated object
##' @param combined_region The name of the aggregated region (e.g. 'england)
##' @param regions Character vector of regions to use
##' @param incidence_states  Character vector of states to calculate incidence
##' @param reset_states Logical, indicating if incidence states should be reset
##' to zero at the start of the simulation
##' @param rm.rtUK Logical, indicating if the aggregated UK trajectory should be
##' removed
##' @param output_region Character vector of regions to output, defaults to
##' `combined_region`
##' @param simulation_start_date If not `NULL` then removes all data before
##'  the given date.
##' @export
spim_simulate_process_output <- function(obj, combined_region, regions,
                                         incidence_states,
                                         reset_states = FALSE,
                                         rm.rtUK = FALSE,
                                         output_region = NULL,
                                         simulation_start_date = NULL) {

  # IF output_region is NULL, outputs for all regions and national aggregation
  # will be processed. Else, only the output_region's results will (e.g.
  # aggregation to the England or UK level)
  output_region_only <- !is.null(output_region)
  output_region <- output_region %||% combined_region
  ret <- spim_simulate_combine_trajectories(obj, combined_region, regions,
                                            rm.rtUK)
  ret <- spim_simulate_simplify_rt(ret)
  ret <- spim_simulate_add_diagnoses_admitted(ret)
  ret <- spim_simulate_add_trajectory_incidence(ret, incidence_states)
  if (reset_states) {
    ret <- spim_simulate_reset_cumulative_states(ret, incidence_states)
  }
  ret$multivariant_Rt_general <- NULL
  ret$multivariant_eff_Rt_general <- NULL

  if (output_region_only) {
    f <- function(x) x[, , output_region, , drop = FALSE]

    ret$summary_state <- f(ret$summary_state)
    ret$state <- f(ret$state)
    ret$state_by_age <- lapply(ret$state_by_age, f)
    ret$n_vaccinated <- f(ret$n_vaccinated)
    ret$n_protected <- lapply(ret$n_protected,
                              function(x) x[, output_region, , drop = FALSE])
    ret$n_doses <- f(ret$n_doses)
  }

  ret <- spim_simulate_remove_dates_to(ret, simulation_start_date)

  ret
}


##' @title Summarise simulations over particles
##' @param x A simulated object
##' @param at Numeric vector of quantiles at which to calculate summary,
##' defaults to 2.5, 50, 97.5
##' @export
spim_simulate_create_summary <- function(x, at = c(1, 20, 39) / 40) {
  f <- function(x) quantile_digest(x, at)

  if ("summary_state" %in% names(x)) {
    summary_state <- aperm(apply(x$summary_state, c(1, 3, 4), f), c(2, 1, 3, 4))
    colnames(summary_state) <- paste0(at * 100, "%")
    x$summary_state <- round(summary_state)
  }

  # return immediately for summary run
  if (!("state" %in% names(x))) return(x)

  state <- aperm(apply(x$state, c(1, 3, 4), f), c(2, 1, 3, 4))
  colnames(state) <- paste0(at * 100, "%")

  # round everything except Rt and prop_strain_2
  nms <- rownames(state)
  not_to_round <- c(grep("Rt_", nms),
                    grep("prop", nms))
  if (any(not_to_round)) {
    nms <- nms[-not_to_round]
  }
  state[nms, , , ] <- round(state[nms, , , ])

  x$state <- state

  x
}


## This is great, but there's a garbage protection bug (possibly in
## Rcpp) that is causing the tdigest object to get get collected. In
## this case we fall back on quantile.

quantile_digest <- function(x, at) {
  if (anyNA(x)) {
    x <- x[!is.na(x)]
  }
  tryCatch(
    tdigest::tquantile(tdigest::tdigest(x), at),
    error = function(e) quantile(x, at, names = FALSE))
}


##' @title Create tidy (long) dataframe of simulated results
##' @param res A list of simulated objects
##' @param run_grid data.frame with same number of rows as `res` giving metadata
##' to attach to results
##' @param combined Ignored
##' @export
spim_simulate_tidy_states <- function(res, run_grid, combined) {
  stopifnot(length(res) == nrow(run_grid))
  ret <- lapply(seq_along(res), function(i)
    tidy_state_one(res[[i]], run_grid[i, ]))

  lapply(list_transpose(ret), dplyr::bind_rows)
}

##' @title Create tidy (long) dataframe of simulated results
##' @param x A simulated object
##' @param common A one-row data.frame giving metadata to attach to results
##' @export
tidy_state_one <- function(x, common) {
  stopifnot(nrow(common) == 1L)
  res <- list()
  if ("summary_state" %in% names(x)) {
    ## This tests if the incoming object is already summarised to
    ## quantiles (we can tell based on the names that would have been
    ## added to the second dimension)
    if (is.null(colnames(x$summary_state))) {
      name_2 <- "particle"
    } else {
      name_2 <- "quantile"
    }
    s <- aperm(x$summary_state, c(4, 2, 3, 1))
    dn <- append(dimnames(s), list("all", "all"), 2)
    names(dn) <- c("date", name_2, "group", "vaccine_status", "region", "state")
    dn$date <- sircovid::sircovid_date_as_date(max(x$date))
    if (is.null(dn[[2]])) {
      dn[[2]] <- seq_len(ncol(s))
    }
    ret <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret$value <- c(s)
    res$summary_state <- ret

  }

  if ("state" %in% names(x)) {
    ## As above, do we have an already-summarised object?
    if (is.null(colnames(x$state))) {
      name_2 <- "particle"
    } else {
      name_2 <- "quantile"
    }
    ## fastest: date, particle, group = all, vaccine = all, region, state

    s <- aperm(x$state, c(4, 2, 3, 1))
    dn <- append(dimnames(s), list("all", "all"), 2)
    names(dn) <- c("date", name_2, "group", "vaccine_status", "region", "state")
    dn$date <- sircovid::sircovid_date_as_date(x$date)
    if (is.null(dn[[2]])) {
      dn[[2]] <- seq_len(ncol(s))
    }
    ret <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret$value <- c(s)

    # date, region, strain, state
    p <- aperm(abind_quiet(x$n_protected, along = 4), c(3, 2, 4, 1))
    dnp <- set_names(dimnames(p), c("date", "region", "strain", "state"))
    dnp$date <- sircovid::sircovid_date_as_date(x$date)
    dnp$quantile <- "mean"

    ret_p <- do.call(expand.grid, dnp)[seq(length(dnp), 1)]
    ret_p$value <- c(p)

    # date, [particle] = mean, group, region, state
    d <- aperm(x$n_doses, c(4, 1, 3, 2))
    dn$group <- c(sircovid:::sircovid_age_bins()$start, "CHW", "CHR")
    dn$state <- colnames(x$n_doses)
    dn$quantile <- "mean"

    ret_d <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret_d$value <- c(d)


    # date, [particle] = mean, group, vaccine_status, region, state
    # 'av' is "age and vaccine"
    av <-
      unlist(lapply(x$state_by_age, aperm, c(4, 1, 2, 3)), use.names = FALSE)
    dn$group <- rownames(x$state_by_age[[1]])
    dn$vaccine_status <- colnames(x$state_by_age[[1]])
    dn$state <- names(x$state_by_age)
    ret_av <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret_av$value <- av

    res$state <- ret
    res$n_protected <- ret_p
    res$n_doses <- ret_d
    res$state_by_age <- ret_av

    rownames(common) <- NULL
  }

  lapply(res, function(x) cbind(common, x))
}
