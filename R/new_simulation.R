##' Control output created as a result of running simulations. Many of
##' the outputs are large, or require additional computation, so
##' reducing this will make your life more pleasant.
##'
##' @title Create output control for simulations
##'
##' @param keep Character vector of outputs to keep. Must be names of
##'   variables as included in the index of the model run.
##'
##' @param time_series Logical, indicating if a time series should be
##'   output (these are quite large)
##'
##' @param rt Logical, indicating if Rt should be calculated. This is
##'   quite slow.
##'
##' @param rt_weighted Logical, indicating if weighted Rt should be
##'   calculated. This is also quite slow. Only allowed if `rt` is
##'   `TRUE`
##'
##' @param rt_type Character vecrtor of Rt types to use, passed to
##'   [sircovid::lancelot_Rt].  Can be any or all of `eff_Rt_all`,
##'   `eff_Rt_general`, `Rt_all` and `Rt_general`.
##'
##' @param state_by_age Logical, indicating if state by age should be
##'   output. These are really large and *will* cause memory issues
##'   for you if you try to run with long time series unless you have
##'   a lot of RAM.
##'
##' @param vaccination Logical, indicating if vaccination information
##'   should be output.
##'
##' @return An object of `spim_simulate_control_output`. Do not modify
##'   this object after creation.
##'
##' @export
spim_simulate_control_output <- function(keep, time_series = TRUE,
                                         rt = TRUE, rt_weighted = rt,
                                         rt_type = NULL,
                                         state_by_age = TRUE,
                                         vaccination = TRUE) {
  assert_character(keep)
  assert_scalar_logical(time_series)
  assert_scalar_logical(rt)
  assert_scalar_logical(rt_weighted)
  assert_scalar_logical(state_by_age)
  assert_scalar_logical(vaccination)

  if (rt_weighted && !rt) {
    stop("If rt_weighted is TRUE, then rt must be TRUE")
  }

  if (is.null(rt_type)) {
    rt_type <- c("eff_Rt_general", "Rt_general")
  }
  valid_rt <- c("eff_Rt_all", "eff_Rt_general", "Rt_all", "Rt_general")
  err <- setdiff(rt_type, valid_rt)
  if (length(err) > 0) {
    stop("Invalid value for 'rt_type': ", paste(squote(err), collapse = ", "))
  }

  ret <- list(
    keep = keep,
    time_series = time_series,
    rt = rt,
    rt_weighted = rt_weighted,
    rt_type = rt_type,
    ## this likely takes up a lot of the memory so switching off
    state_by_age = state_by_age,
    vaccination = vaccination)
  class(ret) <- c("spim_simulate_control_output", "immutable")
  ret
}


##' Create control parameters to run a set of simulations. This
##' includes parameters that vary compared with the upstream
##' simulation and grids of parameter types to loop over.
##'
##' @title Create control parameters
##'
##' @param flavour The simulation flavour; this is a single string
##'   (e.g., "mtp" or "add_omicron")
##'
##' @param regions A vector of regions to simulate over.  These must
##'   all be valid values in [sircovid::regions]
##'
##' @param start The start date, as an R `Date` object. Typically this
##'   is the last date of the data/fits.
##'
##' @param end The end date of the simulation as an R `Date` object
##'
##' @param parameters A list of parameters.  These might be direct
##'   replacements against the baseline or structured lists of
##'   parameters with names that correspond to those within `grid`.
##'   We'll document this more later!
##'
##' @param grid The parameter grid, indicating the set of simulations
##'   to run. It must have at least one row.
##'
##' @param output Output control, created by
##'   [spimalot::spim_simulate_control]
##'
##' @export
spim_simulate_control <- function(flavour, regions, date_start, date_end,
                                  parameters, grid, output) {
  assert_scalar_character(flavour)

  assert_character(regions)
  err <- setdiff(regions, sircovid::regions("all"))
  if (length(err) > 0) {
    stop("Invalid region: ", paste(squote(err), collapse = ", "))
  }

  assert_is(date_start, "Date")
  assert_is(date_end, "Date")
  if (date_end <= date_start) {
    stop("'date_end' must be greater than 'date_start'")
  }

  assert_is(parameters, "list")
  ## Handling of beta_step is special because it is typically added
  ## later on.
  if (!("beta_step" %in% names(parameters))) {
    parameters["beta_step"] <- list(NULL)
  }

  assert_is(grid, "data.frame")
  if (nrow(grid) == 0) {
    stop("At least one row required in 'grid'")
  }
  if (inherits(grid, "tbl")) {
    grid <- as.data.frame(grid)
  }

  assert_is(output, "spim_simulate_control_output")

  ## TODO: we might make this tuneable when configuring the control
  ## object.
  ignore <- c("analysis", "scenario")
  names_variable <- setdiff(names(grid), ignore)

  ## NOTE: only a simple assertion here, because this is/will be
  ## tested properly above.
  stopifnot(all(names_variable %in% names(parameters)))
  names_constant <- setdiff(names(parameters), names_variable)

  ret <- list(flavour = flavour,
              regions = regions,
              date_start = date_start,
              date_end = date_end,
              parameters = parameters,
              grid = grid,
              output = output)

  validate_simulate_parameters(ret, FALSE)

  ## TODO: this one can't be immutable yet because the simulation task
  ## is in enough of a mess that we need to change these parameters
  ## too much. Eventually we'd like to add that here though, but it
  ## will take some time.
  class(ret) <- "spim_simulate_control"
  ret
}


##' Add `beta_step` into a `spim_simulate_control` object.  Typically
##' this happens just before running simulations, after the Rt values
##' have been converted into beta values according to the assumptions
##' of the simulation.
##'
##' @title Add beta_step into control
##'
##' @param control A control object from [spimalot::spim_simulate_control]
##'
##' @param beta_step Beta values; either a named list (if `beta_step`
##'   appears in the run grid, which it probably will) or a single
##'   set. These should be 3d arrays with dimensions corresponding to
##'   1. particle, 2. region, and 3. step (so time multiplied by
##'   `steps_per_day`)
##'
##' @export
spim_simulate_set_beta_step <- function(control, beta_step) {
  assert_is(control, "spim_simulate_control")
  if (!is.null(control$parameters$beta_step)) {
    stop("'beta_step' has already been set")
  }
  ## TODO: validate that the given values are sensible.
  control$parameters$beta_step <- beta_step

  ## TODO: check that beta_step satisfies the above conditions, but
  ## probably below in validate_simulation_parameters?

  validate_simulate_parameters(control, TRUE)
  control
}


##' Convert an [spimalot::spim_simulate_control] object into a list of
##' parmeters for simulation; these will then need to be swapped into
##' the model too.
##'
##' @title Prepare parameter update list
##'
##' @param control A [spimalot::spim_simulate_control] object
##'
##' @export
spim_simulate_parameter_grid <- function(control) {
  grid <- control$grid
  pars <- control$parameters
  nms <- validate_simulate_parameters(control, TRUE)

  f <- function(i) {
    ret <- pars[nms$constant]
    for (v in nms$variable) {
      level <- grid[[v]][[i]]
      value <- pars[[v]][[level]]
      if (is.null(value)) {
        stop(sprintf("Did not find level '%s' in control$parameters$%s",
                     level, variable))
      }
      ## This function needs to account for the bug:
      ##
      ## When variable == "strain_cross_immunity" we have a list of
      ## three levels equal_to_ve, lower_than_ve, higher_than_ve Each
      ## of these has sub-levels of one_, two_, four_, eight_ and
      ## sixteen_fold
      ##
      ## TODO: Why is this a special snowflake? We're in control of
      ## this, why are we making work for ourselves here?  If for
      ## some reason this turns out to be unavoidable, then a
      ## paragraph or two of explanation might help this be less
      ## obscure.  There's a comment above that I believe relates to
      ## this.  Providing useful feedback here will be hard in this
      ## current form, but we could move this *below*
      if (v == "strain_cross_immunity") {
        level_ve <- grid[["strain_vaccine_efficacy"]][[i]]
        value <- value[[level_ve]]
        if (is.null(value)) {
          stop(sprintf("Did not find level '%s' in control$parameters$%s$%s",
                       level, variable, level_ve))
        }
      }
      ret[[v]] <- value
    }
    ret
  }

  lapply(seq_len(nrow(grid)), f)
}


validate_simulate_parameters <- function(control, require_beta_step) {
  grid <- control$grid
  parameters <- control$parameters

  ## TODO: We'll relax this at some point, and likely need to later
  ## for the mtps
  expected_grid <- c("scenario", "vaccine_daily_doses", "booster_daily_doses",
                     "strain_transmission", "strain_cross_immunity",
                     "strain_vaccine_efficacy", "analysis", "beta_step")
  assert_names_setequal(grid, expected_grid)

  ## TODO: We will want to allow some of these to be missing, so that
  ## this represents some maximal set of things to change.
  expected_parameters <- c("vaccine_eligibility_min_age",
                           "vaccine_booster_eligibility",
                           "vaccine_daily_doses",
                           "booster_daily_doses",
                           "strain_vaccine_efficacy",
                           "vaccine_uptake",
                           "strain_cross_immunity",
                           "strain_transmission",
                           "strain_seed_date",
                           "strain_seed_size",
                           "strain_seed_pattern",
                           "beta_step",
                           "rt_sd",
                           "rt_schools_modifier",
                           "rt_scenarios",
                           "rt_seasonality_date_peak",
                           "rt_seasonality_amplitude")
  assert_names_setequal(parameters, expected_parameters)

  ## TODO: we might make this tuneable when configuring the control
  ## object.
  ignore <- c("analysis", "scenario")
  names_variable <- setdiff(names(grid), ignore)

  ## NOTE: only a simple assertion here, because this is/will be
  ## tested properly above.
  msg <- setdiff(names_variable, names(parameters))
  if (length(msg)) {
    stop(sprintf(
      "All parameters in 'grid' must be in 'parameters', but missing: %s",
      paste(squote(msg), collapse = ", ")))
  }

  if (require_beta_step && is.null(parameters$beta_step)) {
    stop("beta_step has not been added yet")
  }

  names_constant <- setdiff(names(parameters), names_variable)

  list(constant = names_constant, variable = names_variable)
}
