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
                                         state_by_age = FALSE,
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
##' @param date_start The start date, as an R `Date` object. Typically this
##'   is the last date of the data/fits.
##'
##' @param date_end The end date of the simulation as an R `Date` object
##'
##' @param expected A character vector of expected parameter names.
##'   These are the names of things that *must* be present within
##'   `parameters`, and this will be checked for you.  It will be
##'   tempting to write `names(parameters)` here, but we suggest not
##'   doing that as this will provide a human readable description of
##'   what you are currently looking to simulate against.
##'
##' @param parameters A list of parameters.  These might be direct
##'   replacements against the baseline or structured lists of
##'   parameters with names that correspond to those within `grid`.
##'   We'll document this more later!
##'
##' @param grid The parameter grid, indicating the set of simulations
##'   to run. It must have at least one row.  Not all variables are
##'   allowed in here; in particular, we disallow most of the `rt_`
##'   variables.  An error will be thrown if anything unexpected is
##'   found.
##'
##' @param output Output control, created by
##'   [spimalot::spim_simulate_control]
##'
##' @export
spim_simulate_control <- function(flavour, regions, date_start, date_end,
                                  expected, parameters, grid, output) {
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

  assert_character(expected)

  assert_is(parameters, "list")
  ## Handling of beta_step is special because it is typically added
  ## later on.
  if (!("beta_step" %in% names(parameters))) {
    parameters["beta_step"] <- list(NULL)
  }

  if (!is.null(parameters$rt_schools_schedule)) {
    parameters$rt_schools_schedule <-
      validate_rt_schools_schedule(parameters$rt_schools_schedule,
                                   regions)
  }

  assert_is(grid, "data.frame")
  if (nrow(grid) == 0) {
    stop("At least one row required in 'grid'")
  }
  if (inherits(grid, "tbl")) {
    grid <- as.data.frame(grid)
  }

  assert_is(output, "spim_simulate_control_output")

  ret <- list(flavour = flavour,
              regions = regions,
              date_start = date_start,
              date_end = date_end,
              expected = expected,
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
                     level, v))
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
                       level, v, level_ve))
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
  expected <- control$expected

  ## These are things that we "know" how to modify.  Most of that
  ## modification still happens in the simulate task but will move
  ## (back) into spimalot fairly shortly.  The provided set of
  ## parameters must be a subset of these.
  allowed_parameters <- c("vaccine_eligibility_min_age",
                          "booster_eligibility_min_age",
                          "vaccine_daily_doses",
                          "booster_daily_doses",
                          "strain_vaccine_efficacy",
                          "vaccine_waning_days",
                          "vaccine_uptake",
                          "serial_interval_scale",
                          "strain_cross_immunity",
                          "strain_transmission",
                          "strain_seed_date",
                          "strain_seed_size",
                          "strain_seed_pattern",
                          "strain_seed_breakdown",
                          "beta_step",
                          "rt_sd",
                          "rt_schools_schedule",
                          "rt_schools_modifier",
                          "rt_scenarios",
                          "rt_seasonality_date_peak",
                          "rt_seasonality_amplitude")

  required_parameters <- "beta_step"

  ## These must not vary across the grid and are constrained to being
  ## shared.  The rt_scenarios one is special as it implies the
  ## scenarios that become beta_step eventually.
  allowed_grid <- setdiff(
    allowed_parameters,
    c("rt_sd", "rt_schools_modifier", "rt_schools_schedule", "rt_scenarios",
      "rt_seasonality_date_peak", "rt_seasonality_amplitude"))

  err <- setdiff(expected, allowed_parameters)
  if (length(err) > 0) {
    stop("Don't know how to work with parameter: ",
         paste(squote(err), collapse = ", "))
  }
  msg <- setdiff(required_parameters, expected)
  if (length(msg) > 0) {
    stop("Your expected list must include required parameter: ",
         paste(squote(msg), collapse = ", "))
  }

  assert_names_setequal(parameters, expected)

  ## TODO: we might make this tuneable when configuring the control
  ## object.
  ignore <- c("analysis", "scenario")
  names_variable <- setdiff(names(grid), ignore)

  err <- setdiff(names_variable, expected)
  if (length(err) > 0) {
    stop("Unexpected parameter in grid (not found in 'parameters'): ",
         paste(squote(err), collapse = ", "))
  }
  err <- setdiff(names_variable, allowed_grid)
  if (length(err) > 0) {
    stop("Disallowed parameter in grid (must be constant across simulations): ",
         paste(squote(err), collapse = ", "))
  }

  err <- setdiff(names_variable, expected)
  if (length(err) > 0) {
    stop("Unknown parameter in grid not found in parameters: ",
         paste(squote(err), collapse = ", "))
  }

  if (require_beta_step && is.null(parameters$beta_step)) {
    stop("beta_step has not been added yet")
  }

  names_constant <- setdiff(names(parameters), names_variable)

  ## TODO: should we check the levels here? Seems sensible.

  list(constant = names_constant, variable = names_variable)
}


validate_rt_schools_schedule <- function(x, regions,
                                         name = deparse(substitute(x))) {
  ## Further checks are of course possible, especially:
  ## * Every closure is followed by an opening
  assert_is(x, "data.frame", name = name)
  assert_has_names(x, c("nation", "year", "month", "day", "schools"),
                   name = name)

  ## Simplify date handling onwards
  x$date <- as.Date(sprintf("%s-%s-%s", x$year, x$month, x$day))
  x$year <- NULL
  x$month <- NULL
  x$day <- NULL

  ## Expand national to regional closures
  stopifnot(is.null(x$region))
  x$region <- x$nation
  i <- which(x$nation == "england")
  regions <- sircovid::regions("england")
  x_eng <- x[rep(i, length(regions)), ]
  x_eng$region <- rep(regions, each = length(i))
  ret <- rbind(x[-i, ], x_eng)

  ## Check we have all wanted regions:
  msg <- setdiff(regions, ret$region)
  if (length(msg) > 0) {
    stop("Missing school closure information for region: ",
         paste(squote(msg), collapse = ", "))
  }

  ## Filter to regions we will simulate with:
  ret <- ret[ret$region %in% regions, ]
  rownames(ret) <- NULL

  ret
}


##' Calculate multiplicative beta scaling factor for schools.  We
##' assume here that schools has come from
##' [spimalot::spim_simulate_control] but we may add validation later.
##'
##' @title Calculate multiplicative beta for schools
##'
##' @param dates A sircovid date, fractional dates (e.g. 100.25) are
##'   allowed.
##'
##' @param schedule A data.frame of school open/closed
##'   information. Required columns are `region`, `date` and
##'   `schools`; really this should have been sorted out by
##'   [spimalot::spim_simulate_control]
##'
##' @param modifier The amount to *decrease* beta by when schools are
##'   closed for holidays etc. A value of 0.15 is a 15% reduction in
##'   contacts, or a multiplicative beta of 0.85
##'
##' @param region The region to subset from the schedule.
##'
##' @return A vector the same length as `dates` with the scaling
##'   factor.
##'
##' @export
spim_beta_mult_schools <- function(dates, schedule, modifier, region) {
  schedule <- schedule[schedule$region == region, ]
  idx <- findInterval(dates, sircovid::sircovid_date(schedule$date))
  i <- schedule$schools[idx] == "closed"
  ret <- rep(1, length(dates))
  ret[i] <- 1 - modifier
  ret
}


##' Calculate multiplicative beta scaling factor due to seasonality.
##'
##' @title Calculate multiplicative beta due to seasonality
##'
##' @inheritParams spim_beta_mult_schools
##'
##' @param date_peak Date of peak as a
##'   [sircovid::sircovid_date] value (a single number)
##'
##' @param value Size of the peak to annual average difference. A
##'   value of 0.1 will produce multiplicative values that range from
##'   0.9 to 1.1.
##'
##' @return A vector the same length as `dates` with the scaling
##'   factor.
##'
##' @export
spim_beta_mult_seasonality <- function(dates, date_peak, value) {
  delta <- ((dates - date_peak) %% 365) / 365
  1 + cos(2 * pi * delta) * value
}
