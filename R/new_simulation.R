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
  assert_is(grid, "data.frame")
  if (nrow(grid) == 0) {
    stop("At least one row required in 'grid'")
  }
  assert_is(output, "spim_simulate_control_output")

  ret <- list(flavour = flavour,
              regions = regions,
              date_start = date_start,
              date_end = date_end,
              parameters = parameters,
              grid = grid,
              output = output)
  ## TODO: this one can't be immutable yet because the simulation task
  ## is in enough of a mess that we need to change these parameters
  ## too much. Eventually we'd like to add that here though, but it
  ## will take some time.
  class(ret) <- "spim_simulate_control"
  ret
}
