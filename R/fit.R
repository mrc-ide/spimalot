##' Construct a particle filter
##'
##' @title Construct a particle filter
##'
##' @param data Data, from [spimalot::spim_data]
##'
##' @param pars Parameters, from [spimalot::spim_pars]
##'
##' @param control A list of control parameters including
##'   `n_particles`, `n_threads` and `compiled_compare`, typically the
##'   `particle_filter` element of the result of
##'   [spimalot::spim_control()]
##'
##' @param deterministic Logical, indicating if the particle filter to built
##'   is to be run deterministically or stochastically
##'
##' @param initial Optionally a matrix of model state or an initial
##'   value function.
##'
##' @param initial_date The initial date to start from.  Typically
##'   this should be 0 for a parent fit, and be the restart date for a
##'   restart date.
##'
##' @return A [mcstate::particle_filter] object
##'
##' @export
spim_particle_filter <- function(data, pars, control,
                                 deterministic = FALSE,
                                 initial = NULL,
                                 initial_date = 0) {
  ## We do need to get the steps per day out regardless.  A lot of
  ## work considering this is always 4!
  p <- pars$model(pars$initial())
  if (inherits(p, "multistage_parameters")) {
    p <- p[[1]]$pars
  }
  steps_per_day <- p$steps_per_day
  initial_date <- sircovid::as_sircovid_date(initial_date)

  ## First check that all columns are as expected, and that no data
  ## double-use is done
  sircovid::lancelot_check_data(data)

  if (nlevels(data$region) > 1 && length(unique(data$region)) > 1) {
    population <- "region"
  } else {
    population <- NULL
  }

  ## Then organise into mcstate format:
  data <- mcstate::particle_filter_data(data, "date", steps_per_day,
                                        initial_date, population)

  ## We might use the built-in compare function, but probably we will
  ## use the R one from the package.
  compare <- if (control$compiled_compare) NULL else sircovid::lancelot_compare

  if (is.null(initial)) {
    initial <- sircovid::lancelot_initial
  } else if (is.matrix(initial)) {
    ## This likely needs a small tweak in mcstate so that the sampling
    ## works as expected.  Probably not very difficult to get right.
    if (!is.null(population)) {
      stop("restart + multiregion filter will need work")
    }
    initial <- mcstate::particle_filter_initial(initial)
  } else {
    ## TODO: is this branch ever used?  I suspect not now that we
    ## accept a matrix initial, and looks that way from a quick look
    ## through ncov
    assert_is(initial, "function")
  }

  if (deterministic) {
    mcstate::particle_deterministic$new(
      data = data, model = sircovid::lancelot, compare = compare,
      index = sircovid::lancelot_index, initial = initial,
      n_threads = control$n_threads)
  } else {
    mcstate::particle_filter$new(
      data = data, model = sircovid::lancelot,
      n_particles = control$n_particles,
      compare = compare, index = sircovid::lancelot_index, initial = initial,
      n_threads = control$n_threads, seed = control$seed)
  }
}


##' Run the fit (using pmcmc); this may take a very long time
##'
##' @title Run fit
##'
##' @param pars Parameters
##'
##' @param filter A particle filter object, typically from
##'   [spimalot::spim_particle_filter]
##'
##' @param control A [mcstate::pmcmc_control] object, typically the
##'   `pmcmc` element of the result of [spimalot::spim_control()]
##'
##' @return A `mcstate_pmcmc` object
##'
##' @export
spim_fit_run <- function(pars, filter, control) {
  message("Running chains - this will take a while!")
  initial <- replicate(control$n_chains,
                       pars$mcmc$propose(pars$mcmc$initial(), 1))
  ret <- mcstate::pmcmc(pars$mcmc, filter, initial = initial, control = control)

  ## Add some additional version information, which will make the
  ## vaccination projection more robust by preventing us mis-aligning
  ## the updated variables. This will propagate through the forecasts

  ## Data here is an odin parameter
  data <- ret$predict$transform(ret$pars[1, ])
  info <- lapply(data,
                 function(d) ret$predict$filter$model$new(d$pars, 0, 1)$info())

  ## Just to make the below a bit less cumbersome.
  base <- pars$base

  ret$info <- list(version = packageVersion("sircovid"),
                   info = info,
                   data = data,
                   ## NOTE: probably we'll change format here soon,
                   ## but these are the "small" inputs to spim_pars
                   date = base$date,
                   region = base$region,
                   beta_date = base$beta_date,
                   epoch_dates = base$epoch_dates,
                   model_type = base$model_type,
                   restart_date = base$restart_date)

  ret
}
