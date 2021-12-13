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
##' @return A [mcstate::particle_filter] object
##'
##' @export
spim_particle_filter <- function(data, pars, control, deterministic = FALSE) {
  p <- pars$model(pars$initial())
  if (inherits(p, "multistage_parameters")) {
    p <- p[[1]]$pars
  }

  initial_step <- 0 # replaced later
  data <- mcstate::particle_filter_data(data, "date", p$steps_per_day,
                                        initial_step)
  ret <- sircovid::lancelot_particle_filter(data, control$n_particles,
                                            control$n_threads, control$seed,
                                            control$compiled_compare)
  if (deterministic) {
    inputs <- ret$inputs()
    ret <- mcstate::particle_deterministic$new(inputs$data,
                                               inputs$model,
                                               inputs$compare,
                                               inputs$index,
                                               inputs$initial,
                                               inputs$n_threads)
  }

  ret
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
