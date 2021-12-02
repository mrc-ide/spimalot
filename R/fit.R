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
##' @param model The name of the sircovid model used. Default is `"carehomes"`
##'
##' @param deterministic Logical, indicating if the particle filter to built
##'   is to be run deterministically or stochastically
##'
##' @return A [mcstate::particle_filter] object
##'
##' @export
spim_particle_filter <- function(data, pars, control, model = "carehomes",
                                 deterministic = FALSE) {
  p <- pars$model(pars$initial())
  if (inherits(p, "multistage_parameters")) {
    p <- p[[1]]$pars
  }

  initial_step <- 0 # replaced later
  data <- mcstate::particle_filter_data(data, "date", p$steps_per_day,
                                        initial_step)
  if (model == "lancelot") {
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
  } else {
    stop(sprintf("Unknown model '%s'", model))
  }
  ret
}


##' Run a pmcmc; this may take a very long time
##'
##' @title Run a pmcmc
##'
##' @param pars Parameters, from [spimalot::spim_pars]
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
spim_pmcmc <- function(pars, filter, control) {
  message("Running chains - this will take a while!")

  ## Each chain starts from the initial() result plus an initial proposal step
  ## Note that it is similar to one MCMC step with p(acceptance) = 1
  initial <- replicate(control$n_chains, pars$propose(pars$initial(), 1))
  ret <- mcstate::pmcmc(pars, filter, initial = initial, control = control)

  ## Later on, we'll need to access a number of inputs
  pars_inputs <- attr(pars, "inputs")

  ## Add some additional version information, which will make the
  ## vaccination projection more robust by preventing us mis-aligning
  ## the updated variables. This will propagate through the forecasts
  data <- ret$predict$transform(ret$pars[1, ])
  info <- ret$predict$filter$model$new(data, 0, 1)$info()
  ret$info <- list(version = packageVersion("sircovid"),
                   info = info,
                   data = data,
                   ## NOTE: probably we'll change format here soon,
                   ## but these are the "small" inputs to spim_pars
                   date = pars_inputs$date,
                   region = pars_inputs$region,
                   beta_date = pars_inputs$beta_date,
                   multistrain = pars_inputs$multistrain,
                   model_type = pars_inputs$model_type)

  ## Copy across the core vaccination data that will be used in any
  ## restart/forward simulation later.

  vaccination_copy <- c("schedule_real", "priority_population",
                        "mean_days_between_doses", "efficacy")
  ret$vaccine <- pars_inputs$vaccination[vaccination_copy]

  ret
}
