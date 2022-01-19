##' High-level control.
##'
##' TODO: document parallel detection
##'
##' @title High-level control
##'
##' @param short_run Logical, indicating if this is a debug run. This
##'   will use freer particles, mcmc steps and a small sample
##'
##' @param n_chains The number of chains to run
##'
##' @param deterministic Logical, indicating if the model to fit to
##'   data is run deterministically or stochastically
##'
##' @param date_restart Optionally, dates save restart data in the
##'   pmcmc (see [mcstate::pmcmc]
##'
##' @param n_particles number of particles to be used in particle filter if
##'   `short_run = FALSE`
##'
##' @param n_mcmc number of steps to be used in PMCMC if
##'   `short_run = FALSE`
##'
##' @param n_sample number of steps to retain (across all chains) if
##'   `short_run = FALSE`
##'
##' @param burnin number of steps out of `n_mcmc` to be used as a burn-in in
##'    PMCMC if `short_run = FALSE`
##'
##' @param workers Logical, indicating if we should enable workers. If
##'   `TRUE`, then a number of workers between 1 and 4 will be used
##'   depending on `n_chains` and the detected number of cores.
##'
##' @param n_threads Explicit number of threads, overriding detection
##'   by [spim_control_cores]
##'
##' @return A list of options
##' @export
spim_control <- function(short_run, n_chains, deterministic = FALSE,
                         date_restart = NULL,
                         n_particles = 192, n_mcmc = 1500, burnin = 500,
                         forecast_days = 57, workers = TRUE, n_sample = 1000,
                         n_threads = NULL, rt = FALSE, cum_admit = FALSE,
                         diagnoses_admitted = FALSE,
                         cum_n_vaccinated = FALSE,
                         cum_infections_disag = FALSE) {

  if (short_run) {
    n_particles <- min(10, n_particles)
    n_mcmc <- min(20, n_mcmc)
    n_sample <- 10
    burnin <- 1
  }

  n_steps_retain <- ceiling(n_sample / n_chains)

  rerun_every <- 100
  if (!is.null(date_restart)) {
    date_restart <- sircovid::sircovid_date(date_restart)
  }

  parallel <- spim_control_parallel(n_chains, workers, n_threads)

  pmcmc <- mcstate::pmcmc_control(n_mcmc,
                                  n_chains = n_chains,
                                  n_threads_total = parallel$n_threads_total,
                                  n_workers = parallel$n_workers,
                                  save_state = TRUE,
                                  save_trajectories = TRUE,
                                  use_parallel_seed = TRUE,
                                  progress = interactive(),
                                  save_restart = date_restart,
                                  nested_step_ratio = 1, # ignored if single
                                  rerun_every = rerun_every,
                                  rerun_random = TRUE,
                                  filter_early_exit = TRUE,
                                  n_burnin = burnin,
                                  n_steps_retain = n_steps_retain)

  if (deterministic) {
    ## Disable early exit, if it's been set up, as we also don't support that
    pmcmc$filter_early_exit <- FALSE

    ## Increase the number of workers because each will be running
    ##   separately. If running on a laptop this probably does not want
    ##   increasing
    if (workers) {
      pmcmc$n_workers <- min(pmcmc$n_chains, pmcmc$n_threads_total)
    }
  }

  particle_filter <- list(n_particles = n_particles,
                          n_threads = parallel$n_threads_total,
                          seed = NULL,
                          compiled_compare = FALSE,
                          rt = rt,
                          cum_admit = cum_admit,
                          diagnoses_admitted = diagnoses_admitted,
                          cum_n_vaccinated = cum_n_vaccinated,
                          cum_infections_disag = cum_infections_disag)

  list(pmcmc = pmcmc,
       particle_filter = particle_filter)
}


## This is only going to be really useful with the single region fit,
## as we don't really use workers otherwise.
spim_control_parallel <- function(n_chains, workers, n_threads = NULL,
                                  verbose = TRUE) {
  n_threads <- n_threads %||% spim_control_cores()
  if (workers) {
    pos <- 1:4
    n_workers <- max(pos[n_threads %% pos == 0 & pos <= n_chains])
  } else {
    n_workers <- 1L
  }
  if (verbose) {
    message(sprintf("Running on %d workers with %d threads",
                    n_workers, n_threads))
  }
  list(n_threads_total = n_threads, n_workers = n_workers)
}


##' Return the number of available cores, looking at `CONTEXT_CORES`,
##' `MC_CORES` and `mc.cores` in turn
##'
##' @title Return the number of cores
##' @return An integer
##' @export
spim_control_cores <- function() {
  as.integer(Sys.getenv("CONTEXT_CORES",
                        Sys.getenv("MC_CORES",
                                   getOption("mc.cores", 1))))
}
