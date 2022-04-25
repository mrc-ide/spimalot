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
##' @param multiregion Logical, indicating if we are fitting multiple
##'   regions at once (in which case even the deterministic model may
##'   benefit from multithreading).
##'
##' @param severity Logical, indicating if we are outputting severity
##'   trajectories (e.g. IFR, IHR, HFR). Default to FALSE.
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
##' @param compiled_compare Use a compiled compare function (rather
##'   than the R version).  This can speed things up with the
##'   deterministic models in particular.
##'
##' @param mcmc_path Path to store the mcmc results in
##'
##' @param verbose Logical, indicating if we should print information
##'   about the parallel configuration
##'
##' @return A list of options
##' @export
spim_control <- function(short_run, n_chains, deterministic = FALSE,
                         multiregion = FALSE, severity = FALSE,
                         date_restart = NULL, n_particles = 192, n_mcmc = 1500,
                         burnin = 500, workers = TRUE, n_sample = 1000,
                         n_threads = NULL, compiled_compare = FALSE,
                         mcmc_path = NULL, verbose = TRUE) {
  if (short_run) {
    n_particles <- min(10, n_particles)
    n_mcmc <- min(20, n_mcmc)
    n_sample <- min(10, n_mcmc)
    n_chains <- min(4, n_chains)
    burnin <- 1
  }

  n_steps_retain <- ceiling(n_sample / n_chains)

  rerun_every <- if (deterministic) Inf else 100
  if (!is.null(date_restart)) {
    date_restart <- sircovid::sircovid_date(date_restart)
  }

  parallel <- spim_control_parallel(n_chains, workers, n_threads,
                                    deterministic, multiregion,
                                    verbose)

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
                                  filter_early_exit = !deterministic,
                                  n_burnin = burnin,
                                  n_steps_retain = n_steps_retain,
                                  path = mcmc_path)

  n_threads <- parallel$n_threads_total / parallel$n_workers

  # Force severity to be FALSE if running multiregion fits, as severity
  # weight calculation does not currently work with multiregion = TRUE
  if (multiregion) {
    severity <- FALSE
  }

  particle_filter <- list(n_particles = n_particles,
                          n_threads = n_threads,
                          seed = NULL,
                          compiled_compare = compiled_compare,
                          severity = severity)

  list(pmcmc = pmcmc,
       particle_filter = particle_filter)
}


spim_control_parallel <- function(n_chains, workers, n_threads,
                                  deterministic, multiregion,
                                  verbose = TRUE) {
  n_threads <- n_threads %||% spim_control_cores()
  max_workers <- 4

  if (workers) {
    if (deterministic && !multiregion) {
      ## Increase the number of workers because each will be running
      ## separately
      n_workers <- min(n_chains, n_threads)
      n_threads <- n_workers
    } else if (deterministic && multiregion) {
      ## In the case of the deterministic multiregion case, we get
      ## more from workers (where available) than from within-particle
      ## parallelisation, so let's try and get these onto up to max
      ## workers, then increase the number of threads a bit.
      n_workers <- min(n_chains, n_threads, max_workers)
      n_threads_given <- n_threads
      n_threads <- ceiling(n_threads / n_workers) * n_workers
      if (verbose && n_threads > n_threads_given) {
        message(sprintf("Increasing total threads from %d to %d",
                        n_threads_given, n_threads))
      }
    } else {
      pos <- seq_len(max_workers)
      n_workers <- max(pos[n_threads %% pos == 0 & pos <= n_chains])
    }
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
