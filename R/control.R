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
##' @param date_restart Optionally, dates save restart data in the
##'   pmcmc (see [mcstate::pmcmc]
##'
##' @param n_particles number of particles to be used in particle filter if
##'   `short_run = FALSE`
##'
##' @param n_mcmc number of steps to be used in PMCMC if
##'   `short_run = FALSE`
##'
##' @param burnin number of steps out of `n_mcmc` to be used as a burn-in in
##'    PMCMC if `short_run = FALSE`
##'
##' @param workers Logical, indicating if we should enable workers. If
##'   `TRUE`, then a number of workers between 1 and 4 will be used
##'   depending on `n_chains` and the detected number of cores.
##'
##' @return A list of options
##' @export
spim_control <- function(short_run, n_chains, date_restart = NULL,
                         n_particles = 192, n_mcmc = 1500, burnin = 500,
                         workers = TRUE) {
  if (short_run) {
    n_particles <- 10
    n_mcmc <- 20
    n_sample <- 10
    burnin <- 1
  } else {
    n_particles <- n_particles
    n_mcmc <- n_mcmc
    n_sample <- 1000
    burnin <- burnin
  }

  rerun_every <- 100
  if (!is.null(date_restart)) {
    date_restart <- sircovid::sircovid_date(date_restart)
  }

  forecast_days <- 57

  parallel <- spim_control_parallel(n_chains, workers)

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
                                  rerun_random = TRUE)

  particle_filter <- list(n_particles = n_particles,
                          n_threads = parallel$n_threads_total,
                          seed = NULL,
                          compiled_compare = FALSE)

  thin <- ceiling(n_chains * (n_mcmc - burnin) / n_sample)

  forecast <- list(n_sample = n_sample, burnin = burnin,
                   forecast_days = forecast_days,
                   thin = thin)

  list(pmcmc = pmcmc,
       particle_filter = particle_filter,
       forecast = forecast)
}


spim_control_parallel <- function(n_chains, workers) {
  nt <- spim_control_cores()
  if (workers) {
    pos <- 1:4
    nw <- max(pos[nt %% pos == 0 & pos <= n_chains])
  } else {
    nw <- 1L
  }
  message(sprintf("Running on %d workers with %d threads", nw, nt))
  list(n_threads_total = nt, n_workers = nw)
}


spim_control_cores <- function() {
  as.integer(Sys.getenv("CONTEXT_CORES",
                        Sys.getenv("MC_CORES",
                                   getOption("mc.cores", 1))))
}
