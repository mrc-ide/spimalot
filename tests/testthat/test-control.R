test_that("Detect cores", {
  withr::with_envvar(
    c(CONTEXT_CORES = NA_character_, MC_CORES = NA_character_),
    withr::with_options(list(mc.cores = NULL), {
      expect_equal(spim_control_cores(), 1)

      options(mc.cores = 5)
      expect_equal(spim_control_cores(), 5)

      Sys.setenv(MC_CORES = 4)
      expect_equal(spim_control_cores(), 4)

      Sys.setenv(CONTEXT_CORES = 6)
      expect_equal(spim_control_cores(), 6)
    }))
})


test_that("Sensible parallel control", {
  skip_if_not_installed("mockery")
  mock_cores <- mockery::mock(32, cycle = TRUE)
  mockery::stub(spim_control_parallel, "spim_control_cores", mock_cores)
  expect_equal(
    spim_control_parallel(4, FALSE, NULL, FALSE, FALSE, FALSE),
    list(n_threads_total = 32, n_workers = 1))
  mockery::expect_called(mock_cores, 1)
  expect_equal(mockery::mock_args(mock_cores)[[1]], list())

  expect_equal(
    spim_control_parallel(4, FALSE, 32, FALSE, FALSE, FALSE),
    list(n_threads_total = 32, n_workers = 1))
  expect_equal(
    spim_control_parallel(4, TRUE, 32, FALSE, FALSE, FALSE),
    list(n_threads_total = 32, n_workers = 4))
  expect_equal(
    spim_control_parallel(8, TRUE, 32, FALSE, FALSE, FALSE),
    list(n_threads_total = 32, n_workers = 4))
  expect_equal(
    spim_control_parallel(3, TRUE, 32, FALSE, FALSE, FALSE),
    list(n_threads_total = 32, n_workers = 2))
  expect_equal(
    spim_control_parallel(1, TRUE, 32, FALSE, FALSE, FALSE),
    list(n_threads_total = 32, n_workers = 1))

  expect_message(
    spim_control_parallel(8, TRUE, 32, FALSE, FALSE, TRUE),
    "Running on 4 workers with 32 threads")
})


test_that("Overall spim control", {
  expect_message(
    ctl <- spim_control(TRUE, 4, n_threads = 16),
    "Running on 4 workers with 16 threads")
  expect_setequal(names(ctl), c("pmcmc", "particle_filter"))
  expect_s3_class(ctl$pmcmc, "pmcmc_control")
})


test_that("spim control can contain restart dates", {
  ctl <- suppressMessages(spim_control(TRUE, 4, n_threads = 16))
  expect_null(ctl$pmcmc$save_restart)

  dates <- c("2021-01-01", "2021-02-01")
  ctl <- suppressMessages(
    spim_control(TRUE, 4, n_threads = 16, date_restart = dates))
  expect_equal(ctl$pmcmc$save_restart, sircovid::sircovid_date(dates))
})


test_that("spim control short run is shorter", {
  ctl_long <- suppressMessages(spim_control(FALSE, 4, n_threads = 16))
  ctl_short <- suppressMessages(spim_control(TRUE, 4, n_threads = 16))

  expect_lt(ctl_short$pmcmc$n_steps, ctl_long$pmcmc$n_steps)

  expect_equal(ctl_short$pmcmc$n_burnin, 7)
  expect_equal(ctl_long$pmcmc$n_burnin, 503)

  expect_equal(ctl_short$pmcmc$n_steps_retain, 3)
  expect_equal(ctl_long$pmcmc$n_steps_retain, 250)

  expect_equal(ctl_short$particle_filter$n_particles, 10)
  expect_equal(ctl_long$particle_filter$n_particles, 192)
})


test_that("Allow disabling workers for deterministic fit", {
  suppressMessages(
    withr::with_envvar(c(MC_CORES = 2), {
      control1 <- spimalot::spim_control(
        TRUE, 2, TRUE, n_mcmc = 100,
        burnin = 5, workers = TRUE)
    }))
  expect_equal(control1$pmcmc$n_workers, 2)
  suppressMessages(
    withr::with_envvar(c(MC_CORES = 2), {
    control2 <- spimalot::spim_control(
      TRUE, 4, TRUE, n_mcmc = 100,
      burnin = 5, workers = FALSE)
    }))
  expect_equal(control2$pmcmc$n_workers, 1)
})


test_that("Can change number of samples", {
  suppressMessages(
    ctrl <- spim_control(FALSE, 8, n_sample = 200, n_mcmc = 1000, burnin = 500))
  expect_equal(ctrl$pmcmc$n_steps_retain, 25) # i.e., 25 * 8 == 200
})


test_that("parallel control", {
  expect_equal(
    spim_control_parallel(8, TRUE, 16, FALSE, FALSE, FALSE),
    list(n_threads_total = 16, n_workers = 4))
  expect_equal(
    spim_control_parallel(8, TRUE, 16, TRUE, FALSE, FALSE),
    list(n_threads_total = 8, n_workers = 8))
  expect_equal(
    spim_control_parallel(8, TRUE, 16, TRUE, TRUE, FALSE),
    list(n_threads_total = 16, n_workers = 4))
  expect_equal(
    spim_control_parallel(8, TRUE, 10, TRUE, TRUE, FALSE),
    list(n_threads_total = 12, n_workers = 4))

  expect_equal(
    spim_control_parallel(8, FALSE, 16, FALSE, FALSE, FALSE),
    list(n_threads_total = 16, n_workers = 1))
  expect_equal(
    spim_control_parallel(8, FALSE, 16, TRUE, FALSE, FALSE),
    list(n_threads_total = 16, n_workers = 1))
  expect_equal(
    spim_control_parallel(8, FALSE, 16, TRUE, TRUE, FALSE),
    list(n_threads_total = 16, n_workers = 1))
  expect_equal(
    spim_control_parallel(8, FALSE, 10, TRUE, TRUE, FALSE),
    list(n_threads_total = 10, n_workers = 1))
})


test_that("save path into control", {
  suppressMessages(ctl <- spim_control(TRUE, 4, mcmc_path = "mcmc"))
  expect_equal(ctl$pmcmc$path, "mcmc")
})


test_that("Don't rerun deterministic models", {
  expect_equal(
    suppressMessages(
      spim_control(FALSE, 4, deterministic = FALSE)$pmcmc$rerun_every),
    100)
  expect_equal(
    suppressMessages(
      spim_control(FALSE, 4, deterministic = TRUE)$pmcmc$rerun_every),
    Inf)
})


test_that("Control compiled compare", {
  ctl <- spim_control(TRUE, 4, verbose = FALSE)
  expect_false(ctl$particle_filter$compiled_compare)

  ctl <- spim_control(TRUE, 4, compiled_compare = TRUE, verbose = FALSE)
  expect_true(ctl$particle_filter$compiled_compare)
})
