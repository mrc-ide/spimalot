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
  mock_cores <- mockery::mock(32, cycle = TRUE)
  mockery::stub(spim_control_parallel, "spim_control_cores", mock_cores)
  expect_equal(
    spim_control_parallel(4, FALSE, verbose = FALSE),
    list(n_threads_total = 32, n_workers = 1))
  mockery::expect_called(mock_cores, 1)
  expect_equal(mockery::mock_args(mock_cores)[[1]], list())

  expect_equal(
    spim_control_parallel(4, FALSE, 32, FALSE),
    list(n_threads_total = 32, n_workers = 1))
  expect_equal(
    spim_control_parallel(4, TRUE, 32, FALSE),
    list(n_threads_total = 32, n_workers = 4))
  expect_equal(
    spim_control_parallel(8, TRUE, 32, FALSE),
    list(n_threads_total = 32, n_workers = 4))
  expect_equal(
    spim_control_parallel(3, TRUE, 32, FALSE),
    list(n_threads_total = 32, n_workers = 2))
  expect_equal(
    spim_control_parallel(1, TRUE, 32, FALSE),
    list(n_threads_total = 32, n_workers = 1))

  expect_message(
    spim_control_parallel(8, TRUE),
    "Running on 4 workers with 32 threads")
})


test_that("Overall spim control", {
  expect_message(
    ctl <- spim_control(TRUE, 4, n_threads = 16),
    "Running on 4 workers with 16 threads")
  expect_setequal(names(ctl), c("pmcmc", "particle_filter", "forecast"))
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

  expect_equal(ctl_short$forecast$burnin, 1)
  expect_equal(ctl_long$forecast$burnin, 500)

  expect_equal(ctl_short$forecast$n_sample, 10)
  expect_equal(ctl_long$forecast$n_sample, 1000)

  expect_equal(ctl_short$particle_filter$n_particles, 10)
  expect_equal(ctl_long$particle_filter$n_particles, 192)
})
