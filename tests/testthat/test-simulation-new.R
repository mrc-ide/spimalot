test_that("can construct a control output object", {
  control <- spim_simulate_control_output(c("deaths", "admissions"))
  expect_s3_class(control, "spim_simulate_control_output")

  expect_equal(control$keep, c("deaths", "admissions"))
  expect_true(control$time_series)
  expect_true(control$rt)
  expect_true(control$rt_weighted)
  expect_equal(control$rt_type, c("eff_Rt_general", "Rt_general"))
  expect_true(control$state_by_age)
  expect_true(control$vaccination)
})


test_that("Can create minimal control object", {
  control <- spim_simulate_control_output(character(0), time_series = FALSE,
                                          rt = FALSE, rt_weighted = FALSE,
                                          rt_type = "Rt_general",
                                          state_by_age = FALSE,
                                          vaccination = FALSE)
  expect_s3_class(control, "spim_simulate_control_output")

  expect_equal(control$keep, character(0))
  expect_false(control$time_series)
  expect_false(control$rt)
  expect_false(control$rt_weighted)
  expect_equal(control$rt_type, "Rt_general")
  expect_false(control$state_by_age)
  expect_false(control$vaccination)
})


test_that("Disallow invalid values", {
  expect_error(
    spim_simulate_control_output(TRUE),
    "'keep' must be character")
  expect_error(
    spim_simulate_control_output(character(0), rt_type = "yes, please"),
    "Invalid value for 'rt_type': 'yes, please'")
  expect_error(
    spim_simulate_control_output(character(0), rt = FALSE, rt_weighted = TRUE),
    "If rt_weighted is TRUE, then rt must be TRUE")
})
