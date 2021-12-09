test_that("can construct a control output object", {
  output <- spim_simulate_control_output(c("deaths", "admissions"))
  expect_s3_class(output, "spim_simulate_control_output")

  expect_equal(output$keep, c("deaths", "admissions"))
  expect_true(output$time_series)
  expect_true(output$rt)
  expect_true(output$rt_weighted)
  expect_equal(output$rt_type, c("eff_Rt_general", "Rt_general"))
  expect_true(output$state_by_age)
  expect_true(output$vaccination)
})


test_that("Can create minimal output object", {
  output <- spim_simulate_control_output(character(0), time_series = FALSE,
                                         rt = FALSE, rt_weighted = FALSE,
                                         rt_type = "Rt_general",
                                         state_by_age = FALSE,
                                         vaccination = FALSE)
  expect_s3_class(output, "spim_simulate_control_output")

  expect_equal(output$keep, character(0))
  expect_false(output$time_series)
  expect_false(output$rt)
  expect_false(output$rt_weighted)
  expect_equal(output$rt_type, "Rt_general")
  expect_false(output$state_by_age)
  expect_false(output$vaccination)
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


test_that("Can construct a basic simulation control object", {
  output <- spim_simulate_control_output("deaths")
  grid <- data.frame(vaccine_daily_doses = "mtp")
  parameters <- list(vaccine_daily_doses = 10000)
  control <- spim_simulate_control("flavour",
                                   sircovid::regions("england"),
                                   as.Date("2022-01-01"),
                                   as.Date("2022-03-01"),
                                   parameters,
                                   grid,
                                   output)
  expect_s3_class(control, "spim_simulate_control")
  expect_equal(control$flavour, "flavour")
  expect_equal(control$regions, sircovid::regions("england"))
  expect_equal(control$date_start, as.Date("2022-01-01"))
  expect_equal(control$date_end, as.Date("2022-03-01"))
  expect_equal(control$parameters, parameters)
  expect_equal(control$grid, grid)
  expect_equal(control$output, output)
})


test_that("Disallow invalid values in control", {
  output <- spim_simulate_control_output("deaths")
  grid <- data.frame(vaccine_daily_doses = "mtp")
  parameters <- list(vaccine_daily_doses = 10000)
  start <- as.Date("2022-01-01")
  end <- as.Date("2022-03-01")
  regions <- sircovid::regions("england")
  flavour <- "flavour"

  expect_error(
    spim_simulate_control(TRUE, sircovid::regions("england"), start, end,
                          parameters, grid, output),
    "'flavour' must be character")

  expect_error(
    spim_simulate_control(flavour, "denmark", start, end,
                          parameters, grid, output),
    "Invalid region: 'denmark'")

  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), end, start,
                          parameters, grid, output),
    "'date_end' must be greater than 'date_start'")
  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), TRUE, end,
                          parameters, grid, output),
    "'date_start' must be a Date")
  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), start, TRUE,
                          parameters, grid, output),
    "'date_end' must be a Date")
  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), start, end,
                          NULL, grid, output),
    "'parameters' must be a list")
  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), start, end,
                          parameters, NULL, output),
    "'grid' must be a data.frame")
  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), start, end,
                          parameters, grid[integer(0), , drop = FALSE], output),
    "At least one row required in 'grid'")
  expect_error(
    spim_simulate_control(flavour, sircovid::regions("england"), start, end,
                          parameters, grid, NULL),
    "'output' must be a spim_simulate_control_output")
})
