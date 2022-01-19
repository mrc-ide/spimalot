test_that("check_model_type identifies correct models", {
  expect_silent(spim_check_model_type("BB"))
  expect_silent(spim_check_model_type("NB"))
  expect_error(spim_check_model_type("BN"),
               "Expected 'BN' to be one of 'BB' or 'NB'")
})


test_that("check_region identifies correct regions", {
  expect_silent(check_region("london"))
  expect_error(check_region("islington"),
               "region must be one of")
})


test_that("spim_check_region works", {
  expect_equal(spim_check_region("london", FALSE), "london")
  expect_error(spim_check_region("islington", FALSE),
               "region must be one of")
  expect_error(spim_check_region("all", FALSE),
               "region must be one of")

  expect_equal(spim_check_region("all", TRUE),
               sircovid::regions("all"))
  expect_equal(spim_check_region("england", TRUE),
               sircovid::regions("england"))
  expect_equal(spim_check_region("test", TRUE),
               c("london", "south_west"))
  expect_error(spim_check_region("london", TRUE),
               "Invalid value 'london' for 'region' when multiregion = TRUE")
})
