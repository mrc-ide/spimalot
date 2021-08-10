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
