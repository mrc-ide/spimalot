test_that("spim_region_name translates regions", {
  expect_equal(spim_region_name(c("london", "east_of_england")),
               c("London", "East of England"))
  expect_equal(spim_region_name(c("london", "east_of_england"), "upper"),
               c("LONDON", "EAST OF ENGLAND"))
  expect_equal(spim_region_name(c("london", "east_of_england"), "code"),
               c("LON", "EE"))
  expect_error(spim_region_name(c("london", "east_of_england"), "other"),
               "Unknown region name type 'other'")
  expect_error(spim_region_name(c("london", "islington")),
               "Invalid region")
})
