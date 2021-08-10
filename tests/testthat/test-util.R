test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("data_frame does not mangle names", {
  expect_equal(names(data_frame("name with spaces" = 1:10, b = 1)),
               c("name with spaces", "b"))
})
