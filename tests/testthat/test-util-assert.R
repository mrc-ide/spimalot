test_that("assert_has_names", {
  object <- list(a = 1, b = 2, c = 3)
  expect_silent(assert_has_names(object, "a"))
  expect_silent(assert_has_names(object, c("a", "b")))
  expect_error(assert_has_names(object, c("a", "x", "y")),
               "Required names missing from 'object': 'x', 'y'")
})

test_that("assert_scalar", {
  expect_error(assert_scalar(NULL), "must be a scalar")
  expect_error(assert_scalar(numeric(0)), "must be a scalar")
  expect_error(assert_scalar(1:2), "must be a scalar")
})


test_that("assert_character", {
  expect_error(assert_character(1), "must be character")
  expect_error(assert_character(TRUE), "must be character")
})


test_that("assert_numeric", {
  expect_error(assert_numeric("one"), "must be numeric")
  expect_error(assert_numeric(TRUE), "must be numeric")
})


test_that("assert_logical", {
  expect_silent(assert_logical(TRUE))
  expect_silent(assert_logical(c(TRUE, FALSE)))
  expect_error(assert_logical("one"), "must be a logical")
  expect_error(assert_logical(1), "must be a logical")
})


test_that("assert_scalar_logical", {
  expect_silent(assert_scalar_logical(TRUE))
  expect_error(assert_scalar_logical(c(TRUE, FALSE)))
})


test_that("assert_scalar_numeric", {
  expect_silent(assert_scalar_numeric(1))
  expect_error(assert_scalar_numeric(c(1, 2)))
})


test_that("assert_is", {
  expect_error(assert_is("x", "foo"), "must be a foo")
  expect_silent(assert_is(structure("x", class = "foo"), "foo"))
})


test_that("match_value", {
  expect_error(match_value("foo", letters), "must be one of")
  expect_silent(match_value("a", letters))
})


test_that("assert_file_exists", {
  path <- tempfile(tmpdir = normalizePath(tempdir(), mustWork = TRUE))
  expect_error(assert_file_exists(path), "File '.+' does not exist$")
  expect_error(assert_file_exists(basename(path), dirname(path)),
               "File '.+' does not exist \\(relative to '.+'\\)")

  writeLines(character(0), path)
  expect_silent(assert_file_exists(path))
  expect_silent(assert_file_exists(basename(path), dirname(path)))
})


test_that("assert_length", {
  object <- 1:5
  expect_silent(assert_length(object, 5))
  expect_error(assert_length(object, 3),
               "'object' must have length 3")
})


test_that("assert_scalar_positive_integer", {
  expect_silent(assert_scalar_positive_integer(10))
  expect_error(assert_scalar_positive_integer(-10),
               "'-10' must be at least 1")
  expect_error(assert_scalar_positive_integer(0),
               "'0' must be at least 1")
  expect_error(assert_scalar_positive_integer(NA_integer_),
               "must be non-NA")
  expect_error(assert_scalar_positive_integer(1:2),
               "must be a scalar")
  expect_error(assert_scalar_positive_integer(1.5),
               "must be an integer")
})


test_that("assert_unique", {
  expect_silent(assert_unique(1:10))
  object <- rep(1:10, 3)
  expect_error(assert_unique(object),
               "'object' contains duplicate elements")
})


test_that("assert setequal", {
  x <- 1:3
  y <- 2:4
  expect_silent(assert_setequal(x, rev(x)))
  expect_error(assert_setequal(x, y),
               "'x' and 'y' are not setequal")
})
