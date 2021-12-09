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


test_that("immutable objects cannot be assigned", {
  x <- list(a = 1, b = 2)
  class(x) <- c("myobj", "immutable")

  expect_error(
    x$a <- 10,
    "Objects of class 'myobj' are immutable and you may not alter them")
  expect_error(
    x[["a"]] <- 10,
    "Objects of class 'myobj' are immutable and you may not alter them")
  expect_error(
    x[[1]] <- 10,
    "Objects of class 'myobj' are immutable and you may not alter them")
  expect_error(
    x[c("b", "d")] <- list(10, 20),
    "Objects of class 'myobj' are immutable and you may not alter them")
  expect_error(
    x[1:2] <- list(10, 20),
    "Objects of class 'myobj' are immutable and you may not alter them")
})
