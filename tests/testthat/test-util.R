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



test_that("add_dimension", {
  dold <- c(2, 5, 8)
  old <- array(runif(80), dold)

  ## check dimensions
  expect_error(add_dimension(old, 0), "d < 1")
  expect_equal(dim(add_dimension(old, 1)), c(1, 2, 5, 8))
  expect_equal(dim(add_dimension(old, 2)), c(2, 1, 5, 8))
  expect_equal(dim(add_dimension(old, 3)), c(2, 5, 1, 8))
  expect_equal(dim(add_dimension(old, 4)), c(2, 5, 8, 1))
  expect_error(add_dimension(old, 5), "d > dim(x) + 1", fixed = TRUE)

  ## check array
  expect_equal(drop(add_dimension(add_dimension(old, 3), 5)),
               old)

  ## can add name if none before
  expect_equal(dimnames(add_dimension(old, 3, "new")),
               list(NULL, NULL, "new", NULL))

  ## can add name if some before
  nms <- list(rep("a", 2), rep("b", 5), rep("c", 8))
  old <- array(runif(80), dold, nms)
  expect_equal(dimnames(add_dimension(old, 3, "new")),
               list(nms[[1]], nms[[2]], "new", nms[[3]]))
})
