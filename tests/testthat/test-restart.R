test_that("Can create restart function", {
  m <- matrix(1:50, 5)
  p <- list(steps_per_day = 4)
  f <- spim_restart_initial(m, 100, FALSE)
  expect_true(is.function(f))
  res <- f(NULL, 3, p)

  expect_equal(dim(res$state), c(5, 3))
  expect_equal(res$step, 400)
  expect_true(all(res$state[1, ] %in% m[1, ]))
  expect_true(all(diff(res$state) == 1))
})
