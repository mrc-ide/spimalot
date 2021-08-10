test_that("add_alpha", {
  expect_equal(add_alpha("blue", 0.5),
               "#0000FF80")
  expect_equal(add_alpha(c("red", "blue"), 0.2),
               c("#FF000033", "#0000FF33"))
})


test_that("mix", {
  expect_equal(mix_cols("white", "green", 0.5),
               "#7FFF7F")
})


test_that("spim handy colours", {
  cols <- spim_colours()
  expect_type(cols, "list")
  expect_named(cols)
  expect_match(
    unlist(cols, FALSE, FALSE),
    "^#[0-9A-F]{6}$")
})
