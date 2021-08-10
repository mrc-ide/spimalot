add_alpha <- function(col, alpha = 1) {
  unname(
    apply(sapply(col, grDevices::col2rgb) / 255, 2,
          function(x)
            grDevices::rgb(x[1], x[2], x[3], alpha = alpha)))
}


mix_cols <- function(a, b, weight = 0.5) {
  assert_length(b, length(a))
  assert_scalar_numeric(weight)
  a_rgb <- vapply(a, grDevices::col2rgb, numeric(3))
  b_rgb <- vapply(b, grDevices::col2rgb, numeric(3))
  ret <- a_rgb * (1 - weight) + b_rgb * weight
  grDevices::rgb(ret[1, ], ret[2, ], ret[3, ], maxColorValue = 255)
}


spim_colours <- function() {
  cols <- list(
    green = "#58A449",
    green2 = "#228B22",
    blue = "#278B9A",
    sky_blue = "#C0DDE1",
    sky_blue2 = "#7EBAC2",
    orange = "#E48C2A",
    brown = "#724615",
    puce = "#AF9699",
    purple = "#5C5992",
    cyan = "#699196",
    cyan2 = "#405D61")
  cols$now <- cols$blue
  cols$forecast <- cols$green
  cols
}
