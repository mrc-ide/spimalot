add_alpha <- function(col, alpha = 1) {
  apply(sapply(col, grDevices::col2rgb) / 255, 2,
        function(x)
          grDevices::rgb(x[1], x[2], x[3], alpha = alpha))
}


mix_cols <- function(a, b, weight = 0.5) {
  a_rgb <- vapply(a, grDevices::col2rgb, numeric(3))
  b_rgb <- vapply(b, grDevices::col2rgb, numeric(3))
  ret <- a_rgb * (1 - weight) + b_rgb * weight
  grDevices::rgb(ret[1, ], ret[2, ], ret[3, ], maxColorValue = 255)
}


spim_colours <- function() {
  cols <- list(
    green = "#58A449FF",
    green2 = "#228B22FF",
    blue = "#278B9AFF",
    sky_blue = "#C0DDE1FF",
    sky_blue2 = "#7EBAC2FF",
    orange = "#E48C2AFF",
    brown = "#724615FF",
    puce = "#AF9699FF",
    purple = "#5C5992FF",
    cyan = "#699196FF",
    cyan2 = "#405D61FF")
  cols$now <- cols$blue
  cols$forecast <- cols$green
  cols
}
