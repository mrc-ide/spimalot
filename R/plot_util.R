add_alpha <- function(col, alpha = 1) {
  apply(sapply(col, grDevices::col2rgb) / 255, 2,
        function(x)
          grDevices::rgb(x[1], x[2], x[3], alpha = alpha))
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
