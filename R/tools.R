##' Convert regions to different types of names
##'
##' @title Convert region names
##'
##' @param region Vector of standard region names (london, scotland,
##'   england, uk)
##'
##' @param type Conversion type. Current can be one of "name", "code" or
##'   "upper"
##'
##' @return A vector of new names
##' @export
spim_region_name <- function(region, type = "name") {
  if (type %in% c("name", "upper")) {
    map <- c(
      london = "London",
      east_of_england = "East of England",
      midlands = "Midlands",
      north_east_and_yorkshire = "North East and Yorkshire",
      north_west = "North West",
      south_east = "South East",
      south_west = "South West",
      scotland = "Scotland",
      wales = "Wales",
      northern_ireland = "Northern Ireland",
      england = "England",
      uk = "United Kingdom")

    if (type == "upper") {
      map <- toupper(map)
    }
  } else if (type == "code") {
    map <- c(
      london = "LON",
      east_of_england = "EE",
      midlands = "MID",
      north_east_and_yorkshire = "NE",
      north_west = "NW",
      south_east = "SE",
      south_west = "SW",
      scotland = "SCO",
      wales = "WAL",
      northern_ireland = "NI",
      england = "ENG",
      uk = "UK")
  } else {
    stop(sprintf("Unknown region name type '%s'", type))
  }

  res <- map[region]
  if (any(is.na(res))) {
    stop("Invalid region")
  }
  unname(res)
}


##' Finds an rrq controller if available
##'
##' @title Find rrq controller
##'
##' @param root Root of the orderly project (used to anchor the rrq
##'   file store).
##'
##' @return Returns an rrq controller object if found, otherwise errors
##' @export
spim_rrq_controller <- function(root = here::here()) {
  queue_id <- Sys.getenv("CONTEXT_ID", "")
  if (queue_id == "") {
    stop("No rrq controller found")
  } else {
    message(sprintf("Found rrq controller for queue '%s'", queue_id))
    message(sprintf("Using root directory '%s'", root))
    if (packageVersion("rrq") < "0.5.0") {
      withr::with_dir(root, rrq::rrq_controller(queue_id))
    } else {
      withr::with_dir(root, rrq::rrq_controller$new(queue_id))
    }
  }
}


##' Round to nearest given number
##' @title Round to given number
##' @param x Number to round
##' @param y Number to round to
##' @export
##' @examples
##'  # Round to nearest 100
##'  round(146, 100)
##'  # Round to nearest 1000
##'  round(7849, 1000)
round_to <- function(x, y) {
  stopifnot(y >= 1)
  round(x / y) * y
}