##' Convert regions to different types of names
##'
##' @title Convert region names
##'
##' @param region Vector of standard region names (london, scotland,
##'   england, uk)
##'
##' @param type Convertion type. Current can be one of "name", "code" or
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
##' @return Returns an rrq controller object if found, otherwise errors
##' @export
spim_rrq_controller <- function(root = here::here()) {
  queue_id <- Sys.getenv("CONTEXT_ID", "")
  if (queue_id == "") {
    stop("No rrq controller found")
  } else {
    message(sprintf("Found rrq controller for queue '%s'", queue_id))
    message(sprintf("Using root directory '%s'", root))
    withr::with_dir(root, rrq::rrq_controller(queue_id))
  }
}
