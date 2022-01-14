##' Check that a model type is valid (i.e., one of `NB` or `BB`)
##'
##' @title Check a model type
##' @param model_type String
##' @export
spim_check_model_type <- function(model_type) {
  if (!(model_type %in% c("BB", "NB"))) {
    stop(sprintf("Expected '%s' to be one of 'BB' or 'NB'", model_type),
         call. = FALSE)
  }
}


##' Check that a model type is valid (i.e., one of `NB` or `BB`)
##'
##' @title Check a model type
##' @param sircovid_model String
##' @export
spim_check_sircovid_model <- function(sircovid_model) {
  if (!(sircovid_model %in% c("carehomes", "lancelot"))) {
    stop(sprintf("Expected '%s' to be one of 'carehomes' or 'lancelot'",
                 sircovid_model),
         call. = FALSE)
  }
}


check_region <- function(region) {
  match_value(region, sircovid::regions("all"))
}

##' Validate (and possibly expand) a region for a fit.
##'
##' @title Validate and expand regions
##'
##' @param region A single string indicating a region or group of
##'   regions.  For `multiregion = FALSE` this must be a suitable NHS
##'   region.  When `multiregion = TRUE` this must be one of `all`,
##'   `england` or `test`
##'
##' @param multiregion Logical, indicating if this is a multiregion fit
##'
##' @return A character vector of regions
##' @export
spim_check_region <- function(region, multiregion) {
  assert_scalar_character(region)
  if (multiregion) {
    if (region %in% c("all", "england")) {
      region <- sircovid::regions(region)
    } else if (region == "test") {
      ## We could use a different set here really, or something more
      ## clever
      region <- c("london", "south_west")
    } else {
      stop(sprintf("Invalid value '%s' for 'region' when multiregion = TRUE",
                   region))
    }
  } else {
    match_value(region, sircovid::regions("all"))
  }
  region
}
