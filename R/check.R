check_model_type <- function(model_type, name = deparse(substitute(name))) {
  if (!(model_type %in% c("BB", "NB"))) {
    stop(sprintf("Expected '%s' to be one of 'BB' or 'NB'", model_type),
         call. = FALSE)
  }
}


check_region <- function(region) {
  if (region %in% sircovid::regions("all")) {
    return()
  }
}
