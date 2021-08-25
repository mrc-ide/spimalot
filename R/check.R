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


check_region <- function(region) {
  match_value(region, sircovid::regions("all"))
}
