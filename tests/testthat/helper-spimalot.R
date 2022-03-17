test_dummy_pars_pmcmc <- function(vars, regions = NULL) {
  regions <- regions %||% sircovid::regions("all")
  info <- lapply(regions, function(r)
    data_frame(region = r, name = vars, initial = runif(length(vars)),
               min = 0, max = 1, integer = FALSE, include = TRUE))
  info <- dplyr::bind_rows(info)

  prior <- lapply(regions, function(r)
    data_frame(region = r, type = "null", name = vars,
               gamma_scale = NA, gamma_shape = NA,
               beta_shape1 = NA, beta_shape2 = NA))
  prior <- dplyr::bind_rows(prior)

  m <- diag(length(vars))
  colnames(m) <- vars
  proposal <- lapply(regions, function(r)
    data_frame(region = r, name = vars, m * runif(length(vars))))
  proposal <- dplyr::bind_rows(proposal)

  ret <- list(info = info, prior = prior, proposal = proposal)
  class(ret) <- "spim_pars_pmcmc"
  ret
}
