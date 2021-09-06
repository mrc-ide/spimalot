## FIXME - DOCUMENT
set_nested_pars <- function(pars, fixed, varied) {
  if (is.null(fixed) && is.null(varied)) {
    stop("At least one of 'fixed' and 'varied' should be non-NULL")
  }

  parameters <- unique(pars$info$name)

  if (is.null(fixed)) {
    fixed <- setdiff(parameters, varied)
  } else if (is.null(varied)) {
    varied <- setdiff(parameters, fixed)
  }

  stdif <- setdiff(varied, parameters)

  if (length(stdif) > 0) {
    stop(sprintf("Parameters set as 'varied' don't exist in pars$info: %s",
                 str_collapse(stdif)))
  }

  stdif <- setdiff(fixed, parameters)

  if (length(stdif) > 0) {
    stop(sprintf("Parameters set as 'fixed' don't exist in pars$info: %s",
                 str_collapse(stdif)))
  }

  stdif <- setdiff(c(varied, fixed), parameters)

  if (length(stdif) > 0) {
    stop(sprintf("Parameters in pars$info not set as fixed or varied: %s",
                 str_collapse(stdif)))
  }


  ret <- average_nested_parameters(pars, fixed, varied)
  class(ret) <- "spim_pars_pmcmc"
  list(parameters = ret, fixed = fixed, varied = varied)
}


average_nested_parameters <- function(pars, fixed, varied) {

  proposal <- pars$proposal
  info <- pars$info
  prior <- pars$prior

  ## zero out relationship between fixed/varied parameters in proposal
  proposal[proposal$name %in% fixed, varied] <- 0
  proposal[proposal$name %in% varied, fixed] <- 0

  for (f in fixed) {
    ## take the mean of the fixed with itself in proposal
    proposal[proposal$name == f, f] <- mean(proposal[proposal$name == f, f])
    ## and zero out the others in proposal
    proposal[proposal$name == f, setdiff(fixed, f)] <- 0
    proposal[proposal$name %in% setdiff(fixed, f), f] <- 0

    ## take mean of fixed over regions in info
    info[info$name == f, "initial"] <- mean(info[info$name == f, "initial"])

    cols <- c("gamma_scale", "gamma_shape", "beta_shape1", "beta_shape2")
    n_regions <- length(unique(prior$region))
    prior[prior$name == f, cols] <-
      matrix(colMeans(prior[prior$name == f, cols], na.rm = TRUE),
             n_regions, 4, TRUE)
  }

  list(info = info, prior = prior, proposal = proposal)
}