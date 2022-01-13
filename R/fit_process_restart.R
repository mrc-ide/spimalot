fit_process_restart <- function(samples, parameters) {
  if (is.null(samples$restart)) {
    return(NULL)
  }

  pars <- spim_fit_parameters(samples, parameters)
  pars$prior <- fit_process_restart_priors(samples$pars, pars)
  pars$sample <- samples$pars
  class(pars) <- "spim_pars_pmcmc"

  pars$base <- parameters$base

  list(state = samples$restart,
       info = samples$info,
       pars = pars)
}


fit_process_restart_priors <- function(values, parameters) {
  nms <- colnames(values)
  stopifnot(all(nms %in% parameters$info$name))

  wrapper <- function(nm) {
    x <- values[, nm]
    info <- parameters$info[parameters$info$name == nm, ]
    prior <- parameters$prior[parameters$prior$name == nm, ]
    if (length(unique(x)) > 1) {
      new_prior <- fit_prior(x, info, prior)
    } else {
      ## if all values of the parameter are the same, keep the original prior
      new_prior <- prior
    }
    new_prior
  }

  res <- lapply(nms, wrapper)
  do.call("rbind", res)
}


fit_prior <- function(x, info, prior) {
  ## Can't update null priors; these will eventually be replaced in
  ## the underlying priors task, but this requires work from Marc.
  if (prior$type == "null") {
    return(prior)
  }

  if (prior$type == "gamma") {
    ll <- function(theta) {
      -suppressWarnings(
        sum(dgamma(x, scale = theta[1], shape = theta[2], log = TRUE)))
    }
    cols <- c("gamma_scale", "gamma_shape")
  } else if (prior$type == "beta") {
    ll <- function(theta) {
      suppressWarnings(
        -sum(dbeta(x, shape1 = theta[1], shape2 = theta[2], log = TRUE)))
    }
    cols <- c("beta_shape1", "beta_shape2")
  } else {
    stop(sprintf("Unsupported prior type '%s'", prior$type))
  }

  fit <- optim(unlist(prior[cols]), ll)
  prior[cols] <- fit$par
  prior
}
