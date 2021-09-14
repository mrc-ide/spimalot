fit_process_restart <- function(samples, parameters, data, control) {
  if (is.null(samples$restart)) {
    return(NULL)
  }

  ## Add 1 to burnin to account for removal of initial parameters
  samples <- mcstate::pmcmc_thin(samples, control$burnin + 1L, control$thin)

  pars <- spim_fit_parameters(samples, parameters)
  pars$prior <- fit_process_restart_priors(samples$pars, pars)

  pars$sample <- drop(samples$pars)
  class(pars) <- "spim_pars_pmcmc"

  samples$restart$state <- drop_dimension(samples$restart$state, 3)

  list(state = samples$restart,
       info = samples$info,
       pars = pars,
       data = data)
}


fit_process_restart_priors <- function(values, parameters) {
  nms <- unique(rownames(values))
  stopifnot(all(nms %in% parameters$info$name))

  wrapper <- function(nm, reg) {
    x <- values[nm, reg, ]
    info <- parameters$info[parameters$info$region == reg &
                              parameters$info$name == nm, ]
    prior <- parameters$prior[parameters$prior$region == reg &
                              parameters$prior$name == nm, ]
    fit_prior(x, info, prior)
  }

  res <- lapply(colnames(values), function(x) {
    out <- lapply(nms, wrapper, x)
    do.call(rbind, out)
  })

  do.call(rbind, res)
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
