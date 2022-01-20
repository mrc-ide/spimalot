fit_process_restart <- function(samples, parameters) {
  if (is.null(samples$restart)) {
    return(NULL)
  }

  ## TODO: this is now done here and also later, but it's fairly fast
  pars <- spim_fit_parameters(samples, parameters)
  ## TODO: could pass in samples$info$pars, but somehow feels more
  ## awkward to do so
  pars$prior <- fit_process_restart_priors(samples$pars, pars)
  pars$sample <- samples$pars
  class(pars) <- "spim_pars_pmcmc"

  pars$base <- parameters$base

  list(state = samples$restart,
       info = samples$info,
       pars = pars)
}


fit_process_restart_priors <- function(values, parameters, nms) {
  nms <- colnames(values)
  stopifnot(all(nms %in% parameters$info$name))

  info <- parameters$info
  prior <- parameters$prior

  if (length(dim(values)) == 3) {
    wrapper_multiregion <- function(nm, region) {
      if (is.na(region)) {
        x <- values[, nm, 1]
        prior <- prior[prior$name == nm & is.na(prior$region), ]
        info <- info[info$name == nm & is.na(info$region), ]
      } else {
        x <- values[, nm, region]
        prior <- prior[prior$name == nm & prior$region == region, ]
        info <- info[info$name == nm & info$region == region, ]
      }
      if (length(unique(x)) > 1) {
        prior <- fit_prior(x, info, prior)
      }
      prior
    }

    nms_fixed <- info$name[is.na(info$region)]
    nms_varied <- setdiff(info$name, nms_fixed)
    region <- last(dimnames(values))
    prior_fixed <- lapply(nms_fixed, wrapper_multiregion, NA)
    prior_varied <- unlist(lapply(region, function(r)
      lapply(nms_varied, wrapper_multiregion, r)),
      FALSE, FALSE)
    res <- c(prior_fixed, prior_varied)
  } else {
    wrapper_single <- function(nm) {
      x <- values[, nm]
      info <- info[info$name == nm, ]
      prior <- prior[prior$name == nm, ]
      if (length(unique(x)) > 1) {
        prior <- fit_prior(x, info, prior)
      }
      prior
    }
    res <- lapply(nms, wrapper_single)
  }

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
