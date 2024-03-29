##' Load the pmcmc parameters from disk. We expect three files; one
##' for the overall parameters (`info`), one with details of the priors
##' (`prior`) and one describing the proposal kernel (`proposal`).
##'
##' @title Load a set of parameters for the pmcmc
##'
##' @param path Directory where the csv files are found.
##'
##' @param info Filename for the parameter info, relative to `path`
##'
##' @param prior Filename for the parameter priors, relative to `path`
##'
##' @param proposal Filename for the parameter proposal, relative to `path`
##'
##' @return A spim_pars_pmcmc object
##' @export
spim_pars_pmcmc_load <- function(path, info = "info.csv", prior = "prior.csv",
                                 proposal = "proposal.csv") {
  assert_file_exists(path)
  assert_file_exists(info, path)
  assert_file_exists(prior, path)
  assert_file_exists(proposal, path)
  ret <- list(info = read_csv(file.path(path, info)),
              prior = read_csv(file.path(path, prior)),
              proposal = read_csv(file.path(path, proposal)))
  ## TODO: at this point we might split into regions and verify
  ## everything, rather than the fiddle done in spim_pars_single. Easy
  ## to tweak later though.
  class(ret) <- "spim_pars_pmcmc"
  ret
}


##' Write pmcmc parameters out after fitting and aggregation (inverse
##' of [spimalot::spim_pars_pmcmc_load])
##'
##' @title Write out parameters
##'
##' @param p A list with elements `info`, `prior` and `proposal`
##'
##' @param path Directory to write to
##'
##' @return Nothing, called for side effects
##' @export
spim_pars_pmcmc_save <- function(p, path) {
  stopifnot(all(c("info", "prior", "proposal") %in% names(p)))
  dir.create(path, FALSE, TRUE)
  write_csv(p$info, file.path(path, "info.csv"))
  write_csv(p$prior, file.path(path, "prior.csv"))
  write_csv(p$proposal, file.path(path, "proposal.csv"))
}


##' Validate vector of beta dates
##'
##' @title Validate vector of beta dates
##'
##' @param beta_date Vector of beta dates
##'
##' @export
spim_pars_check_beta_date <- function(beta_date) {
  assert_date_string(beta_date)
  assert_increasing(as.Date(beta_date), name = "beta_date")
  beta_date
}


spim_pars_info_single <- function(region, info) {
  assert_has_names(info, c("region", "include",
                           "name", "initial", "max", "integer"))
  if (!(region %in% info$region)) {
    stop(sprintf("Did not find region '%s' in parameter info", region))
  }
  info <- info[info$region == region & info$include, ]
  info <- info[setdiff(names(info), c("region", "include"))]
  rownames(info) <- NULL
  assert_unique(info$names)
  info
}


spim_pars_info_nested <- function(region, info) {
  assert_has_names(info, c("region", "include", "name", "initial",
                           "max", "integer"))
  msg <- setdiff(region, info$region)
  if (length(msg) > 0) {
    stop(sprintf("Did not find region %s in parameter info",
                 paste(squote(msg), collapse = ", ")))
  }

  nms_varied <- unique(info$name[!is.na(info$region)])
  nms_fixed <- info$name[is.na(info$region)]

  if (length(nms_fixed) == 0) {
    stop("Did not find any fixed parameters, seems unlikely")
  }

  info_fixed <- info[is.na(info$region) & info$include, ]
  info_fixed$region <- NULL
  info_fixed$include <- NULL

  ## This won't include great error messages if things go wrong.
  f <- function(x) {
    assert_setequal(x$region, region)
    x <- x[match(region, x$region), setdiff(names(x), "include")]
    ret <- as.list(x)
    ret$name <- ret$name[[1]]
    ret
  }
  info_varied <- info[info$region %in% region & info$include, ]
  info_varied <- lapply(split(info_varied, info_varied$name), f)

  list(fixed = info_fixed, varied = info_varied)
}


spim_pars_prior_single <- function(region, info, prior) {
  prior_cols <- c("region", "type", "name", "gamma_scale", "gamma_shape",
                  "beta_shape1", "beta_shape2")
  assert_has_names(prior, prior_cols)
  prior <- prior[prior$region == region & prior$name %in% info$name, ]
  prior <- prior[setdiff(names(prior), "region")]
  rownames(prior) <- NULL
  assert_setequal(prior$name, info$name)
  prior
}


spim_pars_prior_nested <- function(region, info, prior) {
  prior_cols <- c("region", "type", "name", "gamma_scale",
                  "gamma_shape", "beta_shape1", "beta_shape2")
  assert_has_names(prior, prior_cols)

  prior_fixed <- prior[is.na(prior$region), ]
  assert_setequal(prior_fixed$name, info$fixed$name)
  prior_fixed$region <- NULL

  f <- function(x) {
    assert_setequal(x$region, region)
    lapply(split(x, x$region), as.list)
  }
  prior_varied <- prior[prior$region %in% region, ]
  prior_varied <- lapply(split(prior_varied, prior_varied$name), f)

  list(fixed = prior_fixed, varied = prior_varied)
}


spim_pars_proposal_nested <- function(region, info, proposal, kernel_scaling) {
  proposal_fixed <- proposal[is.na(proposal$region), ]
  assert_setequal(proposal_fixed$name, info$fixed$name)
  i <- match(info$fixed$name, proposal_fixed$name)
  proposal_fixed <- as.matrix(proposal_fixed[i, info$fixed$name]) *
    kernel_scaling
  rownames(proposal_fixed) <- info$fixed$name

  f <- function(x) {
    assert_setequal(x$name, names(info$varied))
    i <- match(names(info$varied), x$name)
    m <- as.matrix(x[i, names(info$varied)]) * kernel_scaling
    rownames(m) <- names(info$varied)
    m
  }
  proposal_varied <- proposal[proposal$region %in% region, ]
  proposal_varied <- lapply(split(proposal_varied, proposal_varied$region), f)
  nms <- c(dimnames(proposal_varied[[1]]), list(region))
  proposal_varied <-
    array(unlist(proposal_varied, FALSE, FALSE), lengths(nms), nms)

  list(fixed = proposal_fixed, varied = proposal_varied)
}


spim_pars_proposal_single <- function(region, info, proposal, kernel_scaling) {
  proposal <- proposal[proposal$region == region, ]
  assert_setequal(proposal$name, info$name)
  assert_setequal(setdiff(names(proposal), c("name", "region")), info$name)
  proposal <- as.matrix(proposal[match(info$name, proposal$name), info$name])
  rownames(proposal) <- info$name
  proposal * kernel_scaling
}


make_prior <- function(d) {
  if (d$type == "gamma") {
    ## TODO: as_duration was droppd from here as never used, but if it
    ## is, then we'd transform p to 1/p
    shape <- d$gamma_shape
    scale <- d$gamma_scale
    function(p) {
      dgamma(p, shape = shape, scale = scale, log = TRUE)
    }
  } else if (d$type == "beta") {
    shape1 <- d$beta_shape1
    shape2 <- d$beta_shape2
    function(p) {
      dbeta(p, shape1 = shape1, shape2 = shape2, log = TRUE)
    }
  } else if (d$type == "null") {
    NULL
  } else {
    stop("Unknown prior type")
  }
}


##' Add a parameter to a `spim_pars_mcmc` object
##'
##' @title Add parameter
##'
##' @param pars Original `spim_pars_pmcmc` object
##'
##' @param name Name of the new parameter
##'
##' @param initial Initial value for the new parameter
##'
##' @param min Minimum value
##'
##' @param max Maximum value
##'
##' @param proposal_variance Variance for the proposal VCV
##'
##' @param prior List of values for the prior (must include "type",
##'   then other values for the prior data.frame)
##'
##' @param integer Logical, indicates if integer
##'
##' @param include Logical, indicates if fitted
##'
##' @export
spim_add_par <- function(pars, name, initial, min, max, proposal_variance,
                         prior, integer = FALSE, include = TRUE) {
  assert_is(pars, "spim_pars_pmcmc")

  new_par <- data_frame(
    region = unique(pars$info$region),
    name = name,
    initial = initial,
    min = min,
    max = max,
    integer = integer,
    include = include)
  pars$info <- rbind(pars$info, new_par)
  pars$info <- pars$info[order(pars$info$region, pars$info$name), ]

  proposal <- pars$proposal
  proposal[[name]] <- 0
  ## next line expands regions
  new_prop <- proposal[proposal$name == proposal$name[1], ]
  new_prop$name <- name
  new_prop[, -c(1, 2)] <- 0
  new_prop[[name]] <- proposal_variance

  proposal <- rbind(proposal, new_prop)
  proposal <- proposal[order(proposal$region, proposal$name), ]

  nms <- c(names(proposal)[1:2], sort(names(proposal)[-c(1, 2)]))
  pars$proposal <- proposal[nms]

  new_prior <- pars$prior[pars$prior$name == pars$prior$name[1], ]
  new_prior[setdiff(names(new_prior), c("region", "name"))] <- NA
  stopifnot(all(names(prior) %in% names(new_prior)))
  new_prior$name <- name
  for (i in names(prior)) {
    new_prior[[i]] <- prior[[i]]
  }

  pars$prior <- rbind(pars$prior, new_prior)
  pars$prior <- pars$prior[order(pars$prior$region, pars$prior$name), ]

  class(pars) <- "spim_pars_pmcmc"

  pars
}


##' @export
##' @rdname spim_add_par
spim_add_par_beta <- function(pars) {
  re <- "^beta([0-9]+)$"
  n <- max(as.integer(
    sub(re, "\\1", grep(re, pars$info$name, value = TRUE))))
  last <- paste0("beta", n)
  name <- paste0("beta", n + 1)
  info <- pars$info[match(last, pars$info$name), ]
  proposal_variance <- pars$proposal[[last]][match(last, pars$proposal$name)]
  prior <- as.list(pars$prior[match(last, pars$prior$name), ])
  prior <- prior[setdiff(names(prior), c("region", "name"))]
  spim_add_par(pars, name, info$initial, info$min, info$max,
               proposal_variance, prior)
}


spim_pars_mcmc_single <- function(info, prior, proposal, transform) {
  pars_mcmc <- Map(
    mcstate::pmcmc_parameter,
    name = info$name,
    initial = info$initial,
    min = info$min,
    max = info$max,
    integer = info$integer,
    prior = lapply(split(prior, prior$name), make_prior))

  ret <- mcstate::pmcmc_parameters$new(pars_mcmc, proposal, transform)

  ## Try and transform a single case and see if it works:
  ret$model(ret$initial())

  ret
}


spim_pars_mcmc_nested <- function(info, prior, proposal, transform) {
  prior_fixed <- lapply(split(prior$fixed, prior$fixed$name),
                        make_prior)
  pars_fixed <- Map(
    mcstate::pmcmc_parameter,
    name = info$fixed$name,
    initial = info$fixed$initial,
    min = info$fixed$min,
    max = info$fixed$max,
    integer = info$fixed$integer,
    prior = prior_fixed)

  ## These could be done per region, or could be done separately.  I
  ## don't really care.  The current approach is clearly pretty messed
  ## here.
  prior_varied <- lapply(prior$varied, function(x)
    lapply(x, make_prior))
  regions <- info$varied[[1]]$region
  pars_varied <- lapply(names(info$varied), function(i)
    mcstate::pmcmc_varied_parameter(
      name = i,
      populations = regions,
      initial = info$varied[[i]]$initial,
      min = info$varied[[i]]$min,
      max = info$varied[[i]]$max,
      integer = info$varied[[i]]$integer[[1]],
      prior = prior_varied[[i]]))
  names(pars_varied) <- names(info$varied)

  ret <- mcstate::pmcmc_parameters_nested$new(
    parameters = c(pars_varied, pars_fixed),
    proposal_varied = proposal$varied,
    proposal_fixed = proposal$fixed,
    transform = transform)

  ## Try and transform a single case and see if it works:
  ret$model(ret$initial())

  ret
}
