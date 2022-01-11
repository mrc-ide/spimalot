##' Load parameters for fitting.
##'
##' There are lots of types of parameters here, depending on the
##' perspective
##'
##' * This function deals with loading in the "logical" parameters;
##'   things like "vaccine efficacy" or "vaccination doses".  These
##'   can be any R object.
##' * Our aim is to create `sircovid` parameters, the result of
##'   `sircovid::lancelot_parameters`.  This will be a list with
##'   structured parmeters that correspond directly to dust/odin
##'   inputs.  Most likely these will be
##'   [mcstate::multistage_parameters] objects that cover multiple
##'   epochs with different model parameters, and these might be
##'   grouped together in a list so that a multiparameter (nested)
##'   filter can be run.
##' * We'll end up working with small(ish) numeric vectors of
##'   parameters from the pmcmc process.  Our transform function will
##'   convert this vector of parameters, using the logical parameters,
##'   into baseline parameters.
##' * The mcmc parameters object that this function will return (as
##'   `$mcmc`) is a [mcstate::pmcmc_parameters] or
##'   [mcstate::pmcmc_parameters_nested] object, which contains all
##'   the information for the above.
##'
##' @title Load parameters
##'
##' @param path The path to load from.  This will be a directory
##'   containing files "info.csv", "prior.csv" and "proposal.csv" (for
##'   [spimalot::spim_pars_pmcmc_load]) but also the baseline
##'   parameters ("base.rds") and a transformation function
##'   ("transform.R")
##'
##' @param region A single region, or vector of regions, to load.
##'
##' @param assumptions The name of any assumptions to apply to filter
##'   the baseline parmeters (e.g., "central")
##'
##' @param kernel_scaling The scaling coefficient for loading the
##'   proposal variance-covariance matrix.
##'
##' @return A list.
##' @export
spim_fit_pars_load <- function(path, region, assumptions, kernel_scaling) {
  parameters <- spimalot::spim_pars_pmcmc_load(path)
  if (length(region) == 1) {
    spim_fit_pars_load_single(parameters, region, assumptions, kernel_scaling)
  } else {
    spim_fit_pars_load_nested(parameters, region, assumptions, kernel_scaling)
  }
}


spim_fit_pars_load_single <- function(path, region, assumptions,
                                      kernel_scaling) {
  info <- spim_pars_info_single(region, parameters$info)
  prior <- spim_pars_prior_single(region, info, parameters$prior)
  proposal <- spim_pars_proposal_single(region, info, parameters$proposal,
                                        kernel_scaling)

  dat <- load_transform(path, region, assumptions)
  base <- dat$base
  transform <- dat$transform

  mcmc <- spim_pars_mcmc_single(info, prior, proposal, transform)

  list(region = region,
       assumptions = assumptions,
       info = info,
       prior = prior,
       proposal = proposal,
       transform = transform,
       raw = parameters,
       base = base,
       mcmc = mcmc)
}


spim_fit_pars_load_nested <- function(path, region, assumptions,
                                      kernel_scaling) {
  info <- spim_pars_info_nested(region, parameters$info)
  prior <- spim_pars_prior_nested(region, info, parameters$prior)
  proposal <- spim_pars_proposal_nested(region, info, parameters$proposal,
                                        kernel_scaling)

  dat <- lapply(region, function(r) load_transform(path, r, assumptions))
  names(dat) <- region
  base <- lapply(dat, "[[", "base")
  transform <- lapply(dat, "[[", "transform")

  mcmc <- spim_pars_mcmc_nested(info, prior, proposal, transform)

  ret <- list(region = region,
              assumptions = assumptions,
              info = info,
              prior = prior,
              proposal = proposal,
              transform = transform,
              raw = parameters,
              base = base,
              mcmc = mcmc)
  ret
}


load_transform <- function(path, region, assumptions) {
  assert_file_exists("transform.R", path)
  assert_file_exists("base.rds", path)

  e <- new.env()
  sys.source(file.path(path, "transform.R"), e)
  stopifnot(is.function(e$make_transform),
            is.function(e$apply_assumptions))
  make_transform <- e$make_transform
  base <- e$apply_assumptions(
    readRDS(file.path(path, "base.rds"))[[region]], assumptions)
  transform <- make_transform(base)
  list(base = base, transform = transform)
}
