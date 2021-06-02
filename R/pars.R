## Our main parameter function is a bit of a beast because there are
## many many things that go into it. It's not totally obvious that
## this can be greatly simplified in the interface as these really are
## separate bits of inputs.

##' Create parameters object for use with [mcstate::pmcmc]. This
##' function is very data-hungry as it (alongside
##' [spimalot::spim_data()]) is the main point at which data enters
##' into the model. Here we end up setting a lot of data that are are
##' not fitting to, but taking as fixed inputs (vaccination over time)
##' as well as information on the parameters that the model will vary
##' (their names, ranges, priors and their proposal distributions).
##'
##' @title Create parameters object
##'
##' @inheritParams spim_data
##'
##' @param multistrain Logical, indicating if the model is a
##'   "multistrain" model allowing for mulitiple competing strains.
##'
##' @param beta_date A vector of date (strings) for the beta
##'   parameters. Must align with parameters
##'
##' @param vaccination Vacination data, from
##'   [spimalot::spim_vaccination_data]
##'
##' @param parameters Parameter information, from
##'   [spimalot::spim_pars_pmcmc_load]
##'
##' @return An [mcstate::pmcmc_parameters] object which can be used
##'   with [mcstate::pmcmc]
##'
##' @export
spim_pars <- function(date, region, model_type, multistrain,
                      beta_date, vaccination, parameters) {
  assert_is(parameters, "spim_pars_pmcmc")

  if (length(region) == 1) {
    ret <- spim_pars_single(date, region, model_type, multistrain,
                            beta_date, vaccination, parameters)
  } else {
    stop("writeme")
  }

  ## This will allow us to recreate things later in the restart
  inputs <- list(date = date,
                 region = region,
                 model_type = model_type,
                 multistrain = multistrain,
                 beta_date = beta_date,
                 vaccination = vaccination,
                 parameters = parameters)

  attr(ret, "inputs") <- inputs

  ret
}


##' Load the pmcmc parameters from disk. We expect three files; one
##' for the overal parameters (`info`), one with details of the priors
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


##' Write pmcmc paramters out after fitting and aggregation (inverse
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
  stopifnot(setequal(names(p), c("info", "prior", "proposal")))
  dir.create(path, FALSE, TRUE)
  write_csv(p$info, file.path(path, "info.csv"))
  write_csv(p$prior, file.path(path, "prior.csv"))
  write_csv(p$proposal, file.path(path, "proposal.csv"))
}


spim_pars_single <- function(date, region, model_type, multistrain,
                             beta_date, vaccination, parameters) {
  ## We take 'info' as the canonical source of names, then check that
  ## prior and proposal align correctly.
  info <- spim_pars_info(region, parameters$info)
  prior <- spim_pars_prior(region, info, parameters$prior)
  proposal <- spim_pars_proposal(region, info, parameters$proposal)

  pars <- Map(
    mcstate::pmcmc_parameter,
    name = info$name,
    initial = info$initial,
    min = info$min,
    max = info$max,
    discrete = info$discrete,
    prior = lapply(split(prior, prior$name), make_prior))

  transform <- spim_transform(region, model_type, multistrain, beta_date,
                              vaccination)

  mcstate::pmcmc_parameters$new(pars, proposal, transform)
}


##' Create vector of beta dates
##'
##' @title Create vector of beta dates
##'
##' @param date Todays dates (last date is set two weeks ago)
##'
##' @export
spim_pars_beta <- function(date) {
  ## Dates are as follows
  ##  1. 2020-03-16 - PM advises WFH, against non-essential travel etc
  ##  2. 2020-03-23 - PM announces full lockdown
  ##  3. 2020-03-25 - lockdown into full effect
  ##  4. 2020-05-11 - initial easing of lockdown
  ##  5. 2020-06-15 - non-essential shops can open
  ##  6. 2020-07-04 - restaurants, pubs etc can open
  ##  7. 2020-08-01 - "Eat out to help out" scheme starts
  ##  8. 2020-09-01 - Schools reopen
  ##  9. 2020-09-14 - "Rule of six" introduced
  ## 10. 2020-10-14 - Tiered system introduced
  ## 11. 2020-10-31 - lockdown announced
  ## 12. 2020-11-05 - lockdown 2 starts
  ## 13. 2020-12-02 - lockdown 2 ends
  ## 14. 2020-12-18 - school Christmas holidays
  ## 15. 2021-01-05 - Lockdown 3 starts
  ## 16. 2021-03-08 - Step 1 of roadmap: schools reopen
  ## 17. 2021-04-12 - Step 2 of roadmap: outdoors hospitality
  ## 18. Two weeks ago
  ## TODO: recent dates to consider: schools return (2021-04-19)
  c("2020-03-16", "2020-03-23", "2020-03-25",
    "2020-05-11", "2020-06-15", "2020-07-04",
    "2020-08-01", "2020-09-01", "2020-09-14",
    "2020-10-14", "2020-10-31", "2020-11-05",
    "2020-12-02", "2020-12-18", "2021-01-05",
    "2021-03-08", "2021-04-12",
    as.character(as.Date(date) - 21))
}


spim_pars_pmcmc <- function(region, info, prior, proposal) {
  if (length(region) == 1) {
    spim_pars_pmcmc_single(region, info, prior, proposal)
  } else {
    stop("writeme")
  }
}


spim_pars_pmcmc_single <- function(region, info, prior, proposal) {

  list(region = region, info = info, prior = prior, proposal = proposal)
}


spim_pars_info <- function(region, info) {
  assert_has_names(info, c("name", "initial", "max", "discrete"))
  info <- info[info$region == region & info$include, ]
  info <- info[setdiff(names(info), c("region", "include"))]
  rownames(info) <- NULL
  assert_unique(info$names)
  info
}


spim_pars_prior <- function(region, info, prior) {
  prior_cols <- c("region", "type", "name", "gamma_scale", "gamma_shape",
                  "beta_shape1", "beta_shape2")
  assert_has_names(prior, prior_cols)
  prior <- prior[prior$region == region & prior$name %in% info$name, ]
  prior <- prior[setdiff(names(prior), "region")]
  rownames(prior) <- NULL
  assert_setequal(prior$name, info$name)
  prior
}


spim_pars_proposal <- function(region, info, proposal) {
  proposal <- proposal[proposal$region == region, ]
  assert_setequal(proposal$name, info$name)
  assert_setequal(setdiff(names(proposal), c("name", "region")), info$name)
  proposal <- as.matrix(proposal[match(info$name, proposal$name), info$name])
  rownames(proposal) <- info$name
  proposal
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


parameter_subset_region <- function(pars, region) {
  assert_is(pars, "spim_pars_pmcmc")
  lapply(pars, function(x) x[x$region == region, ])
}
