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
##' @inheritParameters spim_data
##'
##' @param multistrain Logical, indicating if the model is a
##'   "multistrain" model allowing for mulitiple competing strains.
##'
##' @param vaccination Vacination data (TODO: DESCRIBE FORMAT, TODO:
##'   DESCRIBE ORIGIN)
##'
##' @param info Parameter information, describing the parameters that
##'   will be varied and their ranges. TODO: DESCRIBE FORMAT
##'
##' @param proposal Parameter proposal information, used to create the
##'   variance covariance matrix for the mcmc sampler. TODO: DESCRIBE
##'   FORMAT
##'
##' @param prior Parameter prior information TODO: DESCRIBE FORMAT
##'
##' @return An [mcstate::pmcmc_parameters] object which can be used
##'   with [mcstate::pmcmc]
##'
##' @export
spim_pars <- function(date, region, model_type, multistrain,
                      vaccination, info, prior, proposal) {
  severity <- read_csv(spimalot_file("extdata/support_severity.csv"))
  progression <- read_csv(spimalot_file("extdata/support_progression.csv"))
  beta <- spim_pars_beta(date)

  if (length(region) == 1) {
    spim_pars_single(date, region, model_type, multistrain,
                     beta, severity, progression, vaccination,
                     info, prior, proposal)
  } else {
    stop("writeme")
  }
}


spim_pars_single <- function(date, region, model_type, multistrain,
                             beta, severity, progression, vaccination,
                             info, prior, proposal) {
  ## We take 'info' as the canonical source of names, then check that
  ## prior and proposal align correctly.
  info <- spim_pars_info(region, info)
  prior <- spim_pars_prior(region, info, prior)
  proposal <- spim_pars_prior(region, info, proposal)

  browser()

  pars <- Map(
    mcstate::pmcmc_parameter,
    name = info$name,
    initial = info$initial,
    min = info$min,
    max = info$max,
    discrete = info$discrete,
    prior = lapply(split(prior, prior$name), make_prior))

  transform <- spim_transform(region, model_type, beta_date,
                              severity, progression,
                              vaccination, multistrain)

  mcstate::pmcmc_parameters$new(pars, proposal, transform)


  browser()
}



## TODO: a version of this which allows tweaking around recent dates
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
  ## 16. 2021-03-08 - Step 1 of lockdown end roadmap: schools reopen
  ## 17. Two weeks ago
  ## TODO: recent dates to consider: easter hols (2021-04-01), step 2 roadmap
  ## (2021-04-12) and schools return (2021-04-19)
  c("2020-03-16", "2020-03-23", "2020-03-25",
    "2020-05-11", "2020-06-15", "2020-07-04",
    "2020-08-01", "2020-09-01", "2020-09-14",
    "2020-10-14", "2020-10-31", "2020-11-05",
    "2020-12-02", "2020-12-18", "2021-01-05",
    "2021-03-08", as.character(as.Date(date) - 21))
}



## spim_pars_single <- function(region, ...) {
##   stopifnot(
##     identical(parameters$info$name, parameters$prior$name),
##     identical(parameters$info$name, rownames(parameters$proposal)))

##   prior <- lapply(split(parameters$prior, parameters$prior$name),
##                   make_prior)

##   pars <- Map(
##     mcstate::pmcmc_parameter,
##     name = parameters$info$name,
##     initial = parameters$info$initial,
##     min = parameters$info$min,
##     max = parameters$info$max,
##     discrete = parameters$info$discrete,
##     prior = prior)

##   transform <- carehomes_spim_transform(type, beta_date, region,
##                                         severity,
##                                         progression,
##                                         support_vaccine,
##                                         multistrain)

##   mcstate::pmcmc_parameters$new(pars, parameters$proposal, transform)
## }


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
