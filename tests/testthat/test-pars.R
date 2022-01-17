test_that("spim_pars_load throws useful error messages", {
  dat <- test_dummy_pars_pmcmc(c("a", "b", "c"))
  path <- tempfile()
  expect_error(
    spim_pars_pmcmc_load(path),
    sprintf("File '.+' does not exist", path))

  dir.create(path, FALSE, TRUE)
  expect_error(
    spim_pars_pmcmc_load(path),
    "File 'info.csv' does not exist")
  write_csv(dat$info, file.path(path, "info.csv"))
  expect_error(
    spim_pars_pmcmc_load(path),
    "File 'prior.csv' does not exist")
  write_csv(dat$prior, file.path(path, "prior.csv"))
  expect_error(
    spim_pars_pmcmc_load(path),
    "File 'proposal.csv' does not exist")
  write_csv(dat$proposal, file.path(path, "proposal.csv"))

  pars <- spim_pars_pmcmc_load(path)

  expect_s3_class(pars, "spim_pars_pmcmc")
  expect_equal(pars$info, dat$info)
  expect_equal(pars$prior, dat$prior)
  expect_equal(pars$proposal, dat$proposal)

  path_out <- tempfile()
  spim_pars_pmcmc_save(pars, path_out)
  expect_setequal(
    dir(path_out),
    c("info.csv", "prior.csv", "proposal.csv"))
})


test_that("filter info by region", {
  dat <- test_dummy_pars_pmcmc(c("a", "b", "c"))
  info <- spim_pars_info("london", dat$info)
  expect_equal(nrow(info), 3)
  expect_setequal(names(info), c("name", "initial", "min", "max", "discrete"))
  expect_error(
    spim_pars_info("islington", dat$info),
    "Did not find region 'islington' in parameter info")
  expect_error(
    spim_pars_info("london", dat$info[-3]),
    "Required names missing from 'info'")
})


test_that("filter prior by region", {
  dat <- test_dummy_pars_pmcmc(c("a", "b", "c"))
  info <- spim_pars_info("london", dat$info)
  prior <- spim_pars_prior("london", info, dat$prior)
  expect_equal(nrow(prior), 3)
  expect_setequal(
    names(prior),
    c("type", "name", "gamma_scale", "gamma_shape", "beta_shape1",
      "beta_shape2"))
  expect_error(
    spim_pars_prior("london", info, dat$prior[dat$prior$name != "b", ]),
    "'prior$name' and 'info$name' are not setequal",
    fixed = TRUE)
})


test_that("filter proposal by region", {
  dat <- test_dummy_pars_pmcmc(c("a", "b", "c"))
  info <- spim_pars_info("london", dat$info)
  proposal <- spim_pars_proposal("london", info, dat$proposal, 1)
  expect_equal(spim_pars_proposal("london", info, dat$proposal, 2),
               proposal * 2)

  expect_equal(nrow(proposal), 3)
  expect_equal(ncol(proposal), 3)
  expect_equal(colnames(proposal), c("a", "b", "c"))
  expect_equal(rownames(proposal), colnames(proposal))

  expect_error(
    spim_pars_proposal("london", info,
                       dat$proposal[dat$proposal$name != "b", ]),
    "'proposal$name' and 'info$name' are not setequal",
    fixed = TRUE)
  expect_error(
    spim_pars_proposal("london", info,
                       dat$proposal[names(dat$proposal) != "b"]),
    "are not setequal",
    fixed = TRUE)
})


test_that("make priors", {
  d <- list(gamma_shape = 2,
            gamma_scale = 3,
            beta_shape1 = 0.4,
            beta_shape2 = 0.7)
  p <- make_prior(modifyList(d, list(type = "gamma")))
  expect_equal(p(0.5), dgamma(0.5, shape = 2, scale = 3, log = TRUE))

  p <- make_prior(modifyList(d, list(type = "beta")))
  expect_equal(p(0.5), dbeta(0.5, shape1 = 0.4, shape2 = 0.7, log = TRUE))

  expect_null(make_prior(modifyList(d, list(type = "null"))))

  expect_error(
    make_prior(modifyList(d, list(type = "normal"))),
    "Unknown prior type")
})
