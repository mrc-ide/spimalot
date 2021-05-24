fit1 <- function(region, path) {
  message(sprintf("**** %s ****", region))

  date <- "2021-05-17"
  model_type <- "BB"
  multistrain <- FALSE
  short_run <- TRUE
  n_chains <- 3

  trim_deaths <- 4
  forecast_days <- 57

  control <- spimalot::spim_control(short_run, n_chains)

  ## TODO: this is still not quite what we will need for the restart
  ## because that needs to be filtered by region; we're going to hit
  ## this later too
  parameters <- spimalot::spim_pars_pmcmc_load("example/bb")

  data_vaccination <- read_csv("example/data_vaccination.csv")
  data_rtm <- read_csv("example/data_rtm.csv")
  data_serology <- read_csv("example/data_serology.csv")

  ## This is new, and used only in sorting out the final outputs. Some
  ## explanation would be useful.
  data_admissions <- read_csv("example/admissions_age_sitrep.csv")
  data_admissions <- spim_data_admissions(data_admissions, region)

  beta_date <- spimalot::spim_pars_beta(date)

  ## Vaccination "parameters" inputs; we'll need some real tidyup around
  ## this soon.
  efficacy <- readRDS("example/vaccine_efficacy.rds")
  uptake <- c(0, 0, 0, 0.32, 0.8, 0.8, 0.9, 0.9, 0.9, 0.9, 0.95, 0.95, 0.95,
              0.95, 0.95, 0.95, 0.95, 0.85, 0.95)
  mean_days_between_doses <- 77
  end_date <- as.Date(date) + forecast_days
  set.seed(1) # this has stochastic tiebreaks
  vaccination <- spimalot::spim_vaccination_data(
    date, region, uptake, end_date,
    mean_days_between_doses, efficacy,
    data_vaccination)

  pars <- spimalot::spim_pars(date, region, model_type, multistrain,
                              beta_date, vaccination, parameters)

  data <- spimalot::spim_data(date, region, model_type, data_rtm, data_serology,
                              trim_deaths, full_data = FALSE)
  ## This is the data set including series that we do not fit to, and
  ## with the full series of carehomes deaths.
  data_full <- spimalot::spim_data(date, region, model_type, data_rtm,
                                   data_serology, trim_deaths, full_data = TRUE)

  filter <- spimalot::spim_particle_filter(data, pars, control$particle_filter)
  samples <- spimalot::spim_pmcmc(pars, filter, control$pmcmc)

  ## To run the model at this point, we just need to run:
  ##
  ## > filter$run(pars$model(pars$initial()))

  dat <- spimalot::spim_fit_process(samples, control$forecast,
                                    data_admissions, data_rtm)

  ## One more thing; added here rather than in the processing because
  ## otherwise the arg list is a bit weird with all the odd bits of
  ## data being passed in. But this gives us the data set we fitted to
  ## along with a full version.
  dat$data <- list(fitted = data,
                   full = data_full)
  dat$vacination <- vaccination

  dest <- file.path(path, region)
  dir.create(dest, FALSE, TRUE)

  saveRDS(dat, file.path(dest, "fit.rds"))
}

pkgload::load_all()
for (region in sircovid::regions("all")) {
  fit1(region, "example/fits")
}
