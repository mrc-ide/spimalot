## All the stuff here will disappear once we move entirely to the
## multi-region method.

spim_combined_load <- function(path) {
  regions <- sircovid::regions("all")

  files <- file.path(path, regions, "fit.rds")
  ## TODO: better error message if not found
  stopifnot(all(file.exists(files)))

## -  types <- c(samples = "sample_pmcmc_results.rds",
## -             pmcmc = "pmcmc_results.rds",
## -             # restart = "sample_pmcmc_restart.rds",
## -             simulate = "sample_pmcmc_simulate.rds",
## -             data = "data.rds",
## -             vaccine = "vaccine.rds",
## -             rt = "Rt.rds",
## -             ifr_t = "ifr_t.rds",
## -             deaths = "deaths_by_age.rds",
## -             admissions = "admissions_by_age.rds")


  read_rds <- function(filename, region) {
    message(sprintf("Reading fit for %s", region))
    readRDS(filename)
  }
  dat <- Map(read_rds, files, regions)
  names(dat) <- regions
  ret <- list_transpose(dat)

  message("Reordering trajectories")
  ## reorder by increasing cumulative incidence:
  rank_cum_inc <- lapply(ret$samples, sircovid::get_sample_rank)
  ret$samples <- Map(sircovid::reorder_sample, ret$samples, rank_cum_inc)
  ret$rt <- Map(sircovid::reorder_rt_ifr, ret$rt, rank_cum_inc)
  ret$ifr_t <- Map(sircovid::reorder_rt_ifr, ret$ifr_t, rank_cum_inc)

  ## NOTE: have not ported the "randomise trajectory order" bit over,
  ## but I do not think that we need to.

  message("Aggregating England/UK")
  ## Aggregate some of these to get england/uk entries
  ret$samples <- combined_aggregate_samples(ret$samples)
  ret$data <- combined_aggregate_data(ret$data)
  ret$rt <- combined_aggregate_rt(ret$rt, ret$samples)
  ret$ifr_t <- combined_aggregate_rt(ret$ifr_t, ret$samples)

  ## Copy over some core bits of data; we should check here that these
  ## all match though.
  ret$info <- ret$samples[[1]]$info[c("date", "multistrain", "model_type")]

  ret
}


combine_simulate <- function(simulate, samples, Rt_outputs, date) {

  x <- switch_levels(simulate)
  rt <- switch_levels(Rt_outputs)[c("Rt_general", "eff_Rt_general")]

  dates <- simulate[[1]]$date
  idx_dates <- samples[[1]]$trajectories$date %in% dates

  rt_combined <- lapply(rt, function(x) {
    ret <- aperm(abind::abind(x, along = 3), c(2, 3, 1))
    ret[, , idx_dates]
  })

  state <- lapply(samples, function(x) {
    x$trajectories$state[rownames(simulate[[1]]$state), , idx_dates]
  })

  state_by_age <-  switch_levels(x$state_by_age)

  ret <- list(date = dates,
              state = aperm(abind::abind(state, along = 4), c(1, 2, 4, 3)),
              state_by_age = lapply(state_by_age, abind::abind, along = 3),
              n_protected = abind::abind(x$n_protected, along = 2),
              n_doses = abind::abind(x$n_doses, along = 3))

  c(ret, rt_combined)
}


combined_aggregate_samples <- function(samples) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  samples$england$trajectories <-
    sircovid::combine_trajectories(samples[england], rank = FALSE)
  samples$uk$trajectories <-
    sircovid::combine_trajectories(samples[nations], rank = FALSE)

  samples
}


combined_aggregate_data <- function(data) {
  aggregate1 <- function(region_names, type = "full") {
    d <- lapply(data[region_names], "[[", type)
    x <- d[[1]]
    date_names <- grep("date", names(x), value = TRUE)
    data_names <- setdiff(names(x), date_names)
    dim_t <- min(sapply(d, nrow)) # TODO
    d <- lapply(d, function(x) x[seq_len(dim_t), data_names])
    ret <- Reduce( '+', d)
    death_nms <- c("deaths_hosp", "deaths_carehomes","deaths_comm")
    ret$deaths <- pmax(ret$deaths,
                       rowSums(ret[, death_nms], na.rm = TRUE),
                       na.rm = TRUE)
    ret[-seq_len(nrow(ret) - 5), ] <- NA
    data_frame(ret, x[, date_names])
  }

  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  data$england$full <- aggregate1(england, "full")
  data$england$fitted <- aggregate1(england, "fitted")
  data$england$fitted$deaths <- NA

  data$uk$full <- aggregate1(nations, "full")
  data$uk$fitted <- aggregate1(nations, "fitted")
  data$uk$fitted$deaths <- NA

  data
}


combined_aggregate_rt <- function(rt, samples) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")
  rt$england <- sircovid:::combine_rt(rt[england], samples[england],
                                      rank = FALSE)
  rt$uk <- sircovid:::combine_rt(rt[nations], samples[nations],
                                 rank = FALSE)
  rt
}


## Used to invert regions/variable in rt data
combined_switch_levels <- function(x) {
  nms <- names(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}


combined_spim_summary <- function(dat, type, what, ts_date = "2020-10-01") {
  ## date we can get
  res <- Map(region_to_template,
             Map(c, dat$samples, dat$rt),
             c(regions$name, "England", "United Kingdom"),
             MoreArgs = list(date = date,
                             type = type,
                             what = what,
                             ts_date = ts_date))
  res <- do.call(rbind, res)
  rownames(res) <- NULL

  res
}


region_to_template <- function(samples, name, date, type, what, ts_date) {
  if (what == "forecast_only") {
    ret <- rbind(
      region_to_template_forecast(samples, name, date, type))
  }
  if (what == "nowcast_forecast") {
    ret <- rbind(
      region_to_template_rt(samples, name, date, type, "nowcast", ts_date),
      region_to_template_growth_rate(samples, name, date, type, "nowcast", ts_date),
      region_to_template_incidence(samples, name, date, type, "nowcast", ts_date),
      region_to_template_prevalence(samples, name, date, type, "nowcast", ts_date),
      region_to_template_forecast(samples, name, date, type))
  }
  if (what == "time-series") {
    ret <- rbind(
      region_to_template_rt(samples, name, date, type, "time-series", ts_date),
      region_to_template_growth_rate(samples, name, date, type, "time-series", ts_date),
      region_to_template_incidence(samples, name, date, type, "time-series", ts_date),
      region_to_template_prevalence(samples, name, date, type, "time-series", ts_date))
  }

  ret
}


region_to_template_forecast <- function(samples, name, date, type) {
  rbind(
    region_to_template1(samples, name, "deaths_inc",
                        "type28_death_inc_line", date, type),
    region_to_template1(samples, name, "infections_inc",
                        "infections_inc", date, type),
    region_to_template1(samples, name, "react_pos",
                        "prevalence_mtp", date, type),
    region_to_template1(samples, name, "all_admission_inc",
                        "hospital_inc", date, type),
    region_to_template1(samples, name, "hosp",
                        "hospital_prev", date, type))
}
