##' Create summaries
##'
##' @title Create summaries
##'
##' @param dat A combined data object (loaded by
##'   [spimalot::spim_combined_load]
##'
##' @return A [data.frame] of results
##'
##' @rdname spim_summary
##' @export
spim_summary_nowcast <- function(dat) {
  message("Creating time series")
  f <- function(region) {
    message(paste("  -", region))
    rbind(
      spim_summary_region_rt(region, dat, NULL),
      spim_summary_region_growth_rate(region, dat, NULL),
      spim_summary_region_incidence(region, dat, NULL),
      spim_summary_region_prevalence(region, dat, NULL),
      spim_summary_region_forecast(region, dat))
  }

  regions <- names(dat$samples)
  res <- do.call("rbind", lapply(regions, f))
  rownames(res) <- NULL
  res
}


##' @rdname spim_summary
##' @param time_series_date The start date for the preceeding time series
##' @export
spim_summary_time_series <- function(dat, time_series_date) {
  message("Creating time series")
  f <- function(region) {
    message(paste("  -", region))
    rbind(
      spim_summary_region_rt(region, dat, time_series_date),
      spim_summary_region_growth_rate(region, dat, time_series_date),
      spim_summary_region_incidence(region, dat, time_series_date),
      spim_summary_region_prevalence(region, dat, time_series_date))
  }

  regions <- names(dat$samples)
  res <- do.call("rbind", lapply(regions, f))
  rownames(res) <- NULL
  res
}


##' @export
##' @rdname spim_summary
##'
##' @param result Result of `spim_summary_nowcast` or
##'   `spim_summary_time_series`
##'
##' @param path_template Path to the template file
##'
##' @param path_save Path to save the result at
spim_summary_write <- function(result, path_template, path_save) {
  sheets <- readxl::excel_sheets(path_template)
  template <- setNames(
    lapply(sheets, readxl::read_xlsx, path = path_template,
           .name_repair = "minimal"),
    sheets)

  stopifnot(identical(names(template[[2]]), names(result)))
  template[[2]] <- result
  writexl::write_xlsx(template, path_save)
}


##' Extract admissions by age
##'
##' @title Extract admissions by age
##'
##' @param dat Combined data set
##'
##' @export
spim_extract_admissions_by_age <- function(dat) {

  admissions_demo <- lapply(dat$samples, spim_extract_admissions_by_age_region)

  vapply(admissions_demo, "[[", numeric(19L), "mean_prop_total_admissions")

}

spim_extract_admissions_by_age_region <- function(sample) {

  trajectories <- sample$trajectories$state
  cum_admissions <- trajectories[grep("^cum_admit", rownames(trajectories)), , ]

  total_admissions <- cum_admissions[, , dim(cum_admissions)[3]]
  prop_admissions <- t(total_admissions) / colSums(total_admissions) * 100

  admissions <- apply(cum_admissions, 1:2, diff)
  mean_admissions <- apply(admissions, 1:2, mean)

  list(prop_total_admissions = prop_admissions,
       mean_prop_total_admissions = colMeans(prop_admissions),
       mean_admissions_t = mean_admissions)

}


spim_summary_region_rt <-  function(region, dat, time_series_date) {
  date <- dat$info$date
  model_type <- dat$info$model_type
  qs <- spim_summary_quantiles()

  rt <- dat$rt[[region]]$eff_Rt_general
  rt_date <- dat$rt[[region]]$date[, 1]

  ## This can be considerably simplified
  if (is.null(time_series_date)) {
    i <- which(rt_date == sircovid::sircovid_date(date))
    d <- quantile(rt[i, , drop = FALSE], probs = qs, na.rm = TRUE)
    d <- array(d, c(1, length(d)))
    sample_date <- as.Date(date)
  } else {
    i <- which(rt_date < sircovid::sircovid_date(date) &
               rt_date >= sircovid::sircovid_date(time_series_date))
    d <- t(apply((rt[i, ]), 1, quantile, qs, na.rm = TRUE, names = FALSE))
    sample_date <- sircovid::sircovid_date_as_date(rt_date[i])
  }
  colnames(d) <- names(qs)

  spim_template(
    date, model_type, "Nowcast", "R", region, sample_date, d)
}


spim_summary_region_growth_rate <-  function(region, dat, time_series_date) {
  ## Leaving the calculation part here because I am not sure if we
  ## want it elsewhere
  calculate <- function(date, trajectories, lag) {
    ## determine which dates are within the time_window
    tt <- trajectories$date
    sircovid_date <- sircovid::sircovid_date(date)
    tw <- sircovid_date - seq(lag - 1, 0)

    ## extract incidence in window of (lag) days before the data
    incid <- t(trajectories$state["sympt_cases_inc", , tt %in% tw])

    ## fit linear model to log incidence over time window
    ## to estimate r, where incid[t + 1] = incid[t] * exp(r)
    lm(log(incid + 1) ~ seq_len(lag))$coefficients[2, ]
  }

  date <- dat$info$date
  model_type <- dat$info$model_type
  qs <- spim_summary_quantiles()

  trajectories <- dat$samples[[region]]$trajectories

  if (is.null(time_series_date)) {
    growth_rate <- calculate(date, trajectories, lag = 10)
    d <- signif(quantile(growth_rate, qs), 3)
    d <- array(d, c(1, length(d)))
    sample_date <- date
  } else {
    sample_date <- as.character(seq.Date(as.Date(time_series_date),
                                         as.Date(date) - 1, 1))
    growth_rate <- t(vapply(sample_date,
                            function(x) calculate(x, trajectories, lag = 10),
                            numeric(ncol(trajectories$state))))
    d <- signif(t(apply(growth_rate, 1, quantile, qs, names = FALSE)), 3)
  }

  colnames(d) <- names(qs)

  spim_template(
    date, model_type, "Nowcast", "growth_rate", region, sample_date, d)
}


spim_summary_region_incidence <-  function(region, dat, time_series_date) {
  date <- dat$info$date
  model_type <- dat$info$model_type
  qs <- spim_summary_quantiles()

  trajectories <- dat$samples[[region]]$trajectories

  incid <- t(trajectories$state["infections_inc", , ])

  if (is.null(time_series_date)) {
    i <- which(trajectories$date == sircovid::sircovid_date(date))
    d <- signif(quantile(incid[i, ], qs), 3)
    d <- array(d, c(1, length(d)))
    sample_date <- as.Date(date)
  } else {
    i <- which(trajectories$date < sircovid::sircovid_date(date) &
               trajectories$date >= sircovid::sircovid_date(time_series_date))
    d <- signif(t(apply(incid[i, ], 1, quantile, qs, names = FALSE)), 3)
    sample_date <- sircovid::sircovid_date_as_date(trajectories$date[i])
  }

  colnames(d) <- names(qs)

  spim_template(
    date, model_type, "Nowcast", "incidence", region, sample_date, d)
}


spim_summary_region_prevalence <-  function(region, dat, time_series_date) {
  date <- dat$info$date
  model_type <- dat$info$model_type
  qs <- spim_summary_quantiles()

  trajectories <- dat$samples[[region]]$trajectories

  prev <- t(trajectories$state["react_pos", , ])

  if (is.null(time_series_date)) {
    i <- which(trajectories$date == sircovid::sircovid_date(date))
    d <- signif(quantile(prev[i, ], qs), 3)
    d <- array(d, c(1, length(d)))
    sample_date <- as.Date(date)
  } else {
    i <- which(trajectories$date < sircovid::sircovid_date(date) &
               trajectories$date >= sircovid::sircovid_date(time_series_date))
    d <- signif(t(apply(prev[i, ], 1, quantile, qs, names = FALSE)), 3)
    sample_date <- sircovid::sircovid_date_as_date(trajectories$date[i])
  }

  colnames(d) <- names(qs)

  ## Scale by total population
  pop_size <- sum(sircovid::carehomes_parameters(1, region)$N_tot[2:18])
  d <- d / pop_size * 100

  spim_template(
    date, model_type, "Nowcast", "prevalence", region, sample_date, d)
}


spim_summary_region_forecast <- function(region, dat) {
  rbind(
    spim_summary_region_forecast_trajectory(
      region, dat, "deaths_inc", "type28_death_inc_line"),
    spim_summary_region_forecast_trajectory(
      region, dat, "infections_inc", "infections_inc"),
    spim_summary_region_forecast_trajectory(
      region, dat, "react_pos", "prevalence_mtp"),
    spim_summary_region_forecast_trajectory(
      region, dat, "all_admission_inc", "hospital_inc"),
    spim_summary_region_forecast_trajectory(
      region, dat, "hosp", "hospital_prev"))
}


spim_summary_region_forecast_trajectory <- function(region, date, name, as) {
  date <- dat$info$date
  model_type <- dat$info$model_type
  qs <- spim_summary_quantiles()

  trajectories <- dat$samples[[region]]$trajectories

  if (name == "deaths_inc") {
    value <-
      trajectories$state["deaths_hosp_inc", , ] +
      trajectories$state["deaths_comm_inc", , ] +
      trajectories$state["deaths_carehomes_inc", , ]
  } else if (name == "all_admission_inc") {
    value <-
      trajectories$state["diagnoses_inc", , ] +
      trajectories$state["admitted_inc", , ]
  } else {
    value <- trajectories$state[name, , ]
  }

  i <- which(trajectories$predicted)
  d <- round(t(apply(t(value[, i]), 1, quantile, qs, names = FALSE)))
  colnames(d) <- names(qs)

  sample_date <- sircovid::sircovid_date_as_date(trajectories$date[i])

  if (name == "react_pos") {
    pop_size <- sum(sircovid::carehomes_parameters(1, region)$N_tot[2:18])
    d <- d / pop_size * 100
  }

  spim_template(
    date, model_type, "MTP", as, region, sample_date, d)
}


spim_template <- function(date, model_type, scenario, value_type, region,
                          sample_date, data) {
  version <- packageVersion("sircovid")
  version <- substr(version, 3, nchar(version))
  if (model_type == "BB") {
    model <- "Stochastic Compartmental Positivity"
  } else {
    model <- "Stochastic Compartmental Cases"
  }

  ## TODO: translate region into the correct value (london => London
  ## etc)

  ret <- data_frame(Group = "Imperial",
                    Model = model,
                    Scenario = scenario,
                    ModelType = "Multiple",
                    Version = version,
                    "Creation Day" = lubridate::day(date),
                    "Creation Month" = lubridate::month(date),
                    "Creation Year" = lubridate::year(date),
                    "Day of Value" = lubridate::day(sample_date),
                    "Month of Value" = lubridate::month(sample_date),
                    "Year of Value" = lubridate::year(sample_date),
                    AgeBand = "All",
                    Geography = spim_region_name(region),
                    ValueType = value_type,
                    data)
  rownames(ret) <- NULL
  ret
}


spim_summary_quantiles <- function() {
  qs <- seq(from = 0.05, to = 0.95, by = 0.05)
  names(qs) <- sprintf("Quantile %s", qs)
  c(c(Value = 0.5), qs)
}
