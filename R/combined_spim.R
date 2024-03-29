##' Create summaries
##'
##' @title Create summaries
##'
##' @param dat A combined data object (loaded by
##'   [spimalot::spim_combined_load]
##'
##' @param sircovid_model The name of the sircovid model used.
##'   Default is `"carehomes"`
##'
##' @return A [data.frame] of results
##'
##' @rdname spim_summary
##' @export
spim_summary_nowcast <- function(dat, sircovid_model = "carehomes") {
  spim_check_sircovid_model(sircovid_model)
  message("Creating time series")
  f <- function(region) {
    message(paste("  -", region))
    rbind(
      spim_summary_region_rt(region, dat, NULL),
      spim_summary_region_growth_rate(region, dat, NULL),
      spim_summary_region_incidence(region, dat, NULL),
      spim_summary_region_prevalence(region, dat, NULL, sircovid_model),
      spim_summary_region_forecast(region, dat, sircovid_model))
  }

  regions <- names(dat$samples)
  res <- do.call("rbind", lapply(regions, f))
  rownames(res) <- NULL
  res
}


##' @rdname spim_summary
##' @param time_series_date The start date for the preceding time series
##'
##' @param sircovid_model The name of the sircovid model used.
##'   Default is `"carehomes"`
##'
##' @export
spim_summary_time_series <- function(dat, time_series_date,
                                     sircovid_model = "carehomes") {
  message("Creating time series")
  f <- function(region) {
    message(paste("  -", region))
    rbind(
      spim_summary_region_rt(region, dat, time_series_date),
      spim_summary_region_growth_rate(region, dat, time_series_date),
      spim_summary_region_incidence(region, dat, time_series_date),
      spim_summary_region_prevalence(region, dat, time_series_date,
                                     sircovid_model))
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


##' Extract current Rt by variant
##'
##' @title Extract multivariant Rt
##'
##' @param dat Combined data set
##' @param type_rt Either 'Rt_general' or 'eff_Rt_general'
##'
##' @export
spim_extract_variants_rt <- function(dat, type_rt) {

  stopifnot(type_rt %in% c("Rt_general", "eff_Rt_general"))

  rt <- dat$variant_rt

  date_rt <- dat$info$date

  out <- NULL
  for (i in c(sircovid::regions("england"), "england")) {
    curr_rt <- rt[[i]][[type_rt]][
      rt[[i]]$date == sircovid::sircovid_date(date_rt)]
    curr_rt_alpha <- curr_rt[c(TRUE, FALSE)]
    curr_rt_delta <- curr_rt[c(FALSE, TRUE)]
    alpha <- c(mean(curr_rt_alpha), quantile(curr_rt_alpha, c(0.025, 0.975)))
    delta <- c(mean(curr_rt_delta), quantile(curr_rt_delta, c(0.025, 0.975)))
    epsilon <- delta / alpha
    ret <- rbind(alpha, delta, epsilon)
    colnames(ret) <- c("mean", "LB", "UB")
    ret <- as.data.frame(ret)
    ret$variant <- c("alpha", "delta", "epsilon")
    ret$region <- i

    out <- rbind(out, ret)
  }
  out
}


##' Extract current ALOS by region
##'
##' @title Extract average length of hospital stay
##'
##' @param dat Combined data set
##' @param regions String or vector of strings with region names for output
##' @param date Latest data point date
##'
##' @export
spim_extract_alos <- function(dat, regions, date) {

  tmp <- dat$ifr_t

  out <- NULL
  for (i in regions) {
    curr <- tmp[[i]]$ALOS[tmp[[i]]$date == sircovid::sircovid_date(date)]
    curr <- c(i, mean(curr), quantile(curr, c(0.025, 0.975)))
    out <- rbind(out, curr)
  }
  out
}



spim_summary_region_rt <-  function(region, dat, time_series_date) {
  date <- dat$info$date
  model_type <- dat$info$model_type
  qs <- spim_summary_quantiles()

  rt <- dat$rt[[region]]$eff_Rt_general[, "weighted", ]
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


spim_summary_region_prevalence <-  function(region, dat, time_series_date,
                                            sircovid_model) {
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
  if (sircovid_model == "lancelot") {
    pop_size <- sum(sircovid::lancelot_parameters(1, region)$N_tot[2:18])
  }

  d <- d / pop_size * 100

  spim_template(
    date, model_type, "Nowcast", "prevalence", region, sample_date, d)
}


spim_summary_region_forecast <- function(region, dat, sircovid_model) {
  rbind(
    spim_summary_region_forecast_trajectory(
      region, dat, "deaths_inc", "type28_death_inc_line", sircovid_model),
    spim_summary_region_forecast_trajectory(
      region, dat, "infections_inc", "infections_inc", sircovid_model),
    spim_summary_region_forecast_trajectory(
      region, dat, "react_pos", "prevalence_mtp", sircovid_model),
    spim_summary_region_forecast_trajectory(
      region, dat, "all_admission_inc", "hospital_inc", sircovid_model),
    spim_summary_region_forecast_trajectory(
      region, dat, "hosp", "hospital_prev", sircovid_model))
}


spim_summary_region_forecast_trajectory <- function(region, dat, name, as,
                                                    sircovid_model) {
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
    if (sircovid_model == "lancelot") {
      pop_size <- sum(sircovid::lancelot_parameters(1, region)$N_tot[2:18])
    }
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
