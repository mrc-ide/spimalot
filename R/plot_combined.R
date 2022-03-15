##' Plot trajectories
##'
##' @title Plot trajectories
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param what Vector of traces to add
##'
##' @param date_min Starting date for plot. If not given, defaults to
##'   the beginning of the data.
##'
##' @param age_band Age band to plot
##'
##' @param with_forecast Logical, indicating if we should add the forecast
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_trajectories <- function(dat, regions, what, date_min = NULL,
                                   age_band = NULL, with_forecast = TRUE,
                                   add_betas = FALSE) {

  n_regions <- length(regions)
  op <- par(mfcol = c(length(what), n_regions),
            oma = c(1, 1, 4, 1),
            mar = c(3, 3, 0.5, 0.5))
  on.exit(par(op))

  for (r in regions) {
    spim_plot_trajectories_region(r, dat, what, date_min, age_band,
                                  with_forecast = with_forecast,
                                  add_betas = add_betas)
  }

  mtext(side = 3, text = toupper(spim_region_name(regions)),
        line = 0.5, outer = TRUE, cex = 0.8,
        at = seq(1 / n_regions / 2, by = 1 / n_regions, length.out = n_regions))
}


##' Plot trajectories by age
##'
##' @title Plot trajectories by age
##'
##' @param dat Combined data set
##'
##' @param what Trajectory to plot
##'
##' @param regions Vector of regions to plot
##'
##' @param date_min Starting date for plot
##'
##' @param age_band Age band to plot
##'
##' @param with_forecast Logical, indicating if we should add the forecast
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_trajectories_by_age <- function(dat, regions, what, date_min = NULL,
                                          age_band = NULL, with_forecast = TRUE,
                                          add_betas = FALSE) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    spim_plot_trajectories_region1(
      what, r, dat, date_min, age_band,
      with_forecast = with_forecast, add_betas = add_betas,
      main = toupper(spim_region_name(r)))
  }
}


##' Plot Rt
##'
##' @title Plot Rt
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param rt_type One of the valid Rt types (e.g., eff_Rt_all)
##'
##' @param forecast_until Optional date to forecast till
##'
##' @param variant Variant type (or `NULL`)
##'
##' @param add_betas Logical, indicating if betas should be added to the plot
##'
##' @param multistrain Logical, indicating if this is a multistrain fit
##'
##' @export
spim_plot_Rt <- function(dat, regions, rt_type, forecast_until = NULL,
                         variant = NULL, add_betas = FALSE,
                         multistrain = TRUE) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))
  if (is.null(forecast_until)) {
    forecast_until <- dat$info$date
  }
  if (is.null(variant)) {
    variant <- "weighted"
  }
  for (r in regions) {
    spim_plot_Rt_region(r, dat, rt_type, forecast_until, variant, add_betas,
                        multistrain)
  }
}

##' Plot IFR over time
##'
##' @title Plot IFR over time
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param ifr_t_type One of the valid IFR_t types (e.g., IFR_t_all)
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @param forecast_until Optional date to forecast till
##'
##' @export
spim_plot_ifr_t <- function(dat, regions, ifr_t_type, ymax = 2.5,
                            forecast_until = NULL) {
  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))
  if (is.null(forecast_until)) {
    forecast_until <- dat$info$date
  }
  for (r in regions) {
    spim_plot_ifr_t_region(r, dat, ifr_t_type, ymax, forecast_until)
  }
}


##' Plot average length of stay over time
##'
##' @title Plot average length of stay over time
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param ymin Minimum percentage on y-axis
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @param forecast_until Optional date to forecast till
##'
##' @export
spim_plot_alos <- function(dat, regions, ymin = NULL, ymax = NULL,
                           forecast_until = NULL) {
  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))
  if (is.null(forecast_until)) {
    forecast_until <- dat$info$date
  }
  for (r in regions) {
    spim_plot_alos_region(r, dat, ymin, ymax, forecast_until)
  }
}


##' Plot serology
##'
##' @title Plot serology
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param sero_flow Number identifying which serology flow to use (1 or 2)
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @export
spim_plot_serology <- function(dat, regions, sero_flow, ymax) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_serology_region(r, dat, sero_flow, ymax, TRUE)
    } else {
      spim_plot_serology_region(r, dat, sero_flow, ymax)
    }
  }
}


##' Plot effective susceptibles
##'
##' @title Plot effective susceptibles
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param strain_names Names of strains
##'
##' @param strain_dates Dates of strains introduced into model
##'
##' @export
spim_plot_effective_susceptible <- function(dat, regions,
                                            strain_names, strain_dates) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_effective_susceptible_region(r, dat, strain_names,
                                             strain_dates, TRUE)
    } else {
      spim_plot_effective_susceptible_region(r, dat, strain_names, strain_dates)
    }
  }
}


##' Plot infections per strain
##'
##' @title Plot infections per strain
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param strain_names Names of strains
##'
##' @param strain_dates Dates of strains introduced into model
##'
##' @export
spim_plot_infections_per_strain <- function(dat, regions,
                                            strain_names, strain_dates) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_infections_per_strain_region(r, dat, strain_names,
                                             strain_dates, TRUE)
    } else {
      spim_plot_infections_per_strain_region(r, dat, strain_names, strain_dates)
    }
  }
}


##' Plot vaccine status
##'
##' @title Plot vaccine status
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param vacc_names Character vector of names of vaccine statuses
##'
##' @export
spim_plot_vaccine_status <- function(dat, regions, vacc_names) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_vaccine_status_region(r, dat, vacc_names, TRUE)
    } else {
      spim_plot_vaccine_status_region(r, dat, vacc_names)
    }
  }
}


##' Plot infection status
##'
##' @title Plot infection status
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @export
spim_plot_infection_status <- function(dat, regions) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_infection_status_region(r, dat, TRUE)
    } else {
      spim_plot_infection_status_region(r, dat)
    }
  }
}


##' Plot cumulative attack rate
##'
##' @title Plot cumulative attack rate
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @export
spim_plot_cumulative_attack_rate <- function(dat, regions) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_cumulative_attack_rate_region(r, dat, TRUE)
    } else {
      spim_plot_cumulative_attack_rate_region(r, dat)
    }
  }

}


##' Plot Pillar 2 positivity
##'
##' @title Plot Pillar 2 positivity
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param age_band Age band to plot
##'
##' @param date_min Starting date for plot
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_pillar2_positivity <- function(dat, regions, age_band, date_min,
                                         ymax, add_betas = FALSE) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    spim_plot_pillar2_positivity_region(r, dat, age_band,
                                        date_min, ymax, add_betas)
  }
}


##' Plot Pillar 2 cases
##'
##' @title Plot Pillar 2 cases
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param age_band Age band to plot
##'
##' @param date_min Starting date for plot
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_pillar2_cases <- function(dat, regions, age_band, date_min,
                                    add_betas = FALSE) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    spim_plot_pillar2_cases_region(r, dat, age_band, date_min, add_betas)
  }
}


##' Plot REACT
##'
##' @title Plot REACT
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param date_min Starting date for plot
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @param age_band String, indicating which age band to plot
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_react <- function(dat, regions, date_min, ymax, age_band = NULL,
                            add_betas = FALSE) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    spim_plot_react_region(r, dat, date_min, ymax, age_band, add_betas)
  }
}


##' Plot incidence
##'
##' @title Plot incidence
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param per_1000 Logical, indicating if we plot incidence per 1000 people
##'
##' @export
spim_plot_incidence <- function(dat, regions, per_1000 = FALSE) {

  oo <- par(mfrow = c(2, length(regions) / 2), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    spim_plot_incidence_region(r, dat, per_1000)
  }
}


##' Plot proportion susceptible
##'
##' @title Plot proportion susceptible
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param ymin Minimum y axis value
##'
##' @export
spim_plot_prop_susceptible <- function(dat, regions, ymin) {

  oo <- par(mfrow = c(2, 6), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    if (r == regions[length(regions)]) {
      spim_plot_prop_susceptible_region(r, dat, ymin, TRUE)
    } else {
      spim_plot_prop_susceptible_region(r, dat, ymin)
    }
  }

}

##' Plot log trajectories by age
##'
##' @title Plot log trajectories by age
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param yield Which kind of age-specific trajectories to add, can be
##'   "pillar2", "admissions" or "deaths"
##'
##' @param model_type Which model is being run, either "BB" or "NB"
##'
##' @export
spim_plot_log_traj_by_age <- function(dat, regions, yield, model_type) {

  if (yield == "pillar2") {
    what <- c("under15", "15_24", "25_49", "50_64", "65_79", "80_plus")
    oo <- par(mfrow = c(length(regions), length(what)), oma = c(1, 1, 4, 1),
              mar = c(3, 3, 0.5, 0.5))
    on.exit(oo)
    for (r in regions) {
      spim_plot_log_traj_by_age_region(r, dat, yield, what, model_type)
    }
    mtext(side = 3,
          text = stringr::str_to_sentence(stringr::str_replace(what, "_", " ")),
          line = 0.5, outer = TRUE, cex = 0.9,
          at = c(0.1, 0.26, 0.43, 0.6, 0.76, 0.93))
  } else {
    oo <- par(mfrow = c(length(regions), 5), oma = c(1, 1, 4, 1),
              mar = c(3, 3, 0.5, 0.5))
    on.exit(oo)
    if (yield == "admissions") {
      what <- c("adm_0", "adm_25", "adm_55", "adm_65", "adm_75")
      what_names <- c("Admissions under 25s", "Admissions 25 to 54",
                      "Admissions 55 to 64", "Admissions 65 to 74",
                      "Admissions 75+")
    }
    if (yield == "deaths") {
      what <- c("death_0", "death_55", "death_65", "death_75", "death_chr")
      what_names <- c("Deaths under 25s", "Deaths 25 to 54", "Deaths 55 to 64",
                      "Deaths 65 to 74", "Deaths 75+")
    }

    for (r in regions) {
      spim_plot_log_traj_by_age_region(r, dat, yield, what, model_type)
    }

    mtext(side = 3, text = what_names,
          line = 0.5, outer = TRUE, cex = 0.9, at = c(0.1, 0.32, 0.5, 0.7, 0.9))
    mtext(side = 2, text = spim_region_name(regions),
          line = -0.5, outer = TRUE, cex = 0.7,
          at = c(0.1, 0.24, 0.38, 0.53, 0.67, 0.81, 0.96))
  }

}


##' Plot variant
##'
##' @title Plot variant
##'
##' @param dat Combined data set
##'
##' @param regions Vector of regions to plot
##'
##' @param voc_name The name of the emerging VOC
##'
##' @param date_min Starting date for plot
##'
##' @param date_max End date for plot
##'
##' @export
spim_plot_variant <- function(dat, regions, voc_name,
                              date_min = NULL, date_max = NULL) {

  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in regions) {
    spim_plot_variant_region(r, dat, voc_name, date_min, date_max)
  }
}


spim_plot_variant_region <- function(region, dat, voc_name, date_min,
                                     date_max) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  pos_col <- cols$blue
  dcols <- c(cols$orange, cols$brown)

  n_non_variant <- data$fitted[, "strain_non_variant"]
  ntot <- data$fitted[, "strain_tot"]
  if (all(is.na(ntot))) {
    ## if na, switch to non-fitted data
    n_non_variant <- data$full[, "strain_non_variant"]
    ntot <- data$full[, "strain_tot"]
    dcols[1] <- cols$green2
  }

  npos <- ntot - n_non_variant
  npos[is.na(npos)] <- 0
  ntot[is.na(ntot)] <- 0
  dx <- as.Date(data$fitted$date_string)

  trajectories <- sample$trajectories$state
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  if ("time_index" %in% names(sample$info)) {
    parent <- sample$info$time_index$parent
    trajectories <- trajectories[, , -parent, drop = FALSE]
    x <- x[-parent]
  }

  tot <- trajectories["sympt_cases_inc", , ]
  pos <- tot - trajectories["sympt_cases_non_variant_inc", , ]

  res <- pos / tot * 100

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  #data
  cis <- Hmisc::binconf(x = npos, n = ntot) * 100
  dy <- cis[, "PointEst"]
  lower <- cis[, "Lower"]
  upper <- cis[, "Upper"]
  dy[ntot == 0] <- NA

  pos_cols <- c(mix_cols(pos_col, "white", 0.7),
                mix_cols(pos_col, "white", 0.495))

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  date_min <- date_min %||% min(x[-1L])
  date_max <- date_max %||% dat$info$date

  xlim <- c(date_min, date_max)
  ylim <- c(0, 100)

  plot(date_min, 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = paste0(voc_name, " proportion (%)"),
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
  segments(x0 = dx, y0 = lower, y1 = upper, col = "grey60")
  points(dx, dy, pch = 23, bg = dcols[1], col = dcols[2], cex = 0.8, lwd = 0.6)

  segments(x0 = max(dx), y0 = 0, y1 = 100, lwd = 3, col = "white")

}


spim_plot_log_traj_by_age_region <- function(region, dat, yield, what,
                                             model_type) {

  for (w in what) {
    spim_plot_log_traj_by_age_region1(region, dat, yield, w, model_type)
  }

}


spim_plot_log_traj_by_age_region1 <- function(region, dat, yield, what,
                                              model_type) {

  if (yield == "pillar2") {

    sample <- dat$samples[[region]]$trajectories
    trajectories <- sample$state
    x <- sircovid::sircovid_date_as_date(sample$date)

    data <- dat$data[[region]]$fitted
    date <- dat$info$date
    cols <- spim_colours()
    pos_col <- cols$blue
    dcols <- c(cols$orange, cols$brown)
    dx <- as.Date(data$date_string)
    date_min <- as.Date("2020-05-15")
    ps <- seq(0.025, 0.975, 0.005)

    if (model_type == "BB") {
      res <- trajectories[paste0("pillar2_positivity_", what), , ]
      ylim <- c(0, 40)

      npos <- data[, paste0("pillar2_", what, "_pos")]
      ntot <- data[, paste0("pillar2_", what, "_tot")]
      npos[is.na(npos)] <- 0
      ntot[is.na(ntot)] <- 0

      cis <- Hmisc::binconf(x = npos, n = ntot) * 100
      dy <- cis[, "PointEst"]
      dy[ntot == 0] <- NA
      qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
    }

    if (model_type == "NB") {
      dy <- data[, paste0("pillar2_", what, "_cases")]
      res <- trajectories[paste0("pillar2_cases_", what), , ]
      qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

      xlim <- c(date_min, dat$info$date)
      ymax <- max(dy[dx >= xlim[1] & dx <= xlim[2]],
                  qs[, x >= xlim[1] & x <= xlim[2]], na.rm = TRUE)
      ylim <- c(0, ymax)
    }

    pos_cols <- c(mix_cols(pos_col, "white", 0.7),
                  mix_cols(pos_col, "white", 0.495))

    if (what == "under15") {
      ylab <- spim_region_name(region)
    } else (
      ylab <- ""
    )

    oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
    on.exit(oo)

    plot(date_min, 0, type = "n",
         xlim = c(date_min, dat$info$date),
         ylim = ylim,
         las = 1,
         xlab = "", ylab = ylab)

    ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
             horiz = FALSE, leg = FALSE)
    lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
    if (model_type == "BB") {
      lines(dx, dy, col = grDevices::grey(0.2), lend = 1)
    }
    if (model_type == "NB") {
      points(dx, dy, pch = 23, bg = dcols[1], col = dcols[2], cex = 0.8,
             lwd = 0.6)
    }
    segments(x0 = max(dx), y0 = 0, y1 = 100, lwd = 3, col = "white")

  } else {
    if (grepl("^death_", what)) {
      samples <- dat$deaths[[region]]
    }
    if (grepl("^adm_", what)) {
      samples <- dat$admissions[[region]]
    }

    date <- dat$info$date

    cols <- spim_colours()
    now_col <- cols$blue
    dcols <- c(cols$orange, cols$green2)

    # Get model results
    res <- data.frame(count = samples$output_t[-1L, what],
                      lb = samples$lower_bound[-1L, what],
                      ub = samples$upper_bound[-1L, what])

    date_res <- dat$samples[[region]]$trajectories$date[-c(1, 2)]
    date_res <- sircovid::sircovid_date_as_date(date_res)
    res$date_res <- date_res
    res <- res %>% dplyr::filter(date_res <= date)

    # Get data
    data <- samples$data %>%
      dplyr::select(date_res = "date", count = get("what")) %>%
      dplyr::mutate(date_res = as.Date(date_res))
    suppressMessages(
      data <- dplyr::left_join(as.data.frame(date_res), as.data.frame(data))
    )

    # Vectors for plotting
    x_nowcast <- res$date_res
    y_nowcast <- res$count + 1
    lb <- res$lb + 1
    ub <- res$ub + 1

    dx <- data$date
    dy <- data$count + 1

    xlim <- c(min(date_res), date)
    ylim <- c(1, max(res$ub, na.rm = TRUE))

    oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
    on.exit(oo)

    plot(xlim[1], 1, type = "n",
         xlim = xlim,
         ylim = ylim,
         log = "y",
         xlab = "", ylab = "log scale", cex.lab = 0.8)
    now_cols <- c(mix_cols(now_col, "white", 0.7),
                  mix_cols(now_col, "white", 0.495))

    polygon(c(x_nowcast, rev(x_nowcast)), c(lb, rev(ub)), col = now_cols)
    lines(x_nowcast, y_nowcast, col = now_col, lty = 1, lwd = 1.5, lend = 1)
    points(dx, dy, pch = 23, bg = dcols[1], col = dcols[2], cex = 0.7,
           lwd = 0.6)
  }
}


spim_plot_prop_susceptible_region <- function(region, dat, ymin,
                                              plot_legend = FALSE) {

  sample <- dat$samples[[region]]
  cols <- spim_colours()
  S_col <- cols$purple
  eff_S_col <- cols$blue
  alpha <- 0.3

  N0 <- sum(sircovid::lancelot_parameters(1, region)$population)

  trajectories <- sample$trajectories$state

  index_S <- grep("^S_", rownames(sample$trajectories$state))
  S <- apply(trajectories[index_S, , ], c(2, 3), sum)
  prop_S <- S / N0 * 100

  prop_eff_S <- trajectories["eff_S", , ] / N0 * 100

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(prop_S,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
  qs_eff <- apply(prop_eff_S,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  S_cols <- add_alpha(rep(S_col, 2), alpha)
  eff_S_cols <- add_alpha(rep(eff_S_col, 2), alpha)


  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  xlim <- c(min(x[-1L]), max(x[-1L]))
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = c(ymin, 100),
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "Date", ylab = "Proportion susceptible (%)")

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = S_cols,
           horiz = FALSE, leg = FALSE)
  ci_bands(qs_eff[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = eff_S_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = S_col, lty = 1, lwd = 1.5, lend = 1)
  lines(x, qs_eff["50.0%", ], col = eff_S_col, lty = 1, lwd = 1.5, lend = 1)

  if (plot_legend) {
    leg_cols <- c(S_col, eff_S_col)
    if (plot_legend) {
      leg_cols <- c(cols$cyan, cols$purple)
      legend("topright", legend = c("All", "Effective"),
             cex = 1, x.intersp = 2, ncol = 1,
             fill = add_alpha(leg_cols, alpha * 2),
             border = leg_cols,
             bty = "n")
    }
  }

}


spim_plot_incidence_region <- function(region, dat, per_1000 = FALSE) {

  sample <- dat$samples[[region]]
  title <- spim_region_name(region)
  cols <- spim_colours()
  traj_cols <- c(cols$sky_blue1, cols$sky_blue2)
  line_col <- cols$cyan2

  res <- sample$trajectories$state["infections_inc", , -1L]
  if (per_1000) {
    res <- res / sum(sircovid::lancelot_parameters(1, region)$N_tot) * 1000
  }
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date[-1L])

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)


  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  xlim <- c(min(x), max(x))
  if (per_1000) {
    ylim <- c(0, 12)
    ylab <- "Incidence per day per 1000"
  } else {
    ylim <- c(0, max(res * 1.1,  na.rm = TRUE))
    ylab <- "Incidence per day"
  }

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(title),
       font.main = 1,
       xlab = "Date", ylab = ylab,
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x,
           cols = traj_cols, horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = line_col, lty = 1, lwd = 1.5)

}


spim_plot_pillar2_positivity_region <- function(region, dat, age_band,
                                                date_min, ymax,
                                                add_betas = FALSE,
                                                hard_xlim = FALSE,
                                                data_by = NULL) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  pos_col <- cols$blue
  if (is.null(data_by)) {
    data_by <- "standard"
  }

  if (age_band == "all") {
    npos <- data$fitted[, "pillar2_pos"]
    ntot <- data$fitted[, "pillar2_tot"]
    dcol <- grDevices::grey(0.2)
    if (all(is.na(ntot))) {
      npos <- data$full[, "pillar2_pos"]
      ntot <- data$full[, "pillar2_tot"]
      dcol <- grDevices::grey(0.6)
    }
  } else {
    npos <- data$fitted[, paste0("pillar2_", age_band, "_pos")]
    ntot <- data$fitted[, paste0("pillar2_", age_band, "_tot")]
    dcol <- grDevices::grey(0.2)
    if (all(is.na(ntot))) {
      npos <- data$full[, paste0("pillar2_", age_band, "_pos")]
      ntot <- data$full[, paste0("pillar2_", age_band, "_tot")]
      dcol <- grDevices::grey(0.6)
    }
  }

  npos[is.na(npos)] <- 0
  ntot[is.na(ntot)] <- 0
  dx <- as.Date(data$fitted$date_string)

  trajectories <- sample$trajectories$state
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  if (hard_xlim) {
    predicted <- sample$trajectories$predicted
    trajectories <- trajectories[, , !predicted]
    x <- x[!predicted]
  }

  if (age_band == "all") {
    res <- trajectories["pillar2_positivity", , ]
  } else {
    res <- trajectories[paste0("pillar2_positivity_", age_band), , ]
  }

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  if (data_by == "rolling week") {
    npos <- stats::filter(npos, rep(1 / 7, 7))
    ntot <- stats::filter(ntot, rep(1 / 7, 7))
    dy <- npos / ntot * 100

  } else {
    if (data_by == "week") {
      i_week <- lubridate::week(dx)
      npos <- tapply(npos, i_week, sum, na.rm = TRUE)
      ntot <- tapply(ntot, i_week, sum, na.rm = TRUE)
      dx <- as.Date(sircovid::sircovid_date_as_date(unique(i_week) * 7 - 3.5))

    }

    cis <- Hmisc::binconf(x = npos, n = ntot) * 100
    dy <- cis[, "PointEst"]
    lower <- cis[, "Lower"]
    upper <- cis[, "Upper"]
  }

  #data
  dy[ntot == 0] <- NA

  pos_cols <- c(mix_cols(pos_col, "white", 0.7),
                mix_cols(pos_col, "white", 0.495))

  if (age_band == "all") {
    age_lab <- "all ages"
  } else if (age_band == "over25") {
    age_lab <- "over 25"
  } else if (age_band == "under15") {
    age_lab <- "under15"
  } else if (age_band == "80_plus") {
    age_lab <- "80+"
  } else {
    age_lab <- gsub("_", " to ", age_band)
  }

  if (region %in% c("scotland", "northern_ireland")) {
    ylab <- paste0("Pillar 1 & 2 ", age_lab, " proportion positive (%)")
  } else {
    ylab <- paste0("Pillar 2 ", age_lab, " proportion positive (%)")
  }

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  xlim <- c(date_min, dat$info$date)
  ylim <- c(0, ymax)

  plot(date_min, 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = ylab,
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
  if (data_by == "rolling week") {
    lines(dx, dy, col = dcol)
  } else {
    lines(dx, dy, col = dcol, lend = 1)

  }

  segments(x0 = max(dx), y0 = 0, y1 = 100, lwd = 3, col = "white")

  if (add_betas == TRUE) {
    abline(v = sircovid::sircovid_date_as_date(sample$info$beta_date),
           lty = 2, col = cols$puce)
  }

}


spim_plot_pillar2_cases_region <- function(region, dat, age_band, date_min,
                                           add_betas = FALSE) {
  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  pos_col <- cols$blue
  dcols <- c(cols$orange, cols$brown)


  if (age_band == "all") {
    dy <- data$fitted[, "pillar2_cases"]
    if (all(is.na(dy))) {
      dy <- data$full[, "pillar2_cases"]
      dcols[1] <- cols$green2
    }
  } else {
    dy <- data$fitted[, paste0("pillar2_", age_band, "_cases")]
    if (all(is.na(dy))) {
      dy <- data$full[, paste0("pillar2_", age_band, "_cases")]
      dcols[1] <- cols$green2
    }
  }

  trajectories <- sample$trajectories$state
  if (age_band == "all") {
    res <- trajectories["pillar2_cases", , ]
  } else {
    res <- trajectories[paste0("pillar2_cases_", age_band), , ]
  }

  if (age_band == "all") {
    age_lab <- "all ages"
  } else if (age_band == "over25") {
    age_lab <- "over 25"
  } else if (age_band == "under15") {
    age_lab <- "under15"
  } else if (age_band == "80_plus") {
    age_lab <- "80+"
  } else {
    age_lab <- gsub("_", " to ", age_band)
  }

  if (region %in% c("scotland", "northern_ireland")) {
    ylab <- paste0("Pillar 1 & 2 ", age_lab, " cases")
  } else {
    ylab <- paste0("Pillar 2 ", age_lab, " cases")
  }

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  #data
  dx <- as.Date(data$fitted$date_string)

  xlim <- c(date_min, dat$info$date)
  ymax <- max(dy[dx >= xlim[1] & dx <= xlim[2]],
              qs[, x >= xlim[1] & x <= xlim[2]], na.rm = TRUE)
  ylim <- c(0, ymax)

  pos_cols <- c(mix_cols(pos_col, "white", 0.7),
                mix_cols(pos_col, "white", 0.495))


  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "Date", ylab = ylab)

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
  points(dx, dy, pch = 23, bg = dcols[1], col = dcols[2], cex = 0.8, lwd = 0.6)

  segments(x0 = max(dx), y0 = 0, y1 = 100, lwd = 3, col = "white")

  if (add_betas == TRUE) {
    abline(v = sircovid::sircovid_date_as_date(sample$info$beta_date),
           lty = 2, col = cols$puce)
  }

}


spim_plot_react_region <- function(region, dat, date_min, ymax,
                                   age_band, add_betas = FALSE) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  days_to_agg <- 3
  min_agg_tot <- 200
  cols <- spim_colours()
  pos_col <- cols$blue
  dcols <- c(cols$orange, cols$brown)

  if (is.null(age_band)) {
    npos <- data$fitted[, "react_pos"]
    ntot <- data$fitted[, "react_tot"]
  } else {
    npos <- data$fitted[, paste0("react_", age_band, "_pos")]
    ntot <- data$fitted[, paste0("react_", age_band, "_tot")]
  }

  if (all(is.na(npos))) {
    if (is.null(age_band)) {
      npos <- data$full[, "react_pos"]
      ntot <- data$full[, "react_tot"]
    } else {
      npos <- data$full[, paste0("react_", age_band, "_pos")]
      ntot <- data$full[, paste0("react_", age_band, "_tot")]
    }
    dcols[1] <- cols$green2
  }
  npos[is.na(npos)] <- 0
  ntot[is.na(ntot)] <- 0

  dx <- as.Date(data$fitted$date_string)

  ##aggregate
  agg_dates <- data.frame(start = dx[seq(1, length(dx), days_to_agg)])
  agg_dates$end <- agg_dates$start + days_to_agg - 1
  agg_dates$end[length(agg_dates$end)] <- dx[length(dx)]
  agg_dates$mid <- floor(rowMeans(apply(agg_dates, 2, sircovid::sircovid_date)))
  agg_dates$mid <- sircovid::sircovid_date_as_date(agg_dates$mid)

  agg_dates$end <- as.Date(agg_dates$end)

  aggregate_react <- function(i) {
    agg_pos <- sum(npos[dx >= agg_dates$start[i] & dx <= agg_dates$end[i]])
    agg_tot <- sum(ntot[dx >= agg_dates$start[i] & dx <= agg_dates$end[i]])
    if (agg_tot > min_agg_tot) {
      agg_out <- c(agg_pos, agg_tot)
    } else {
      agg_out <- c(0, 0)
    }
    agg_out
  }

  agg_dates$npos <- rep(0, length(agg_dates$mid))
  agg_dates$ntot <- rep(0, length(agg_dates$mid))
  agg_dates[, c("npos", "ntot")] <- t(sapply(seq_len(length(agg_dates$mid)),
                                             aggregate_react))


  cis <- Hmisc::binconf(x = agg_dates$npos, n = agg_dates$ntot) * 100
  dy <- cis[, "PointEst"]
  lower <- cis[, "Lower"]
  upper <- cis[, "Upper"]
  dy[agg_dates$ntot == 0] <- NA
  dx <- agg_dates$mid

  trajectories <- sample$trajectories$state

  model_params <- sample$predict$transform(sample$pars[1, ])
  model_params <- model_params[[length(model_params)]]$pars

  if (is.null(age_band)) {
    pos <- trajectories["react_pos", , ]
    neg <- (sum(model_params$N_tot_react) - pos)
    ylab <- "REACT proportion positive (%)"
  } else {
    pos <- trajectories[paste0("react_", age_band, "_pos"), , ]
    neg <- (sum(model_params[[paste0("N_", age_band, "_react")]]) - pos)
    age <- gsub("_", " to ", age_band)
    ylab <- paste0("REACT proportion positive ", age, " (%)")
  }

  res <- (pos * model_params$react_sensitivity +
            neg * (1 - model_params$react_specificity)) / (pos + neg) * 100


  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  pos_cols <- c(mix_cols(pos_col, "white", 0.7),
                mix_cols(pos_col, "white", 0.495))


  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  xlim <- c(date_min, dat$info$date)
  ylim <- c(0, ymax)

  plot(date_min, 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "Date", ylab = ylab,
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
  segments(x0 = dx, y0 = lower, y1 = upper, col = "grey60")
  points(dx, dy, pch = 23, bg = dcols[1], col = dcols[2], cex = 0.8, lwd = 0.6)

  segments(x0 = as.Date(date), y0 = 0, y1 = 100, lwd = 3, col = "white")

  if (add_betas == TRUE) {
    abline(v = sircovid::sircovid_date_as_date(sample$info$beta_date),
           lty = 2, col = cols$puce)
  }

}

spim_plot_serology_region <- function(region, dat, sero_flow, ymax,
                                      plot_legend = FALSE) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  alpha <- 0.35

  extract_serodates <- function(data, cap = 150) {
    sero <- data$full[[paste0("sero_tot_15_64_", sero_flow)]] > cap
    sero[is.na(sero)] <- FALSE
    n <- length(sero)

    lag <- seq(4, n)

    start <- sero[lag] & !sero[lag - 1] & !sero[lag - 2] & !sero[lag - 3]
    end <- sero[lag - 3] & !sero[lag - 2] & !sero[lag - 1] & !sero[lag]


    dates <- data$full$date_string
    res <- data.frame(start = dates[lag][start])

    end_dates <- dates[lag - 3][end]
    if (length(end_dates) < length(res$start)) {
      end_dates <- c(end_dates, date)
    }
    res$end <- end_dates
    res
  }

  sero_dates <- extract_serodates(data)
  if (nrow(sero_dates) > 0) {
    tol <- 2

    sero_dates$start <- sero_dates$start - tol
    sero_dates$end <- sero_dates$end + tol

    rownames(data$fitted) <- data$fitted$date_string

    sero_data <-
      data$fitted[, paste0(c("sero_tot_15_64_", "sero_pos_15_64_"), sero_flow)]
    colnames(sero_data) <- c("ntot", "npos")
    sero_data[is.na(sero_data)] <- 0

    summ_serodata <- sapply(X = seq_rows(sero_dates), function(i) {
      w <- seq(from = as.Date(sero_dates[i, "start"]),
               to = as.Date(sero_dates[i, "end"]), 1)
      w <- as.character(as.Date(w))
      colSums(sero_data[w, ], na.rm = TRUE)
    })

    summ_serodata <- data.frame(sero_dates,
                                t(summ_serodata),
                                stringsAsFactors = FALSE)

    summ_serodata <- cbind(summ_serodata,
                           with(summ_serodata,
                                Hmisc::binconf(x = npos, n = ntot) * 100))
    summ_serodata$mid <- (as.numeric(as.Date(summ_serodata$start)) +
                            as.numeric(as.Date(summ_serodata$end))) / 2
  }

  p <- sample$predict$transform(sample$pars[1, ])
  p <- p[[length(p)]]$pars
  sero_sensitivity <- p[[paste0("sero_sensitivity_", sero_flow)]]
  sero_specificity <- p[[paste0("sero_specificity_", sero_flow)]]
  sero_pos <- sample$trajectories$state[paste0("sero_pos_", sero_flow), , ]

  res <- (sero_sensitivity * sero_pos +
            (1 - sero_specificity) * (p$N_tot_15_64 - sero_pos)) /
    p$N_tot_15_64 * 100
  res_infs <- sample$trajectories$state["infections", , ] / p$N_tot_all * 100

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
  qs_infs <- apply(res_infs,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  pos_cols <- add_alpha(rep(cols$purple, 2), alpha)
  inf_cols <- add_alpha(rep(cols$cyan, 2), alpha)

  ylim <- c(0, ymax)
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[x <= date]
  xlim <- c(min(x[-1L]), max(x[-1L]))

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = "",
       font.main = 1,
       xlab = "", ylab = "Cumulative proportion (%)",
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  title(main = toupper(spim_region_name(region)), adj = 0, font.main = 1,
        line = 0.5, cex.main = 1)

  p_labels <- c("2.5%", "25.0%", "75.0%", "97.5%")

  ci_bands(qs[p_labels, seq_along(x)], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  ci_bands(qs_infs[p_labels, seq_along(x)], x, cols = inf_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs_infs["50.0%", seq_along(x)], col = cols$cyan,
        lty = 1, lwd = 1.5, lend = 1)
  lines(x, qs["50.0%", seq_along(x)], col = cols$purple,
        lty = 1, lwd = 1.5, lend = 1)
  if (nrow(sero_dates) > 0) {
    with(summ_serodata, {
      segments(x0 = mid, y0 = Lower, y1 = Upper, lty = 1, lend = 1)
      points(x = mid, y = PointEst, pch = 18)
      points(x = rep(mid, 2), y = c(Lower, Upper), pch = "-")
    })
  }

  if (plot_legend) {
    leg_cols <- c(cols$cyan, cols$purple)
    legend("topleft", legend = c("Infected", "Seropositive"),
           cex = 1, x.intersp = 2, ncol = 1,
           fill = add_alpha(leg_cols, alpha * 2),
           border = leg_cols,
           bty = "n")
  }

}


spim_plot_trajectories_region <- function(region, dat, what = NULL,
                                          date_min = NULL, age_band = NULL,
                                          with_forecast = TRUE,
                                          add_betas = FALSE) {
  if (is.null(what)) {
    what <- c("deaths_hosp", "deaths_comm", "deaths_carehomes",
              "icu", "general", "admitted", "diagnoses")
  }
  for (w in what) {
    spim_plot_trajectories_region1(
      w, region, dat, date_min, age_band,
      with_forecast = with_forecast, add_betas = add_betas)
  }
}




spim_plot_ifr_t_region <- function(region, dat, ifr_t_type, ymax,
                                   forecast_until, include_forecast = TRUE,
                                   add = FALSE) {

  ifr_t <- dat$ifr_t[[region]][[ifr_t_type]][-1L, ]
  col <- spim_colours()$blue

  x <- sircovid::sircovid_date_as_date(dat$ifr_t[[region]]$date[-1L, 1])

  if (!include_forecast) {
    ifr_t <- ifr_t[x <= forecast_until, ]
    x <- x[x <= forecast_until]
  }

  rownames(ifr_t) <- as.character(x)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(ifr_t,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)

  ## remove NAs for plotting purposes
  idx_not_na <- which(!is.na(colSums(qs)))
  x <- x[idx_not_na]
  qs <- qs[, idx_not_na]

  if (!add) {
    oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
    on.exit(oo)
    xlim <- c(x[1], as.Date(forecast_until))
    plot(xlim[1], 0, type = "n",
         xlim = xlim,
         ylim = c(0, ymax),
         main = toupper(spim_region_name(region)),
         font.main = 1,
         xlab = "", ylab = "IFR",
         xaxt = "n")
    segments(x0 = xlim[1], x1 = xlim[2], y0 = seq(0, 4, 0.5),
             col = grDevices::grey(0.9))
    axis.Date(1, at = seq(as.Date("2020-04-01"), as.Date(forecast_until),
                          by = "2 month"), format = "%b")
  }
  cols <- c(mix_cols(col, "white", 0.7),
            mix_cols(col, "white", 0.495))

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = col, lty = 1, lwd = 1.5, lend = 1)

}


spim_plot_alos_region <- function(region, dat, ymin, ymax, forecast_until,
                                  include_forecast = TRUE,
                                  add = FALSE) {

  ALOS <- dat$ifr_t[[region]]$ALOS[-1L, ]
  col <- spim_colours()$blue

  x <- sircovid::sircovid_date_as_date(dat$ifr_t[[region]]$date[-1L, 1])

  if (!include_forecast) {
    ALOS <- ALOS[x <= forecast_until, ]
    x <- x[x <= forecast_until]
  }

  rownames(ALOS) <- as.character(x)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(ALOS,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)

  ## remove NAs for plotting purposes
  idx_not_na <- which(!is.na(colSums(qs)))
  x <- x[idx_not_na]
  qs <- qs[, idx_not_na]

  if (is.null(ymax)) {
    ymax <- max(qs, na.rm = TRUE)
  }
  if (is.null(ymin)) {
    ymin <- 0
  }

  if (!add) {
    oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
    on.exit(oo)
    xlim <- c(x[1], as.Date(forecast_until))
    plot(xlim[1], 0, type = "n",
         xlim = xlim,
         ylim = c(ymin, ymax),
         main = toupper(spim_region_name(region)),
         font.main = 1,
         xlab = "", ylab = "Average length of stay",
         xaxt = "n")
    axis.Date(1, at = seq(as.Date("2020-04-01"), as.Date(forecast_until),
                          by = "2 month"), format = "%b")
  }
  cols <- c(mix_cols(col, "white", 0.7),
            mix_cols(col, "white", 0.495))

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = col, lty = 1, lwd = 1.5, lend = 1)

}

spim_plot_Rt_region <- function(region, dat, rt_type, forecast_until,
                                variant, add_betas, multistrain) {

  beta_date <- dat$info$beta_date
  if (variant == "weighted") {
    sample_Rt <- dat$rt[[region]][[rt_type]][-1L, ]
    x <- sircovid::sircovid_date_as_date(dat$rt[[region]]$date[-1L, 1])
  } else {
    if (variant == "delta") {
      sample_Rt <- dat$variant_rt[[region]][[rt_type]][-1L, 2, ]
      x <-
        sircovid::sircovid_date_as_date(dat$variant_rt[[region]]$date[-1L, 1])
    } else if (variant == "all") {
      if (!multistrain) {
        sample_non_variant <- dat$rt[[region]][[rt_type]][-1L, ]
        sample_variant <- dat$rt[[region]][[rt_type]][-1L, ]
        x <- sircovid::sircovid_date_as_date(dat$rt[[region]]$date[-1L, 1])
      } else {
        sample_non_variant <- dat$variant_rt[[region]][[rt_type]][-1L, 1, ]
        sample_variant <- dat$variant_rt[[region]][[rt_type]][-1L, 2, ]
        x <-
          sircovid::sircovid_date_as_date(dat$variant_rt[[region]]$date[-1L, 1])
      }
    } else {
      sample_Rt <- dat$variant_rt[[region]][[rt_type]][-1L, 1, ]
      x <-
        sircovid::sircovid_date_as_date(dat$variant_rt[[region]]$date[-1L, 1])
    }
  }

  if (rt_type == "beta") {
    ylab <- expression(beta[t])
  } else if (grepl("^eff_", rt_type)) {
    ## TODO: this is printing oddly on the graph
    ylab <- paste("eff", expression(R[t]), variant)
  } else {
    ylab <- paste(expression(R[t]), variant)
  }

  col <- spim_colours()$blue

  ps <- seq(0.025, 0.975, 0.005)

  if (variant != "all") {
    rownames(sample_Rt) <- as.character(x)
    qs <- apply(sample_Rt,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)
    ## remove NAs for plotting purposes
    idx_not_na <- which(!is.na(colSums(qs)))
    x <- x[idx_not_na]
    qs <- qs[, idx_not_na]
  } else {
    rownames(sample_non_variant) <- as.character(x)
    rownames(sample_variant) <- as.character(x)
    qs <- NULL
    qs[["non_variant"]] <- apply(sample_non_variant,  MARGIN = 1,
                                 FUN = quantile, ps, na.rm = TRUE)
    qs[["variant"]] <- apply(sample_variant,  MARGIN = 1,
                             FUN = quantile, ps, na.rm = TRUE)
    ## remove NAs for plotting purposes
    idx_not_na <- which(!is.na(colSums(qs[[1]])))
    x <- x[idx_not_na]
    qs[[1]] <- qs[[1]][, idx_not_na]
    qs[[2]] <- qs[[2]][, idx_not_na]
  }

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(par(oo))

  xlim <- c(x[1], as.Date(forecast_until))
  if (rt_type == "beta") {
    ylim <- c(0, 0.15)
  } else {
    ylim <- c(0, 5)
  }
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "Date", ylab = ylab,
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  cols <- c(mix_cols(col, "white", 0.7),
            mix_cols(col, "white", 0.495))

  if (variant != "all") {
    ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = cols,
             horiz = FALSE, leg = FALSE)
    lines(x, qs["50.0%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  } else {
    ci_bands(qs[[1]][c("2.5%", "25.0%", "75.0%", "97.5%"), ], x,
             cols = c(spim_colours()$sky_blue, spim_colours()$sky_blue),
             horiz = FALSE, leg = FALSE)
    lines(x, qs[[1]]["50.0%", ], col = col, lty = 1, lwd = 1.5, lend = 1)

    ci_bands(qs[[2]][c("2.5%", "25.0%", "75.0%", "97.5%"), ], x,
             cols = c(spim_colours()$cyan, spim_colours()$cyan),
             horiz = FALSE, leg = FALSE)
    lines(x, qs[[2]]["50.0%", ], col = spim_colours()$blue, lty = 1,
          lwd = 1.5, lend = 1)
  }
  if (rt_type != "beta") {
    abline(h = 1, lty = 2)
  }
  if (add_betas) {
    abline(v = sircovid::sircovid_date_as_date(beta_date),
           lty = 2, col = spim_colours()$orange)
  }
}


spim_plot_trajectories_region1 <- function(what, region, dat, date_min,
                                           age_band, with_forecast = TRUE,
                                           add_betas = FALSE,
                                           main = NULL) {
  trajectories <- dat$samples[[region]]$trajectories
  date <- dat$info$date
  cols <- spim_colours()
  beta_date <- dat$info$beta_date

  trajnames <- c(deaths = "deaths_inc",
                 deaths_comm = "deaths_comm_inc",
                 deaths_hosp = "deaths_hosp_inc",
                 deaths_carehomes = "deaths_carehomes_inc",
                 icu = "icu",
                 hosp = "hosp",
                 general = "general",
                 diagnoses = "diagnoses_inc",
                 admitted = "admitted_inc",
                 all_admission = "all_admission_inc")
  labs <- c(deaths = "Daily deaths",
            deaths_comm = "Daily community deaths",
            deaths_carehomes = "Daily care home deaths",
            deaths_hosp = "Daily hospital deaths",
            icu = "ICU beds",
            general = "general beds",
            hosp = "Hospital beds",
            diagnoses = "Daily inpatient diagnoses",
            admitted = "Daily admissions",
            all_admission = ifelse(age_band == "all",
                                   "Daily admissions (all)",
                                   "Daily admissions"))

  if (age_band == "all") {
    ## remove first step as the first time period is typically much more
    ## than one day at the moment
    if (what == "deaths") {
      res <-
        trajectories$state["deaths_hosp_inc", , -1L] +
        trajectories$state["deaths_comm_inc", , -1L] +
        trajectories$state["deaths_carehomes_inc", , -1L]
    } else if (what == "all_admission") {
      res <-
        trajectories$state["diagnoses_inc", , -1L] +
        trajectories$state["admitted_inc", , -1L]
    } else {
      res <- trajectories$state[trajnames[what], , -1L]
    }
  } else {
    if (what == "deaths_hosp") {
      res <- trajectories$state[paste0("deaths_hosp_", age_band, "_inc"), , -1L]
    } else if (what == "all_admission") {
      res <-
        trajectories$state[paste0("all_admission_", age_band, "_inc"), , -1L]
    } else {
      stop(message(paste0("Cannot plot ", what, " by age")))
    }
    labs[what] <- paste(labs[what], gsub("_", " to ", age_band))
  }

  x <- sircovid::sircovid_date_as_date(trajectories$date[-1L])
  colnames(res) <- as.character(x)

  x_nowcast <- x[x <= date]
  x_forecast <- x[x > date]
  res_nowcast <- res[, as.character(x_nowcast)]
  res_forecast <- res[, as.character(x_forecast)]

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res, MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
  qs_forecast <- apply(res_forecast,  MARGIN = 2, FUN = quantile, ps,
                       na.rm = TRUE)
  qs_nowcast <- apply(res_nowcast,  MARGIN = 2, FUN = quantile, ps,
                      na.rm = TRUE)

  data <- dat$data[[region]]
  dx <- as.Date(data$fitted$date_string)

  if (age_band == "all") {
    dy <- data$fitted[, what]
    dy_extra <- data$full[, what]
  } else {
    if (what == "all_admission") {
      dy <- data$fitted[, paste0("all_admission_", age_band)]
      dy_extra <- data$full[, paste0("all_admission_", age_band)]
    } else {
      dy <- data$fitted[, paste0(what, "_", age_band)]
      dy_extra <- data$full[, what]
    }
  }
  dy_extra[!is.na(dy)] <- NA_integer_

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  xlim <- c(date_min %||% as.Date("2020-03-16"),
            if (with_forecast) max(x_forecast) else date)

  rn_res <- as.Date(colnames(res))
  ylim <- c(0, max(dy[which(dx >= xlim[1] & dx <= xlim[2])],
                   qs[, which(rn_res >= xlim[1] & rn_res <= xlim[2])],
                   na.rm = TRUE))
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = main,
       font.main = 1,
       xlab = "", ylab = labs[what],
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")
  now_cols <- c(mix_cols(cols$now, "white", 0.7),
                mix_cols(cols$now, "white", 0.495))
  fore_cols <- c(mix_cols(cols$forecast, "white", 0.7),
                 mix_cols(cols$forecast, "white", 0.495))

  #Extract every first date of month from x_nowcast
  firsts <- x_nowcast[!duplicated(substring(x_nowcast, 1, 7))]
  #Extract every first date of year from x_nowcast
  year_firsts <- x_nowcast[!duplicated(substring(x_nowcast, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  ci_bands(qs_nowcast[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x_nowcast,
           cols = now_cols, horiz = FALSE, leg = FALSE)
  lines(x_nowcast, qs_nowcast["50.0%", ], col = cols$now, lty = 1, lwd = 1.5,
        lend = 1)

  if (length(x_forecast) > 0) {
    ci_bands(qs_forecast[c("2.5%", "25.0%", "75.0%", "97.5%"), ],
             x_forecast, cols = fore_cols, horiz = FALSE, leg = FALSE)
    lines(x_forecast, qs_forecast["50.0%", ], col = cols$forecast, lty = 1,
          lwd = 1.5, lend = 1)
  }

  points(dx, dy, pch = 23, bg = cols$orange, col = cols$brown, cex = 0.7,
         lwd = 0.6)
  points(dx, dy_extra, pch = 23, bg = cols$green2, col = cols$brown, cex = 0.7,
         lwd = 0.6)

  if (add_betas) {
    abline(v = sircovid::sircovid_date_as_date(beta_date),
           lty = 2, col = cols$puce)
  }
}



spim_plot_effective_susceptible_region <- function(region, dat,
                                                   strain_names, strain_dates,
                                                   plot_legend = FALSE) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  alpha <- 0.35

  p <- sample$predict$transform(sample$pars[1, ])
  p <- p[[length(p)]]$pars

  state <- sample$trajectories$state
  state <- state[, , -1L]

  strain_cols <- khroma::colour("bright")(length(strain_names))
  strain_dates <- c(strain_dates, as.Date(date) + 1)
  res <- array(NA, dim = c(length(strain_names), dim(state)[2:3]))

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[-1L]

  res[1, ,  x < strain_dates[2]] <-
    state["effective_susceptible_1", , x < strain_dates[2]]
  for (i in seq_along(strain_names)[-1L]) {
    phase_dates <- (x >= strain_dates[i] & x < strain_dates[i + 1])
    res[i - 1, , phase_dates] <-
      state["effective_susceptible_1", , phase_dates]
    res[i, , phase_dates] <-
      state["effective_susceptible_2", , phase_dates]
  }

  res <- res / p$N_tot_all * 100

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = c(1, 3), FUN = quantile, ps, na.rm = TRUE)

  ylim <- c(0, 100)
  xlim <- c(min(x), max(x))

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = "",
       font.main = 1,
       xlab = "", ylab = "Effective susceptible proportion (%)",
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  title(main = toupper(spim_region_name(region)), adj = 0, font.main = 1,
        line = 0.5, cex.main = 1)

  p_labels <- c("2.5%", "25.0%", "75.0%", "97.5%")

  strain_dates <- c(strain_dates, strain_dates[length(strain_dates)])

  for (i in seq_along(strain_names)) {
    strain_span <- (x >= strain_dates[i] & x < strain_dates[i + 2])
      str_col <- strain_cols[i]
    ci_cols <- c(mix_cols(str_col, "white", 0.7),
                 mix_cols(str_col, "white", 0.495))
    ci_bands(qs[p_labels, i, strain_span], x[strain_span], cols = ci_cols,
             horiz = FALSE, leg = FALSE)
    lines(x[strain_span], qs["50.0%", i, strain_span], col = str_col,
          lty = 1, lwd = 1.5, lend = 1)
  }

  if (plot_legend) {
    leg_cols <- strain_cols
    legend("bottomleft", legend = strain_names,
           cex = 1, x.intersp = 2, ncol = 1,
           fill = add_alpha(leg_cols, alpha * 2),
           border = leg_cols,
           bty = "n")
  }

}


spim_plot_infections_per_strain_region <- function(region, dat,
                                                   strain_names, strain_dates,
                                                   plot_legend = FALSE) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  alpha <- 0.35

  p <- sample$predict$transform(sample$pars[1, ])
  p <- p[[length(p)]]$pars

  state <- sample$trajectories$state
  state <- state[, , -1L]

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[-1L]

  inf_names <- c(rbind(strain_names, paste0(strain_names, " (reinfection)")))

  inf_inc <- state[paste0("infections_inc_strain_", seq_len(4)), , ]

  strain_dates <- c(strain_dates, as.Date(date) + 1)
  res <- array(0, dim = c(length(strain_names) * 2, dim(inf_inc)[2:3]))
  row.names(res) <- inf_names

  res[strain_names[1], ,  x < strain_dates[2]] <-
    inf_inc["infections_inc_strain_1", , x < strain_dates[2]]
  for (i in seq_along(strain_names)[-1L]) {
    phase_dates <- (x >= strain_dates[i] & x < strain_dates[i + 1])
    res[strain_names[i - 1], , phase_dates] <-
      state["infections_inc_strain_1", , phase_dates]
    res[strain_names[i], , phase_dates] <-
      state["infections_inc_strain_2", , phase_dates]
    res[paste0(strain_names[i], " (reinfection)"), , phase_dates] <-
      state["infections_inc_strain_3", , phase_dates]
    res[paste0(strain_names[i - 1], " (reinfection)"), , phase_dates] <-
      state["infections_inc_strain_4", , phase_dates]
  }

  res <- apply(res, c(1, 3), mean)

  remove_rows <- rowSums(res) == 0
  inf_names <- inf_names[!remove_rows]
  res <- res[!remove_rows, ]

  xlim <- c(min(x), max(x))
  ylim <- c(0, 100)

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = "Proportion of daily infections (%)",
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  cols <- khroma::colour("bright")(length(inf_names))

  proportion_fill(x, res, cols)

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  #Plot gray line on 1st of every month
  abline(v = firsts[-1], col = "ivory2", lwd = 0.5)
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  if (plot_legend) {
    legend("bottomleft", legend = inf_names,
           cex = 0.8, x.intersp = 2, ncol = 1,
           fill = cols,
           bty = "o",
           bg = "white",
           inset = 0.02)
  }

}


spim_plot_vaccine_status_region <- function(region, dat, vacc_names,
                                            plot_legend = FALSE) {
  sample <- dat$samples[[region]]
  date <- dat$info$date

  p <- sample$predict$transform(sample$pars[1, ])
  p <- p[[length(p)]]$pars

  state <- sample$trajectories$state
  state <- state[, , -1L]

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[-1L]

  vacc_rows <- grep("^vaccine_status", rownames(state))
  res <- state[vacc_rows, , ]
  res <- apply(res, c(1, 3), mean)

  xlim <- c(min(x), max(x))
  ylim <- c(0, 100)

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = "Vaccine status (%)",
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  cols <- khroma::colour("bright")(length(vacc_names))

  proportion_fill(x, res, cols)

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  #Plot gray line on 1st of every month
  abline(v = firsts[-1], col = "ivory2", lwd = 0.5)
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  if (plot_legend) {
    legend("bottomleft", legend = vacc_names,
           cex = 0.8, x.intersp = 2, ncol = 1,
           fill = cols,
           bty = "o",
           bg = "white",
           inset = 0.02)
  }

}


spim_plot_infection_status_region <- function(region, dat,
                                              plot_legend = FALSE) {
  sample <- dat$samples[[region]]
  date <- dat$info$date

  p <- sample$predict$transform(sample$pars[1, ])
  p <- p[[length(p)]]$pars

  state <- sample$trajectories$state
  state <- state[, , -1L]

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[-1L]

  status_rows <- c("susceptible", "recovered_1", "recovered_2", "deaths")
  res <- state[status_rows, , ]
  res <- apply(res, c(1, 3), mean)
  infected <- p$N_tot_all - colSums(res)
  res <- rbind(res, infected)
  res <- res[c("susceptible", "recovered_1", "recovered_2",
               "infected", "deaths"), ]

  xlim <- c(min(x), max(x))
  ylim <- c(0, 100)

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = "Infection status (%)",
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  cols <- khroma::colour("bright")(nrow(res))

  proportion_fill(x, res, cols)

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  #Plot gray line on 1st of every month
  abline(v = firsts[-1], col = "ivory2", lwd = 0.5)
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  if (plot_legend) {
    legend("bottomleft",
           legend = c("Susceptible",
                      "Recovered (strain 1)",
                      "Recovered (strain 2)",
                      "Infected",
                      "Dead"),
           cex = 0.8, x.intersp = 2, ncol = 1,
           fill = cols,
           bty = "o",
           bg = "white",
           inset = 0.02)
  }

}


spim_plot_cumulative_attack_rate_region <- function(region, dat,
                                                    plot_legend = FALSE) {

  sample <- dat$samples[[region]]
  cols <- spim_colours()
  cum_AR_col <- cols$purple
  non_S_col <- cols$blue
  alpha <- 0.3

  N0 <- sum(sircovid::lancelot_parameters(1, region)$population)

  trajectories <- sample$trajectories$state

  cum_AR <- trajectories["infections", , ] / N0 * 100

  prop_non_S <- (N0 - trajectories["susceptible", , ]) / N0 * 100

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- seq(0.025, 0.975, 0.005)
  qs_cum_AR <- apply(cum_AR,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
  qs_non_S <- apply(prop_non_S,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  cum_AR_cols <- add_alpha(rep(cum_AR_col, 2), alpha)
  non_S_cols <- add_alpha(rep(non_S_col, 2), alpha)

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)

  xlim <- c(min(x[-1L]), max(x[-1L]))
  ylim <- c(0, 100)
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = "Infection status (%)",
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")

  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  #Plot gray line on 1st of every month
  abline(v = firsts[-1], col = "ivory2", lwd = 0.5)
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))

  ci_bands(qs_cum_AR[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x,
           cols = cum_AR_cols, horiz = FALSE, leg = FALSE)
  ci_bands(qs_non_S[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x,
           cols = non_S_cols, horiz = FALSE, leg = FALSE)
  lines(x, qs_cum_AR["50.0%", ], col = cum_AR_col, lty = 1, lwd = 1.5, lend = 1)
  lines(x, qs_non_S["50.0%", ], col = non_S_col, lty = 1, lwd = 1.5, lend = 1)

  if (plot_legend) {
    leg_cols <- c(cum_AR_col, non_S_col)
    legend("topright", legend = c("Cumulative AR", "Non-S proportion"),
           cex = 1, x.intersp = 2, ncol = 1,
           fill = add_alpha(leg_cols, alpha * 2),
           border = leg_cols,
           bty = "n")
  }

}


ci_bands <- function(quantiles, y, palette = NULL, cols = NULL, leg = TRUE,
                     leg_y = 0, leg_x = 1, horiz = TRUE, ...) {
  yy <- c(y, rev(y))
  yy <- c(yy, yy[1])
  n_bands <- (nrow(quantiles) - 1) / 2 + 1
  if (!is.null(palette)) {
    cols <- do.call(what = palette,
                    args = list(n = n_bands))
  }

  for (band in seq_len(n_bands)) {
    x1 <- quantiles[band, ]
    x2 <- quantiles[nrow(quantiles) + 1 - band, ]

    x2 <- rev(x2)
    x2 <- c(x2, x1[1])
    if (horiz) {
      polygon(y = yy,
              x = c(x1, x2),
              col = cols[band],
              border = NA)
    } else {
      polygon(x = yy,
              y = c(x1, x2),
              col = cols[band],
              border = NA)
    }

  }
  if (leg) {
    leg_cols <- which(row.names(quantiles) %in% leg)
    leg <- c(row.names(quantiles)[1], seq(5, 50, 5), "%")
    leg[seq(2, 10, 2)] <- ""
    legend(y = leg_y,
           x = leg_x,
           pch = 15,
           col = cols[leg_cols],
           legend = leg,
           border = NA,
           bty = "n", ...)
  }
}


proportion_fill <- function(x, y, cols) {
  n_y <- nrow(y)
  y <- apply(y, 2, cumsum)
  y <- rbind(rep(0, dim(y)[2]), y)
  y <- t(t(y) / y[nrow(y), ]) * 100

  for (i in seq_len(n_y)) {
    y1 <- y[i, ]
    y2 <- y[i + 1, ]
    polygon(x = c(x, rev(x)),
            y = c(y2, rev(y1)),
            col = cols[i])
  }
}
