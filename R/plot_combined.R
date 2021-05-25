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
##' @param with_forecast Logical, indicating if we should add the forecast
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_trajectories <- function(dat, regions, what, date_min = NULL,
                                   with_forecast = TRUE, add_betas = FALSE) {

  n_regions <- length(regions)
  op <- par(mfcol = c(length(what), length(regions)),
            oma = c(1, 1, 4, 1),
            mar = c(3, 3, 0.5, 0.5))
  on.exit(par(op))

  for (r in regions) {
    spim_plot_trajectories_region(r, dat, what, date_min,
                                  with_forecast = with_forecast,
                                  add_betas = add_betas)
  }

  mtext(side = 3, text = toupper(spim_region_name(regions)),
        line = 0.5, outer = TRUE, cex = 0.8,
        at = seq(1 / n_regions / 2, by = 1 / n_regions, length.out = n_regions))
}


##' Plot Rt
##'
##' @title Plot Rt
##'
##' @param dat Combined data set
##'
##' @param rt_type One of the valid Rt types (e.g., eff_Rt_all)
##'
##' @param forecast_until Optional date to forecast till
##'
##' @export
spim_plot_Rt <- function(dat, rt_type, forecast_until = NULL) {
  oo <- par(mfrow = c(2, 6), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))
  if (is.null(forecast_until)) {
    forecast_until <- dat$info$date
  }
  for (r in names(dat$rt)) {
    spim_plot_Rt_region(r, dat, rt_type, forecast_until)
  }
}


##' Plot serology
##'
##' @title Plot serology
##'
##' @param dat Combined data set
##'
##' @param sero_flow Number identifying which serology flow to use (1 or 2)
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @export
spim_plot_serology <- function(dat, sero_flow, ymax) {
  region_names <- sircovid::regions("all")

  oo <- par(mfrow = c(2, 5), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in region_names) {
    if (r == region_names[length(region_names)]) {
      spim_plot_serology_region(r, dat, sero_flow, ymax, TRUE)
    } else {
      spim_plot_serology_region(r, dat, sero_flow, ymax)
    }
  }
}


##' Plot Pillar 2 positivity
##'
##' @title Plot Pillar 2 positivity
##'
##' @param dat Combined data set
##'
##' @param date_min Starting date for plot
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_pillar2_positivity <- function(dat, date_min, ymax,
                                         add_betas = FALSE) {
  region_names <- sircovid::regions("all")

  oo <- par(mfrow = c(2, 5), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in region_names) {
    spim_plot_pillar2_positivity_region(r, dat, date_min, ymax, add_betas)
  }
}


##' Plot Pillar 2 cases
##'
##' @title Plot Pillar 2 cases
##'
##' @param dat Combined data set
##'
##' @param date_min Starting date for plot
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_pillar2_cases <- function(dat, date_min, add_betas = FALSE) {
  region_names <- sircovid::regions("all")

  oo <- par(mfrow = c(2, 5), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in region_names) {
    spim_plot_pillar2_cases_region(r, dat, date_min, add_betas)
  }
}


##' Plot REACT
##'
##' @title Plot REACT
##'
##' @param dat Combined data set
##'
##' @param date_min Starting date for plot
##'
##' @param ymax Maximum percentage on y-axis
##'
##' @param add_betas Logical, indicating if we should add betas
##'
##' @export
spim_plot_react <- function(dat, date_min, ymax, add_betas = FALSE) {
  region_names <- sircovid::regions("all")

  oo <- par(mfrow = c(2, 5), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))

  for (r in region_names) {
    spim_plot_react_region(r, dat, date_min, ymax, add_betas)
  }
}


spim_plot_pillar2_positivity_region <- function(region, dat, date_min, ymax,
                                                add_betas = FALSE,
                                                hard_xlim = FALSE,
                                                data_by = NULL) {
  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  pos_col <- cols$blue
  dcols <- c(cols$orange, cols$brown)
  alpha <- 0.3
  if (is.null(data_by)) {
    data_by <- "standard"
  }

  over25 <- TRUE

  npos <- data$fitted[, "pillar2_over25_pos"]
  ntot <- data$fitted[, "pillar2_over25_tot"]
  if (all(is.na(ntot))) {
    npos <- data$full[, "pillar2_over25_pos"]
    ntot <- data$full[, "pillar2_over25_tot"]

    ## if still na, switch to all ages
    if (all(is.na(ntot))) {
      npos <- data$fitted[, "pillar2_pos"]
      ntot <- data$fitted[, "pillar2_tot"]
      if (all(is.na(ntot))) {
        npos <- data$full[, "pillar2_pos"]
        ntot <- data$full[, "pillar2_tot"]
        dcols[1] <- cols$green2
      }
      over25 <- FALSE
    } else {
      dcols[1] <- cols$green2
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


  model_params <- sample$predict$transform(sample$pars[1, ])

  if ("p_NC" %in% colnames(sample$pars)) {
    p_NC <- sample$pars[, "p_NC"]
  } else {
    p_NC <- model_params$p_NC
  }

  if (over25) {
    pos <- trajectories["sympt_cases_over25_inc", , ]
    neg <- (sum(model_params$N_tot[6:19]) - pos) * p_NC
  } else {
    pos <- trajectories["sympt_cases_inc", , ]
    neg <- (sum(model_params$N_tot) - pos) * p_NC
  }

  res <- (pos * model_params$pillar2_sensitivity +
            neg * (1 - model_params$pillar2_specificity)) / (pos + neg) * 100



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

  pos_cols <- add_alpha(rep(pos_col, 2), alpha)

  if (over25) {
    if (region == "scotland") {
      ylab <- "Pillar 1 & 2 over 25 proportion positive (%)"
    } else {
      ylab <- "Pillar 2 over 25 proportion positive (%)"
    }
  } else {
    if (region == "scotland") {
      ylab <- "Pillar 1 & 2 proportion positive (%)"
    } else {
      ylab <- "Pillar 2 proportion positive (%)"
    }
  }

  par(mgp = c(1.7, 0.5, 0), bty = "n")
  plot(date_min, 0, type = "n",
       xlim = c(date_min, dat$info$date),
       ylim = c(0, ymax),
       las = 1,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = ylab)

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
  if (data_by == "rolling week") {
    lines(dx, dy, col = grey(0.2))
  } else {
    lines(dx, dy, col = grey(0.2), lend = 1)

  }

  segments(x0 = max(dx), y0 = 0, y1 = 100, lwd = 3, col = "white")

  if (add_betas == TRUE) {
    abline(v = as.Date(sample$info$beta_date), lty = 2, col = cols$puce)
  }

}


spim_plot_pillar2_cases_region <- function(region, dat, date_min,
                                           add_betas = FALSE) {
  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  cols <- spim_colours()
  pos_col <- cols$blue
  dcols <- c(cols$orange, cols$brown)
  alpha <- 0.3


  over25 <- TRUE

  dy <- data$fitted[, "pillar2_over25_cases"]
  if (all(is.na(dy))) {
    dy <- data$full[, "pillar2_over25_cases"]

    ## if still na, switch to all ages
    if (all(is.na(dy))) {
      dy <- data$fitted[, "pillar2_cases"]
      if (all(is.na(dy))) {
        dy <- data$full[, "pillar2_cases"]
        dcols[1] <- cols$green2
      }
      over25 <- FALSE
    } else {
      dcols[1] <- cols$green2
    }
  }

  trajectories <- sample$trajectories$state
  model_params <- sample$predict$transform(sample$pars[1, ])

  if ("phi_pillar2_cases" %in% colnames(sample$pars)) {
    phi_pillar2_cases <- sample$pars[, "phi_pillar2_cases"]
  } else {
    phi_pillar2_cases <- model_params$phi_pillar2_cases
  }

  if (over25) {
    pos <- trajectories["sympt_cases_over25_inc", , ] * phi_pillar2_cases
    if (region == "scotland") {
      ylab <- "Pillar 1 & 2 over 25 cases"
    } else {
      ylab <- "Pillar 2 over 25 cases"
    }
  } else {
    pos <- trajectories["sympt_cases_inc", , ] * phi_pillar2_cases
    if (region == "scotland") {
      ylab <- "Pillar 1 & 2 cases"
    } else {
      ylab <- "Pillar 2 cases"
    }
  }

  res <- pos

  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  #data
  dx <- as.Date(data$fitted$date_string)

  xlim <- c(date_min, dat$info$date)
  ymax <- max(dy[dx >= xlim[1] & dx <= xlim[2]],
              qs[, x >= xlim[1] & x <= xlim[2]], na.rm = TRUE)
  ylim <- c(0, ymax)

  pos_cols <- add_alpha(rep(pos_col, 2), alpha)


  par(mgp = c(1.7, 0.5, 0), bty = "n")
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
    abline(v = as.Date(sample$info$beta_date), lty = 2, col = cols$puce)
  }

}


spim_plot_react_region <- function(region, dat, date_min, ymax,
                                   add_betas = FALSE) {

  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date
  days_to_agg <- 3
  min_agg_tot <- 200
  cols <- spim_colours()
  pos_col <- cols$blue
  dcols <- c(cols$orange, cols$brown)
  alpha <- 0.3

  npos <- data$fitted[, "react_pos"]
  ntot <- data$fitted[, "react_tot"]
  if (all(is.na(npos))) {
    npos <- data$full[, "react_pos"]
    ntot <- data$full[, "react_tot"]
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

  pos <- trajectories["react_pos", , ]
  neg <- (sum(model_params$N_tot_react) - pos)

  res <- (pos * model_params$react_sensitivity +
            neg * (1 - model_params$react_specificity)) / (pos + neg) * 100


  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res,  MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)

  pos_cols <- add_alpha(rep(pos_col, 2), alpha)


  par(mgp = c(1.7, 0.5, 0), bty = "n")
  plot(date_min, 0, type = "n",
       xlim = c(date_min, dat$info$date),
       ylim = c(0, ymax),
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "Date", ylab = "REACT proportion positive (%)")

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = pos_cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = pos_col, lty = 1, lwd = 1.5, lend = 1)
  segments(x0 = dx, y0 = lower, y1 = upper, col = "grey60")
  points(dx, dy, pch = 23, bg = dcols[1], col = dcols[2], cex = 0.8, lwd = 0.6)

  segments(x0 = as.Date(date), y0 = 0, y1 = 100, lwd = 3, col = "white")

  if (add_betas == TRUE) {
    abline(v = as.Date(sample$info$beta_date), lty = 2, col = cols$puce)
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

    summ_serodata <- sapply(X = seq_len(nrow(sero_dates)), function(i) {
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
  par(mgp = c(1.7, 0.5, 0), bty = "n")
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[x <= date]
  xlim <- c(min(x[-1L]), max(x[-1L]))
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim, las = 1,
       main = "",
       font.main = 1,
       xlab = "", ylab = "Cumulative proportion (%)")
  title(main = toupper(spim_region_name(region)), adj = 0, font.main = 1,
        line = 0.5, cex.main = 1)

  if (nrow(sero_dates) > 0) {
    lapply(X = seq_len(nrow(summ_serodata)), FUN = function(i) {

      xx <- rep(unlist(summ_serodata[i, c("start", "end")]), each = 2)
      yy <- c(ylim, rev(ylim))

      polygon(x = xx, y = yy, col = grey(0.9), border = NA)
    })
  }

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
    legend("top", legend = c("Infected", "Seropositive"),
           cex = 1, x.intersp = 2, ncol = 2,
           fill = add_alpha(leg_cols, alpha * 2),
           border = leg_cols,
           box.col = "white")
  }

}


spim_plot_trajectories_region <- function(region, dat, what = NULL,
                                          date_min = NULL, with_forecast = TRUE,
                                          add_betas = FALSE) {
  if (is.null(what)) {
    what <- c("deaths_hosp", "deaths_comm", "deaths_carehomes",
              "icu", "general", "admitted", "diagnoses")
  }
  for (w in what) {
    spim_plot_trajectories_region1(
      w, region, dat, date_min,
      with_forecast = with_forecast, add_betas = add_betas)
  }
}





spim_plot_Rt_region <- function(region, dat, rt_type, forecast_until) {
  sample_Rt <- dat$rt[[region]][[rt_type]][-1L, ]

  if (grepl("^eff_", rt_type)) {
    ## TODO: this is printing oddly on the graph
    ylab <- paste("eff", expression(R[t]))
  } else {
    ylab <- expression(R[t])
  }

  ## sample_Rt,
  ## Rt_date,
  col <- spim_colours()$blue
  alpha <- 0.3

  x <- sircovid::sircovid_date_as_date(dat$rt[[region]]$date[-1L, 1])
  rownames(sample_Rt) <- as.character(x)

  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(sample_Rt,  MARGIN = 1, FUN = quantile, ps, na.rm = TRUE)

  ## remove NAs for plotting purposes
  idx_not_na <- which(!is.na(colSums(qs)))
  x <- x[idx_not_na]
  qs <- qs[, idx_not_na]

  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(par(oo))

  xlim <- c(x[1], as.Date(forecast_until))
  ylim <- c(0, 4)
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(spim_region_name(region)),
       font.main = 1,
       xlab = "Date", ylab = ylab)

  cols <- add_alpha(rep(col, 2), alpha)

  ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x, cols = cols,
           horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = col, lty = 1, lwd = 1.5, lend = 1)
  abline(h = 1, lty = 2)
}


spim_plot_trajectories_region1 <- function(what, region, dat, date_min,
                                           with_forecast = TRUE,
                                           add_betas = FALSE) {
  trajectories <- dat$samples[[region]]$trajectories
  date <- dat$info$date
  cols <- spim_colours()
  beta_date <- dat$info[[region]]$beta_date
  alpha <- 0.3

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
            all_admission = "Daily admissions (all)")

  ## currently remove first step as the first time period is typically much more
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
  dy <- data$fitted[, what]

  dy_extra <- data$full[, what]
  dy_extra[!is.na(dy)] <- NA_integer_

  par(mgp = c(1.7, 0.5, 0), bty = "n")
  xlim <- c(date_min %||% as.Date("2020-03-16"),
            if (with_forecast) max(x_forecast) else date)

  rn_res <- as.Date(colnames(res))
  ylim <- c(0, max(dy[which(dx >= xlim[1] & dx <= xlim[2])],
                   qs[, which(rn_res >= xlim[1] & rn_res <= xlim[2])],
                   na.rm = TRUE))
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       xlab = "", ylab = labs[what])
  now_cols <- add_alpha(rep(cols$now, 2), alpha)
  fore_cols <- add_alpha(rep(cols$forecast, 2), alpha)

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
    abline(v = as.Date(beta_date), lty = 2, col = cols$puce)
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
           btx = "n", ...)
  }
}
