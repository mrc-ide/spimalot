##' Plot trajectories for a region
##'
##' @title Plot trajectories
##'
##' @param region Name of the region
##'
##' @param dat Combined data set
##'
##' @param xlim X-limits (must be dates)
##'
##' @param what Vector of names of traces to plot
##'
##' @param with_forecast Logical, add forecast?
##'
##' @param add_betas Logical, add betas?
##'
##' @return Nothing, called for side effects
##'
##' @export
spim_plot_trajectories_region <- function(region, dat, xlim, what = NULL,
                                          with_forecast = TRUE,
                                          add_betas = FALSE) {
  if (is.null(what)) {
    what <- c("deaths_hosp", "deaths_comm", "deaths_carehomes",
              "icu", "general", "admitted", "diagnoses")
  }
  for (w in what) {
    spim_plot_trajectories_region1(
      w, region, dat, xlim,
      with_forecast = with_forecast, add_betas = add_betas)
  }
}


plot_Rt <- function(dat, rt_type, forecast_until = NULL) {
  oo <- par(mfrow = c(2, 6), oma = c(2, 1, 2, 1), mar = c(3, 3, 3, 1))
  on.exit(par(oo))
  if (is.null(forecast_until)) {
    forecast_until <- dat$info$date
  }
  for (r in names(dat$rt)) {
    plot_Rt_region(r, dat, rt_type, forecast_until)
  }
}


plot_Rt_region <- function(region, dat, rt_type, forecast_until) {
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
  qs <- qs[,idx_not_na]

  oo <- par(mgp = c (1.7, 0.5, 0), bty = "n")
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


spim_plot_trajectories_region1 <- function(what, region, dat, xlim,
                                           with_forecast = TRUE,
                                           add_betas = FALSE) {
  trajectories <- dat$samples[[region]]$trajectories
  date <- as.Date(dat$info$date)
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
  res_nowcast <- res[ , as.character(x_nowcast)]
  res_forecast <- res[ , as.character(x_forecast)]

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
  if (with_forecast) {
    xlim[2] <- max(x_forecast)
  }

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

  points(dx, dy, pch =23, bg = cols$orange, col = cols$brown, cex = 0.7,
         lwd = 0.6)
  points(dx, dy_extra, pch =23, bg = cols$green2, col = cols$brown, cex = 0.7,
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
    if(horiz) {
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
           pch =15,
           col = cols[leg_cols],
           legend = leg,
           border = NA,
           btx = "n", ... )
  }
}
