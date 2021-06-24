##' Plot Rt distribution over time for various scenarios
##' @title Plot Rt distribution
##'
##' @param npi_key data.frame type object with columns:
##'   nations, npi (scenario), Rt (mean), Rt_sd
##' @param xlim Passed to `plot`
##' @param ylim Passed to `plot`
##' @param labels Passed to [legend] `legend` parameter
##' @param legend_ncol Passed to [legend] `ncol` parameter
##' @param npi_key2 Optional second npi_key that will be plotted with dashed
##'   lines
##'
##' @export
spim_plot_rt_dist <- function(npi_key, xlim, ylim, labels = NULL,
                              legend_ncol = 1, npi_key2 = NULL) {

  labels <- labels %||% rownames(npi_key)
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n",
       xlab = expression(R[excl_immunity]),
       ylab = "Density",
       las = 1)
  cols <- viridis(nrow(npi_key))
  for (i in seq_rows(npi_key)) {
    dist <- distr6::dstr("Lognormal", mean = npi_key$Rt[i], sd = npi_key$Rt_sd[i])
    curve(dist$pdf(x), add = TRUE, col = cols[i], from = xlim[1], to = xlim[2],
          n = 1e3, lwd = 2)
    legend_lty <- rep(1, nrow(npi_key2))
  }
  if (!is.null(npi_key2)) {
    cols <- viridis(nrow(npi_key2))
    for (i in seq_rows(npi_key)) {
      dist <- distr6::dstr("Lognormal", mean = npi_key2$Rt[i], sd = npi_key2$Rt_sd[i])
      curve(dist$pdf(x),
            add = TRUE, col = cols[i], from = xlim[1], to = xlim[2],
            n = 1e3, lwd = 2, lty = 2
      )
    }
    legend_lty <- c(legend_lty, 1, 2)
    cols <- c(cols, 1, 1)
  }
  legend("topleft", legend = labels, ncol = legend_ncol, bty = "n",
         lty = legend_lty, lwd = 2, col = cols)
}


##' Plot seasonality trends over time
##' @title Plot sesonality over time
##'
##' @param peak_date Date where seasonal multiplier is highest
##' @param seasonality Seasonal multiplier
##'
##' @export
spim_plot_seasonality <- function(peak_date = as.Date("2020-02-15"),
                                  seasonality = 0.1) {
  x <- seq_len(365)
  dx <- sircovid::sircovid_date_as_date(x)
  plot(dx,
       calc_seasonality(x, sircovid::sircovid_date(peak_date), seasonality),
       type = "l", lwd = 2, ylim = c(0.9, 1.1),
       ylab = "Seasonal multiplier", xlab = "", xaxt = "n")
  axis.Date(1, dx, at = seq.Date(from = dx[1], to = as.Date("2021-01-01"),
                                 by = "1 month"))
  abline(v = c(peak_date, round(peak_date + 365 / 2)), lty = 2, col = "grey30")
}


##' Plot relative VOC mean and sd
##' @title Plot relative strain transmissibility ranges
##'
##' @param R1 Lognormal mean for strain 1
##' @param R1_sd Lognormal standard deviation for strain 1
##' @param epsilon_range Relative range for strain 2
##' @param epsilon_central Relative mean for strain 2
##'
##' @export
spim_plot_voc_range <- function(R1, R1_sd, epsilon_range, epsilon_central) {

    R1_range <- as.list(
        distr6::dstrs(
            "Lognormal",
            data.frame(mean = R1, sd = R1_sd)
        )$quantile(c(0.025, 0.975))
    )

    R2_range <- lapply(seq_along(R1_range), function(i) {
        c(min(R1_range[[i]]) * min(epsilon_range), max(R1_range[[i]]) * max(epsilon_range))
    })

    data_frame(do.call(rbind, Map(rbind, R1_range, R2_range))) %>%
       `colnames<-`(c("min", "max")) %>%
        dplyr::mutate(
            central = c(R1[1], (R1 * epsilon_central)[1], R1[2],
                        (R1 * epsilon_central)[2]),
            variant = rep(c("B.1.1.7", "B.1.617.2"), 2),
            NPI = c(rep("central R after NPI lift", 2),
                    rep("high R after NPI lift", 2)),
            NPI = factor(NPI, levels = unique(NPI))
        ) %>%
        ggplot2::ggplot() +
        geom_pointrange(aes(x = NPI, y = central, ymin = min, ymax = max,
                            colour = variant),
            size = 1, position = position_dodge(width = 0.2)
        ) +
        scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 11)) +
        theme_classic() +
        xlab(NULL) +
        ylab("R excluding immunity") +
        theme(
            text = element_text(size = 14), legend.position = c(0.2, 0.9),
            legend.title = element_blank()
        )
}

##' Select colours for plotting simulation scenarios
##' @title Return accessible scenario colours
##'
##' @param scenarios Unique scenario names
##' @param dark_scenarios Optional unique scenario names that should be same
##'  length as `scenarios` as provided and will be the same colours but darker.
##'  Useful if plotting scenarios and their High R counterparts.
##' @param weight If `dark_scenarios` is not `NULL` then `weight` passed to
##'  `mix_cols` to darken colours by mixing with "#000000" (black)
##' @param palette Colour palette, passed to [khroma::colour]
##' @param highR If `TRUE` then `dark_scenarios` is taken to be all scenarios
##'   that contain "High R"; `dark_scenarios` should be `NULL` if `TRUE`
##' @param preview If `TRUE` then plots the final colour scheme with
##'   [khroma::plot_scheme]
##'
##' @export
spim_scenario_cols <- function(scenarios, dark_scenarios = NULL, weight = 0.3,
                               palette = "bright", highR = TRUE,
                               preview = FALSE) {

  stopifnot(all(table(scenarios)) == 1)

  n_scens <- length(scenarios)

  if (highR) {
    if (!is.null(dark_scenarios)) {
      stop("`dark_scenarios` must be `NULL` if `highR` is `TRUE`")
    }
    dark_scenarios <- grep("High R", scenarios, value = TRUE)
    scenarios <- setdiff(scenarios, dark_scenarios)
    n_scens <- length(scenarios)
    if (length(dark_scenarios) == 0) {
      dark_scenarios <- NULL
    }
  }

  if (!is.null(dark_scenarios)) {
    stopifnot(all(table(dark_scenarios)) == 1,
              n_scens == length(dark_scenarios))
  }

  cols <- khroma::colour(palette)(n_scens)
  names(cols) <- scenarios

  if (!is.null(dark_scenarios)) {
    dark_cols <- mix_cols(cols, rep("#000000", length(cols)), 0.3)
    names(dark_cols) <- dark_scenarios
    cols <- c(cols, dark_cols)
  }

  if (preview) {
    khroma::plot_scheme(cols)
  }

  cols
}