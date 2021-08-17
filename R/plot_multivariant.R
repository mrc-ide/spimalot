## Plotting function for vaccine paper figure 1
##'
##' @title Create plot for vaccine paper figure 1
##'
##' @param dat Main fitting outputs object
##'
##' @param date Date of the analysis
##'
##' @param date_restart Date for multivariant model restart
##'
##' @param manuscript Logical indicating whether plots are for manuscript
##'
##' @return A ggplot2 object for fit to VOC proportion data
##'
##' @param manuscript Logical indicating if plotting for manuscript
##'
##' @export
spim_plot_vaccine_figures <- function(dat, date, date_restart,
                                       manuscript = TRUE){


  out <- list()
  out$rt <- spim_multivariant_rt_plot(dat, date, last_beta_days_ago = 8,
                                      rt_type = "eff_Rt_general",
                                      manuscript = manuscript) +
    ggplot2::ggtitle(NULL)

  out$sd <- spim_plot_seeding_date(dat) + ggplot2::ggtitle(NULL)

  for (i in sircovid::regions("england")) {
    out[[i]] <- spim_plot_voc_proportion(dat, date_restart, i)
  }

  row1 <- out$rt
  row2 <- out$sd + out$london +
    patchwork::plot_layout(ncol = 2, nrow = 1)


  g <- row1 / row2 +
    patchwork::plot_layout(heights = c(2, 1), guides = "keep") +
    patchwork::plot_annotation(tag_levels = 'A')

  out$fig_1 <-
    g & ggplot2::theme(plot.margin = ggplot2::unit(rep(1, 4), units = "mm"))

  out
}


## Plotting function for multivariant Rt
##'
##' @title Create multivariant Rt plot
##'
##' @param dat Main fitting outputs object
##'
##' @param date Date of the analysis
##'
##' @param last_beta_days_ago Integer for setting last beta changepoint in
##'   relation to date of analysis
##'
##' @param region A string for the name of the region to plot
##'
##' @param rt_type A string, must be one of eff_Rt_general or Rt_general
##'
##' @param manuscript Logical indicating whether plots are for manuscript
##'
##' @return A ggplot2 object for multivariant Rt
##'
##' @export
spim_multivariant_rt_plot <- function(dat, date, last_beta_days_ago = 21,
                                      region = "england",
                                      rt_type = "eff_Rt_general",
                                      manuscript = FALSE) {
  # Get relevant betas to current date and filter out school holidays
  betas <- data.frame(
    dates = as.Date(tail(dat$samples[[1]]$info$beta_date, 12)),
    label = c(
      "End of 2nd\nLockdown",
      "School Holidays",
      "Holiday \nRestrictions",
      "Start of 3rd\nLockdown",
      "Roadmap\nStep 1",
      "School Holidays",
      "Roadmap\nStep 2",
      "Roadmap\nStep 3",
      "Step 4\nDelayed",
      "Euro 2020\nQtr Final",
      "Euro 2020\nFinal",
      "Roadmap\nStep 4"),
    label_y = c(1.94, NA, 1.67, 1.94, 1.67, NA,
                1.94, 1.67, 1.94, 1.67, 1.94, 1.67)
  ) %>%
    dplyr::filter(!stringr::str_detect(label, "School"))

  if (!manuscript) {
    betas$label_y <- NA_integer_
  }

  rt_plot <- NULL

  mv <- data.frame(
    dates = sircovid::sircovid_date_as_date(
      dat$variant_rt[[region]][["date"]][-1, 1]),
    rt_variant = rowMeans(dat$variant_rt[[region]][[rt_type]][-1, 2, ]),
    lb_variant = matrixStats::rowQuantiles(
      dat$variant_rt[[region]][[rt_type]][-1, 2, ], probs = 0.025),
    ub_variant = matrixStats::rowQuantiles(
      dat$variant_rt[[region]][[rt_type]][-1, 2, ], probs = 0.975),
    rt_non_variant = rowMeans(dat$variant_rt[[region]][[rt_type]][-1, 1, ]),
    lb_non_variant = matrixStats::rowQuantiles(
      dat$variant_rt[[region]][[rt_type]][-1, 1, ], probs = 0.025),
    ub_non_variant = matrixStats::rowQuantiles(
      dat$variant_rt[[region]][[rt_type]][-1, 1, ], probs = 0.975)
  )
  wt <- data.frame(
    dates = sircovid::sircovid_date_as_date(dat$rt[[region]][["date"]][-1, 1]),
    rt_weighted = rowMeans(dat$rt[[region]][[rt_type]][-1, ]),
    lb_weighted = matrixStats::rowQuantiles(
      dat$rt[[region]][[rt_type]][-1, ], probs = 0.025),
    ub_weighted = matrixStats::rowQuantiles(
      dat$rt[[region]][[rt_type]][-1, ], probs = 0.975)
  )

  if (manuscript) {
    date_end <- "2021-07-19"
  } else {
    date_end <- date
  }

  rt_region <- dplyr::left_join(wt, mv) %>% dplyr::left_join(., betas)
  rt_region <- rt_region %>% dplyr::filter(dates >= as.Date("2020-12-01") &
                                             dates <= as.Date(date_end))

  # Only plot variant after date of first reported cases on 2021-03-23
  variant_names <- c("rt_variant", "lb_variant", "ub_variant")
  rt_region[which(as.Date(rt_region$dates) < "2021-03-23"),
            variant_names] <- NA_integer_

  if (rt_type == "eff_Rt_general") {
    ylim <- c(0, 2)
    ylab <- "Effective R(t)"
  } else if (rt_type == "Rt_general") {
    ylim <- c(0, 4)
    ylab <- "R(t) excluding immunity"
  }

  p <- ggplot2::ggplot(rt_region, ggplot2::aes(x = dates)) +
    ggplot2::annotate(geom = "rect",
                      xmin = as.Date("2020-12-18"),
                      xmax = as.Date("2021-01-05"),
                      ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
    ggplot2::annotate(geom = "rect",
                      xmin = as.Date("2021-04-01"),
                      xmax = as.Date("2021-04-19"),
                      ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = rt_weighted, col = "Weighted")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lb_weighted, ymax = ub_weighted,
                                      fill = "Weighted"),
                         alpha = 0.3, show.legend = FALSE) +
    ggplot2::geom_line(ggplot2::aes(y = rt_variant, col = "Delta")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lb_variant,
                                      ymax = ub_variant, fill = "Delta"),
                         alpha = 0.3, show.legend = FALSE) +
    ggplot2::geom_line(ggplot2::aes(y = rt_non_variant, col = "Alpha")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lb_non_variant,
                                      ymax = ub_non_variant, fill = "Alpha"),
                         alpha = 0.3, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 1, lty = 2, col = "black") +
    ggplot2::geom_vline(xintercept = as.Date(betas[, 1]), lty = 3,
                        col = "red4") +
    ggplot2::geom_label(ggplot2::aes(label = label, y = label_y),
                        hjust = 0.5, size = 3,
                        vjust = 0.5, family = "serif",
                        label.padding = ggplot2::unit(0.15, "lines")) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(paste(stringr::str_to_sentence(region))) +
    ggplot2::scale_x_date(date_breaks = "months",
                          date_labels = "%b") +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   text = ggplot2::element_text(
                     family = "serif", size = 10),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.box.margin = ggplot2::margin(
                     -10, 0, 0, 0, unit = "mm"),
                   plot.margin = ggplot2::unit(rep(0, 4), units = "cm")) +
    ggthemes::scale_colour_colorblind() +
    ggthemes::scale_fill_colorblind()
  p
}


spim_plot_seeding_date <- function(dat) {

  regions <- sircovid::regions("england")
  ylabs <- c("EE", "MID", "LON", "NEY", "NW", "SE", "SW")
  samples <- dat$samples[regions]


  strain_seed_date <- sapply(samples, function(x) x$pars[, "strain_seed_date"])

  seed_date <- list(mean = colMeans(strain_seed_date),
                    lb = apply(strain_seed_date, 2, quantile, 0.025),
                    ub = apply(strain_seed_date, 2, quantile, 0.975))
  seed_date <- as.data.frame(lapply(seed_date, sircovid::sircovid_date_as_date))
  seed_date$regions <- regions

  seed_date$region <- factor(ylabs,
                             levels = ylabs[order(seed_date$mean,
                                                  decreasing = TRUE)])

  seed_date %>%
    ggplot2::ggplot(ggplot2::aes(y = region, col = region)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lb, xmax = ub)) +
    ggplot2::geom_point(ggplot2::aes(x = mean, col = region, fill = region),
                        shape = 23, size = 1) +
    ggplot2::labs(x = "", y = "", title = "Delta seeding date") +
    ggplot2::scale_y_discrete() +
    ggplot2::scale_x_date(date_minor_breaks = "1 day",
                          date_breaks = "3 days",
                          date_labels = "%b %d") +
    ggplot2::theme_minimal() +
    ggsci::scale_color_lancet() + ggsci::scale_fill_lancet() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
                   text = ggplot2::element_text(
                     family = "Times New Roman", size = 10),
                   plot.title = ggplot2::element_text(size = 10))
}


spim_plot_variant_transmission <- function(dat) {

  x <- function(x) as.numeric(x)
  y <- function(x) stringr::str_to_title(stringr::str_replace_all(x, "_", " "))

  regions <- sircovid::regions("england")
  samples <- dat$samples[regions]
  variant_trans <- NULL
  for (i in names(samples)) {
    out <- c(i,
             mean(samples[[i]]$pars[, "strain_transmission_2"]),
             quantile(samples[[i]]$pars[, "strain_transmission_2"], 0.025),
             quantile(samples[[i]]$pars[, "strain_transmission_2"], 0.975)
    )
    variant_trans <- rbind(variant_trans, out)
  }
  variant_trans <- as.data.frame(variant_trans)
  rownames(variant_trans) <- NULL
  colnames(variant_trans) <- c("region", "mean", "lb", "ub")
  variant_trans <- variant_trans[order(variant_trans$mean, decreasing = TRUE), ]
  ylabs <- y(variant_trans$region)
  variant_trans$region <- factor(variant_trans$region,
                                 levels = variant_trans$region[order(as.numeric(
                                   variant_trans$mean), decreasing = TRUE)])
  title <- paste("Delta transmission advantage")
  variant_trans$mean <- x(variant_trans$mean)
  variant_trans$lb <- x(variant_trans$lb)
  variant_trans$ub <- x(variant_trans$ub)

  variant_trans %>%
    ggplot2::ggplot(ggplot2::aes(y = region, col = region)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lb, xmax = ub)) +
    ggplot2::geom_point(ggplot2::aes(x = mean, col = region, fill = region),
                        shape = 23, size = 1) +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::scale_y_discrete(labels = ylabs) +
    ggplot2::theme_minimal() +
    ggsci::scale_color_lancet() + ggsci::scale_fill_lancet() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.position = "none",
                   text = ggplot2::element_text(family = "Times New Roman",
                                                size = 10),
                   plot.title = ggplot2::element_text(size = 10))
}


spim_plot_voc_proportion <- function(dat, date_restart, region) {
  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date

  region_labs <- c("EE", "MID", "LON", "NEY", "NW", "SE", "SW")
  names(region_labs) <- sircovid::regions("england")

  df <- NULL
  df <- data.frame(
    dates = data$fitted[, "date_string"],
    region = region,
    n_non_variant = data$fitted[, "strain_non_variant"],
    ntot = data$fitted[, "strain_tot"]
  ) %>%
    dplyr::mutate(npos = ntot - n_non_variant) %>%
    dplyr::mutate(npos = tidyr::replace_na(npos, 0)) %>%
    dplyr::mutate(ntot = tidyr::replace_na(ntot, 0)) %>%
    dplyr::filter(dates >= date_restart)

  cis <- Hmisc::binconf(x = df$npos, n = df$ntot) * 100
  df$PointEst <- cis[, "PointEst"]
  df$lower <- cis[, "Lower"]
  df$upper <- cis[, "Upper"]
  df$ntot[df$ntot == 0] <- NA

  trajectories <- sample$trajectories$state
  prop_pos <- (1 - trajectories["sympt_cases_non_variant_inc", , -1] /
                 trajectories["sympt_cases_inc", , -1]) * 100

  res <- data.frame(
    dates = sircovid::sircovid_date_as_date(sample$trajectories$date[-1]),
    pos_mean = colMeans(prop_pos, na.rm = TRUE),
    pos_lb = apply(prop_pos, 2, quantile, probs = 0.025, na.rm = TRUE),
    pos_ub = apply(prop_pos, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
    dplyr::filter(dates >= min(df$dates), dates <= max(df$dates))

  g <- ggplot2::ggplot(df, ggplot2::aes(x = dates)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = res$pos_lb, ymax = res$pos_ub),
                         fill = "blue4", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = res$pos_mean), color = "blue4",
                       lwd = 0.5) +
    ggplot2::geom_point(ggplot2::aes(y = PointEst),
                        color = "grey30", size = 0.1, shape = 23) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           color = "grey30", width = 1, size = 0.1) +
    ggplot2::ylab("Delta (%)") +
    ggplot2::xlab("") +
    ggplot2::ggtitle(region_labs[region]) +
    ggplot2::scale_x_date(date_breaks = "months",
                          date_labels = "%b") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 10),
                   text = ggplot2::element_text(family = "Times New Roman",
                                                size = 10))
  g
}
