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
##' @return A ggplot2 object for fit to VOC proportion data
##'
##' @export
spim_plot_vaccine_figure_1 <- function(dat, date, date_restart){

  rt <- spim_multivariant_rt_plot(dat, date, last_beta_days_ago = 8,
                                  rt_type = "eff_Rt_general")

  sd <- spim_plot_seeding_date(dat)

  library(patchwork)
  row1 <- rt + sd + plot_layout(widths = c(3, 1))
  row2 <- spim_plot_voc_proportion(dat, date_restart, "south_east") +
    spim_plot_voc_proportion(dat, date_restart, "south_west") +
    spim_plot_voc_proportion(dat, date_restart, "london") +
    spim_plot_voc_proportion(dat, date_restart, "north_west") +
    plot_layout(ncol = 4, nrow = 1)
  row3 <- spim_plot_voc_proportion(dat, date_restart, "midlands") +
    spim_plot_voc_proportion(dat, date_restart, "east_of_england") +
    spim_plot_voc_proportion(dat, date_restart, "north_east_and_yorkshire") +
    plot_spacer() +
    plot_layout(ncol = 4, nrow = 1)

  g <- row1 / row2 / row3 +
    plot_layout(heights = c(2, 1, 1)) +
    plot_annotation(tag_levels = 'A')
  g
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
##' @return A ggplot2 object for multivariant Rt
##'
##' @export
spim_multivariant_rt_plot <- function(dat, date, last_beta_days_ago = 21,
                                      region = "england",
                                      rt_type = "eff_Rt_general") {
  # Get relevant betas to current date and filter out school holidays
  betas <- data.frame(
    dates = as.Date(tail(dat$samples[[1]]$info$beta_date, 10)),
    label = c(
      "End of 2nd Lockdown",
      "School Holidays",
      "End of December NPI Relaxation",
      "Start of 3rd Lockdown",
      "Roadmap Step 1",
      "School Holidays",
      "Roadmap Step 2",
      "Roadmap Step 3",
      "Original planned date of Step 4",
      "Euro 2020 Final")
  ) %>% dplyr::filter(!stringr::str_detect(label, "School"))

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

  rt_region <- dplyr::left_join(wt, mv) %>% dplyr::left_join(., betas)
  rt_region <- rt_region %>% dplyr::filter(dates >= as.Date("2020-12-01") &
                                             dates <= as.Date(date))

  # Only plot variant after date of first reported cases on 2021-03-23
  variant_names <- c("rt_variant", "lb_variant", "ub_variant")
  rt_region[which(as.Date(rt_region$dates) < "2021-03-23"),
            variant_names] <- NA_integer_

  if (rt_type == "eff_Rt_general") {
    ylim <- c(0, 2)
    ylab <- paste("Effective", expression(R[t]))
  } else if (rt_type == "Rt_general") {
    ylim <- c(0, 4)
    ylab <- paste(expression(R[t]), "excluding immunity")
  }

  p <- ggplot2::ggplot(rt_region, ggplot2::aes(x = dates)) +
    ggplot2::annotate(geom = "rect",
                      xmin = as.Date("2020-12-18"),
                      xmax = as.Date("2021-01-05"),
                      ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
    ggplot2::annotate(geom = "rect",
                      xmin = as.Date("2021-04-02"),
                      xmax = as.Date("2021-04-18"),
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
    ggplot2::geom_text(ggplot2::aes(label = label, y = 0.25),
                       vjust = -0.1, hjust = 0.3,
                       angle = 45, size = 3, family = "Times New Roman") +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(paste(stringr::str_to_sentence(region))) +
    ggplot2::scale_x_date(date_breaks = "months",
                          date_labels = "%b") +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   text = ggplot2::element_text(family = "Times New Roman", size=10),
                   legend.title = ggplot2::element_blank()) +
    ggthemes::scale_colour_colorblind() +
    ggthemes::scale_fill_colorblind()
  p
}


spim_plot_seeding_date <- function(dat) {

  x <- function(x) sircovid::sircovid_date_as_date(as.numeric(x))
  y <- function(x) stringr::str_to_title(stringr::str_replace_all(x, "_", " "))

  regions <- sircovid::regions("england")
  samples <- dat$samples[regions]
  seed_date <- NULL
  for (i in names(samples)) {
    out <- c(i,
             round(mean(samples[[i]]$pars[, "strain_seed_date"])),
             round(quantile(samples[[i]]$pars[, "strain_seed_date"], 0.025)),
             round(quantile(samples[[i]]$pars[, "strain_seed_date"], 0.975))
    )
    seed_date <- rbind(seed_date, out)
  }
  seed_date <- as.data.frame(seed_date)
  rownames(seed_date) <- NULL
  colnames(seed_date) <- c("region", "mean", "lb", "ub")
  seed_date <- seed_date[order(seed_date$mean, decreasing = TRUE), ]
  ylabs <- y(seed_date$region)
  seed_date$region <- factor(seed_date$region,
                             levels = seed_date$region[order(as.numeric(
                               seed_date$mean), decreasing = TRUE)])
  title <- paste("Delta variant seeding date")
  seed_date$mean <- x(seed_date$mean)
  seed_date$lb <- x(seed_date$lb)
  seed_date$ub <- x(seed_date$ub)

  seed_date %>%
    ggplot2::ggplot(ggplot2::aes(y = region, col = region)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lb, xmax = ub)) +
    ggplot2::geom_point(ggplot2::aes(x = mean, col = region, fill = region),
                        shape = 23, size = 3) +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::scale_y_discrete(labels = ylabs) +
    ggplot2::scale_x_date(date_minor_breaks = "1 day",
                          date_breaks = "3 days",
                          date_labels = "%b %d") +
    ggplot2::theme_bw() +
    ggsci::scale_color_lancet() + ggsci::scale_fill_lancet() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
                   text = ggplot2::element_text(family = "Times New Roman", size=10),
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
                        shape = 23, size = 3) +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::scale_y_discrete(labels = ylabs) +
    ggplot2::theme_bw() +
    ggsci::scale_color_lancet() + ggsci::scale_fill_lancet() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.position = "none",
                   text = ggplot2::element_text(family = "Times New Roman", size=10),
                   plot.title = ggplot2::element_text(size = 10))
}


spim_plot_voc_proportion <- function(dat, date_restart, region) {
  sample <- dat$samples[[region]]
  data <- dat$data[[region]]
  date <- dat$info$date

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

  res <- data.frame(
    dates = sircovid::sircovid_date_as_date(sample$trajectories$date[-1]),
    tot_mean = colMeans(
      trajectories["sympt_cases_inc", , -1], na.rm = TRUE),
    tot_lb = matrixStats::colQuantiles(
      trajectories["sympt_cases_inc", , -1], na.rm = TRUE, probs = 0.025),
    tot_ub = matrixStats::colQuantiles(
      trajectories["sympt_cases_inc", , -1], na.rm = TRUE, probs = 0.975),
    neg_mean = colMeans(
      trajectories["sympt_cases_non_variant_inc", , -1], na.rm = TRUE),
    neg_lb = matrixStats::colQuantiles(
      trajectories["sympt_cases_non_variant_inc", , -1],
      na.rm = TRUE, probs = 0.025),
    neg_ub = matrixStats::colQuantiles(
      trajectories["sympt_cases_non_variant_inc", , -1],
      na.rm = TRUE, probs = 0.975)
  ) %>%
    dplyr::filter(dates >= min(df$dates) & dates <= max(df$dates)) %>%
    dplyr::mutate(pos_mean = (tot_mean - neg_mean) / tot_mean * 100) %>%
    dplyr::mutate(pos_lb = (tot_lb - neg_lb) / tot_mean * 100) %>%
    dplyr::mutate(pos_ub = (tot_ub - neg_ub) / tot_mean * 100) %>%
    dplyr::mutate(pos_lb = replace(pos_lb, which(pos_lb < 0), 0)) %>%
    dplyr::mutate(pos_ub = replace(pos_ub, which(pos_ub > 100), 100))

  y <- function(x) stringr::str_to_title(stringr::str_replace_all(x, "_", " "))

  g <- ggplot2::ggplot(df, ggplot2::aes(x = dates)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = res$pos_lb, ymax = res$pos_ub),
                         fill = "blue4", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = res$pos_mean), color = "blue4") +
    ggplot2::geom_point(ggplot2::aes(y = PointEst),
                        color = "chocolate3", size = 1) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           color = "chocolate3", width = 1, size = 0.2) +
    ggplot2::ylab("VOC proportion (%)") +
    ggplot2::xlab("") +
    ggplot2::ggtitle(paste(y(region))) +
    ggplot2::scale_x_date(date_breaks = "months",
                          date_labels = "%b") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 10),
                   text = ggplot2::element_text(family = "Times New Roman",
                                                size=10))
  g
}

