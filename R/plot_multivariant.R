## Plotting function for multivariant Rt
##'
##' @title Create multivariant Rt plot
##'
##' @param dat Main fitting outputs object
##'
##' @param date Date of the analysis
##'
##' @param date_restart Date for multivariant model restart
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
spim_multivariant_rt_plot <- function(dat, date, date_restart,
                                      last_beta_days_ago = 21,
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
                       angle = 45, size = 3) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(paste(stringr::str_to_sentence(region))) +
    ggplot2::scale_x_date(date_breaks = "months",
                          date_labels = "%b") +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.title = ggplot2::element_blank()) +
    ggthemes::scale_colour_colorblind() +
    ggthemes::scale_fill_colorblind()
  p
}


## Plotting function for fitted variant seeding date
##'
##' @title Create variant seeding date plot
##'
##' @param dat Main fitting outputs object
##'
##' @return A ggplot2 object for fitted seeding date by region
##'
##' @export
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
                   plot.title = ggplot2::element_text(size = 10))
}


## Plotting function for fitted variant transmission advantage
##'
##' @title Create variant transmission advantage plot
##'
##' @param dat Main fitting outputs object
##'
##' @return A ggplot2 object for fitted seeding date by region
##'
##' @export
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
