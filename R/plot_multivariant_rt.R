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
                                      rt_type = "eff_Rt_general"){
  # Get relevant betas to current date and filter out school holidays
  betas <- data.frame(
    dates = tail(dat$samples[[1]]$info$beta_date, 10),
    label = c(
      "End of 2nd Lockdown",
      "School Holidays",
      "End of December NPI Relaxation",
      "Start of 3rd Lockdown",
      "Roadmap Step 1",
      "School Holidays",
      "Roadmap Step 2",
      "Roadmap Step 3",
      "Date of Step 4 Delay Announcement",
      "Euros 2021 Final")
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
    ylim <- c(0,2)
    ylab <- paste("Effective", expression(R[t]))
  } else if (rt_type == "Rt_general") {
    ylim <- c(0,4)
    ylab <- paste(expression(R[t]), "excluding immunity")
  }

  p <- ggplot2::ggplot(rt_region, ggplot2::aes(x = dates)) +
    ggplot2::annotate(geom = "rect",
                      xmin = as.Date('2020-12-18'),
                      xmax = as.Date('2021-01-05'),
                      ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4) +
    ggplot2::annotate(geom = "rect",
                      xmin = as.Date('2021-04-02'),
                      xmax = as.Date('2021-04-18'),
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
    ggplot2::geom_vline(xintercept = as.Date(betas[ ,1]), lty = 3,
                        col = "red4") +
    ggplot2::geom_text(ggplot2::aes(label = label, y = 0.25),
                       angle = 45, size = 3) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(paste(stringr::str_to_sentence(region))) +
    ggplot2::scale_x_date(date_breaks = "months",
                          date_labels = "%b") +
    ggplot2::theme_bw() + ggplot2::ylim(ylim) +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(),
                   legend.title = ggplot2::element_blank()) +
    ggthemes::scale_colour_colorblind() +
    ggthemes::scale_fill_colorblind()
  p
}