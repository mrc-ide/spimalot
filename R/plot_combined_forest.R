##' Forest plot
##'
##' @title Forest plot
##' @param dat Combined data from [spimalot::spim_combined_load]
##'
##' @param regions Vector of regions to plot. By default, all regions
##'   found in `dat` will be used (except aggregate regions, such as
##'   england/uk)
##'
##' @param plot_type The type of parameters to plot: `all` would plot
##'   all parameters, `betas` would just plot betas, `non_betas` would
##'   plot all non-beta parameters, `subset` allows specifying a subset
##'   of the parameters
##'
##' @param subset A vector of parameter names to plot if `subset` is selected
##'   for `plot_type`. Otherwise can just be set to NULL (the default)
##'
##' @param par_labels A list of labels to use for parameters,
##'   any parameter not listed will just have a label matching its name
##'
##' @return Nothing, called for side effects
##' @export
spim_plot_forest <- function(dat, regions = NULL, plot_type = "all",
                             subset = NULL, par_labels = list()) {
  if (is.null(regions)) {
    regions <- intersect(sircovid::regions("all"), names(dat$samples))
  } else {
    msg <- setdiff(regions, names(dat$samples))
    if (length(msg) > 0) {
      stop("regions missing from 'dat': ", paste(squote(msg), collapse = ", "))
    }
  }

  match_value(plot_type, c("all", "betas", "non_betas", "subset"))

  samples <- dat$samples[regions]
  date <- dat$info$date
  model_type <- dat$info$model_type

  countries <- setdiff(regions, sircovid::regions("england"))

  region_names <- setdiff(regions, sircovid::regions("nations"))

  if ("start_date" %in% colnames(samples[[1]]$pars)) {
    # order by start_date
    mean_start_date <-
      vnapply(samples, function(x) mean(x$pars[, "start_date"]))

    ordered_regions <- c(names(c(sort(mean_start_date[countries],
                                      decreasing = TRUE),
                                 sort(mean_start_date[region_names],
                                      decreasing = TRUE))))
  } else {
    ordered_regions <- c(sort(countries, decreasing = TRUE),
                         sort(region_names, decreasing = TRUE))
  }

  samples <- samples[ordered_regions]

  # formats for region labels
  labels <- spim_region_name(names(samples), "code")
  par_names <- colnames(samples[[1]]$pars)

  beta_date <- sircovid::sircovid_date_as_date(dat$samples[[1]]$info$beta_date)
  beta_names <- par_names[substr(par_names, 1, 4) == "beta"]
  beta_names <- beta_names[order(as.numeric(gsub("beta", "", beta_names)))]

  ## Can set the xmax for specific parameters here. For any left as NA, it will
  ## instead just use the maximum from the parameters info
  par_max <- rep(NA, length(par_names))
  names(par_max) <- par_names
  par_max[c("m_CHW", "m_CHR")] <- 2e-5
  par_max[beta_names] <- 0.25
  if (model_type == "BB") {
    pars_p_NC <- grep("p_NC", names(par_max), value = TRUE)
    par_max[pars_p_NC] <- 0.01
    par_max["rho_pillar2_tests"] <- 0.02
  }

  n_regions <- length(samples)

  extract_sample <- function(par_name) {
    lapply(samples, function(x) as.numeric(x$pars[, par_name]))
  }

  ylim <- c(0.5, n_regions + 0.5)

  plot_axis <- function() {
    yax_pos <- 0.9
    plot(0, 0, type = "n",
         ylab = "",
         xlab = "",
         xlim = c(0, 0),
         ylim = ylim,
         axes = FALSE
    )
    axis(side = 2, at =  at, labels = labels, las = 1, pos = yax_pos)
  }

  col_line <- "red3" # TODO: move into colours

  plot_par <- function(par_name) {
    if (par_name == "start_date") {
      plot_start_date()
    } else {
      par <- extract_sample(par_name)

      par_info <- subset(pars_info, pars_info$name == par_name)
      if (is.na(par_max[[par_name]])) {
        xmax <- max(par_info$max)
      } else {
        xmax <- par_max[[par_name]]
      }

      xmin <- min(par_info$min)

      if (par_name %in% names(par_labels)) {
        xlab <- par_labels[[par_name]]
      } else {
        if (grepl("^beta", par_name)) {
          k <- as.numeric(gsub("beta", "", par_name))
          xlab <- paste0(par_name, " (", beta_date[k], ")")
        } else {
          xlab <- par_name
        }
      }

      plot(0, 0, type = "n",
           ylab = "",
           xlab = xlab,
           xlim = c(xmin, xmax),
           ylim = ylim,
           yaxt = "n"
      )

      jitter <- 0.5
      regions <- names(par)
      hp <- subset(hps, hps$name == par_name)
      if (is.na(hp$region[1])) {
        hp <- hp[rep(1, length(regions)), ]
        hp$region <- regions
        col <- "red"
      } else {
        col <- "grey20"
      }
      rownames(hp) <- hp$region
      hp <- hp[regions, ] # sort in correct order
      if (hp$type[1] == "beta") {
        shape1 <- hp$beta_shape1
        shape2 <- hp$beta_shape2
        if (!(all(shape1 == 1) && all(shape2 == 1))) {
          prior <- mapply(qbeta,
                          shape1 = shape1,
                          shape2 = shape2,
                          MoreArgs = list(p = c(0.025, 0.975)),
                          SIMPLIFY = TRUE)


          segments(x0 = prior[1, ],
                   y0 = at - jitter, y1 = at + jitter,
                   col = col_line, lty = 2, lwd = 1, lend = 2)
          segments(x0 = prior[2, ],
                   y0 = at - jitter, y1 = at + jitter,
                   col = col_line, lty = 2, lwd = 1, lend = 2)
        }
      }
      if (hp$type[1] == "gamma") {
        shape <- hp$gamma_shape
        scale <- hp$gamma_scale
        prior <- mapply(qgamma,
                        shape = shape,
                        scale = scale,
                        MoreArgs = list(p = c(0.025, 0.975)),
                        SIMPLIFY = TRUE)


        segments(x0 = prior[1, ],
                 y0 = at - jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
        segments(x0 = prior[2, ],
                 y0 = at - jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
      }

      mapply(FUN = plot_ci_bar, res = par, at = at, width = 0.1,
             col = col)
    }
  }

  plot_start_date <- function() {
    numeric_start_date <- extract_sample(par_name = "start_date")
    if ("start_date" %in% names(par_labels)) {
      xlab = par_labels$start_date
    } else {
      xlab = "start_date"
    }
    plot(x = sircovid::sircovid_date("2020-01-01"),
         y = 1,
         type = "n",
         ylab = "",
         xlab = xlab,
         xlim = sircovid::sircovid_date(c("2020-01-15", "2020-03-01")),
         ylim = ylim,
         yaxt = "n")
    mapply(FUN = plot_ci_bar, res = numeric_start_date, at = at, width = 0.1)
  }

  hps <- dat$parameters$prior
  pars_info <- dat$parameters$info

  op <- par(bty = "n", mar = c(3, 0, 1, 0), mgp = c(2, 0.75, 0),
            oma = c(2, 0, 3, 3))
  on.exit(par(op))

  if (plot_type == "all") {
    pars_to_plot <- par_names
  } else if (plot_type == "betas") {
    pars_to_plot <- beta_names
  } else if (plot_type == "non_betas") {
    pars_to_plot <- setdiff(par_names, beta_names)
  } else if (plot_type == "subset") {
    if (is.null(subset)) {
      stop("Expected a 'subset' input as subset plot type has been selected")
    }
    missing_pars <- setdiff(subset, par_names)
    if (length(missing_pars > 0)) {
      stop("The following parameters listed in 'subset' are missing from the
           fitted parameters: ", paste(missing_pars, collapse = ", "))
    }
    pars_to_plot <- subset
  }

  if ("start_date" %in% pars_to_plot) {
    pars_to_plot <- c("start_date", pars_to_plot[pars_to_plot != "start_date"])
  }

  npar <- length(pars_to_plot)
  nrow <- 4
  plot_per_row <- ceiling(npar / nrow)
  colwidth <- 64 / plot_per_row

  reps <- rep(c(3.5, rep(colwidth, plot_per_row)), nrow)

  layout(mat = matrix(rep(seq_along(reps), reps),
                      nrow = nrow, byrow = TRUE),
         heights = rep(3, nrow),
         widths = c(4, rep(1, plot_per_row * colwidth))
  )
  at <- seq_len(n_regions)

  for (i in seq_len(nrow)) {
    plot_axis()
    par(mar = c(3, 0, 1, 0.5))
    pars_row <- seq((i - 1) * plot_per_row + 1, min(npar, i * plot_per_row))
    mapply(plot_par, par_name = pars_to_plot[pars_row])
  }

  mtext(text = paste("Inferred epidemic parameters for NHS regions at", date),
        side = 3, line = 0, outer = TRUE, cex = 1.1)
}


plot_ci_bar <- function(res, at, width = 1,
                        min = 0.025, max = 0.975, col = "grey20",
                        segments = FALSE, pt_col = NULL, horiz = TRUE, ...) {
  cols  <- c("grey80", col)
  qs <- quantile(res,
                 probs = seq(min, max, by = 0.005),
                 na.rm = TRUE)

  palette <- grDevices::colorRampPalette(cols)
  if (segments) {
    segments(y0 = at, x0 = min(res), x1 = max(res), col = cols[2])
    points(y = rep(at, 2), x = range(res), col = cols[2], pch = "-")
  }
  ci_bands(quantiles = cbind(qs, qs),
           y = at + c(-1, 1) * width,
           palette = palette, leg = FALSE, horiz = horiz)

  if (is.null(pt_col)) pt_col <- col
  if (horiz) {
    points(y = at, x = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  } else {
    points(x = at, y = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  }
}
