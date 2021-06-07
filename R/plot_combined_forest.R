##' Forest plot
##'
##' @title Forest plot
##' @param dat Combined data from [spimalot::spim_combined_load]
##'
##' @param plot_betas Logical, indicating if we should plot betas
##'   (rather than all other parameters)
##'
##' @return Nothing, called for side effects
##' @export
spim_plot_forest <- function(dat, plot_betas) {
  samples <- dat$samples[sircovid::regions("all")]
  date <- dat$info$date
  model_type <- dat$info$model_type

  # order by start_date
  mean_start_date <- vnapply(samples, function(x) mean(x$pars[, "start_date"]))

  countries <- c("scotland", "wales", "northern_ireland")
  region_names <- setdiff(names(samples), c(countries, "england", "uk"))

  ordered_regions <- c(names(c(sort(mean_start_date[countries],
                                    decreasing = TRUE),
                               sort(mean_start_date[region_names],
                                    decreasing = TRUE))))
  samples <- samples[ordered_regions]

  # formats for region labels
  labels <- spim_region_name(names(samples), "code")
  par_names <- colnames(samples[[1]]$pars)
  ch_names <- c("eps", "m_CHW", "m_CHR")
  mu_names <- par_names[substr(par_names, 1, 3) == "mu_"]
  alpha_names <- par_names[substr(par_names, 1, 6) == "alpha_"]
  beta_names <- par_names[substr(par_names, 1, 4) == "beta"]
  beta_names <- beta_names[order(as.numeric(gsub("beta", "", beta_names)))]
  p_names <- par_names[substr(par_names, 1, 2) == "p_"]

  p_max <- rep(1, length(p_names))
  if (model_type == "BB") {
    p_max[p_names == "p_NC"] <- 0.01
  }

  n_regions <- length(samples)

  extract_sample <- function(par_name) {
    lapply(samples, function(x) as.numeric(x$pars[, par_name]))
  }

  numeric_start_date <- extract_sample(par_name = "start_date")

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

  plot_par <- function(par_name, xmin = 0, xmax = 0.1) {
    par <- extract_sample(par_name)

    plot(0, 0, type = "n",
         ylab = "",
         xlab = par_name,
         xlim = c(xmin, xmax),
         ylim = ylim,
         yaxt = "n"
    )

    jitter <- 0.5
    regions <- names(par)
    hp <- subset(hps, hps$name == par_name)
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

    mapply(FUN = plot_ci_bar, res = par, at = at, width = 0.1)
  }

  plot_start_date <- function() {
    plot(x = sircovid::sircovid_date("2020-01-01"),
         y = 1,
         type = "n",
         ylab = "",
         xlab = "start_date",
         xlim = sircovid::sircovid_date(c("2020-01-15", "2020-03-01")),
         ylim = ylim,
         yaxt = "n")
    mapply(FUN = plot_ci_bar, res = numeric_start_date, at = at, width = 0.1)
  }

  hps <- dat$parameters$prior

  op <- par(bty = "n", mar = c(3, 0, 1, 0), mgp = c(2, 0.75, 0),
            oma = c(2, 0, 3, 3))
  on.exit(par(op))

  if (plot_betas) {
    npar <- length(beta_names)
  } else {
    npar <- length(par_names) - length(beta_names)
  }
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

  if (plot_betas) {
    ## 1st row
    plot_axis()
    par(mar = c(3, 0, 1, 0.5))
    mapply(FUN = plot_par, par_name = beta_names[1:5], xmax = 0.15)

    ## 2nd row
    plot_axis()
    mapply(FUN = plot_par, par_name = beta_names[6:10], xmax = 0.15)

    ## 3rd row
    plot_axis()
    mapply(FUN = plot_par, par_name = beta_names[11:15], xmax = 0.15)

    ## 4th row
    plot_axis()
    mapply(FUN = plot_par, par_name = beta_names[16:18], xmax = 0.15)

  } else {
    ## 1st row
    plot_axis()
    par(mar = c(3, 0, 1, 0.5))
    plot_start_date()
    mapply(FUN = plot_par, par_name = p_names[1:5], xmax = p_max[1:5])

    if (model_type == "BB") {
      ## 2nd row
      plot_axis()
      mapply(FUN = plot_par, par_name = p_names[6:10], xmax = p_max[6:10])
      plot_par("rho_pillar2_tests", xmax = 0.02)

      ## 3rd row
      plot_axis()
      mapply(FUN = plot_par, par_name = mu_names[1:2], xmax = 1)
      mapply(FUN = plot_par, par_name = mu_names[3:4], xmax = 2)
      mapply(FUN = plot_par, par_name = alpha_names[1:2], xmax = 1)
    } else {
      ## 2nd row
      plot_axis()
      mapply(FUN = plot_par, par_name = p_names[6:9], xmax = p_max[6:9])
      mapply(FUN = plot_par, par_name = mu_names[1:2], xmax = c(1, 1))

      ## 3rd row
      plot_axis()
      mapply(FUN = plot_par, par_name = mu_names[3:4], xmax = c(2, 2))
      plot_par("phi_pillar2_cases", xmax = 1)
      mapply(FUN = plot_par, par_name = alpha_names[1:3], xmax = 1)
    }

    ## 4th row
    plot_axis()
    mapply(FUN = plot_par, par_name = ch_names[1:3], xmax = c(1, 2e-5, 2e-5))
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
