## NOTE: throughout this file the use of region is quite peciular,
## designed to work only(?!) with the multi-region fits, and if you
## give it for single region fits bad things will happen.

##' Plots for fits
##'
##' @title Plots for fit
##'
##' @param samples The samples
##'
##' @param region The region. Leave this `NULL` except for
##'   multi-region fits
##'
##' @return Nothing, called for side effect
##' @rdname spim_plot_fit
##' @export
spim_plot_fit_forecasts <- function(samples, region = NULL) {
  data_date <- samples$info$date
  ## NOTE: The first date will include the whole pre-data counts, so
  ## drop here.
  step <- samples$trajectories$step[-1]
  date <- sircovid::sircovid_date_as_date(step / samples$predict$rate)

  info <- list(
    list(index = "deaths_hosp",
         ylab = "Deaths",
         lag = TRUE),
    list(index = "icu",
         ylab = "Confirmed covid-19 patients in ICU",
         lag = FALSE),
    list(index = "general",
         ylab = "Confirmed covid-19 patients in general beds",
         lag = FALSE),
    list(index = "diagnoses",
         ylab = "Inpatients newly-diagnosed with covid-19",
         lag = TRUE),
    list(index = "admitted",
         ylab = "Patients admitted with covid-19",
         lag = TRUE))

  plot1 <- function(x, region = NULL) {
    if (is.null(region)) {
      y <- t(samples$trajectories$state[x$index, , -1])
    } else {
      which_lay <- match(region, colnames(samples$probabilities))
      which_row <- match(x$index, names(samples$predict$index))
      y <- t(samples$trajectories$state[which_row, , which_lay, -1])
    }
    if (x$lag) {
      y <- rbind(0, diff(y))
    }
    ## matplot can't do date axes, so we pre-plot with plot()
    plot(date, y[, 1], type = "n", ylab = x$ylab, ylim = range(y, na.rm = TRUE))
    matlines(date, y, lty = 1, col = "#cccccc80")
    quantiles <- t(apply(y, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE))
    matlines(date, quantiles, col = "black", lty = 2)
  }

  par(mfrow = c(1, length(info)),
      mar = c(4, 4, 3, 1),
      oma = c(1, 1, 1, 1),
      mgp = c(2, 0.5, 0))
  for (x in info) {
    plot1(x, region)
  }
  if (is.null(region)) {
    text <- sprintf("Forecasts run on %s", data_date)
  } else {
    text <- sprintf("%s forecasts run on %s", region, data_date)
  }
  mtext(text, 3, outer = TRUE, line = -1)
}


##' @rdname spim_plot_fit
##' @export
spim_plot_fit_traces <- function(samples) {
  multiregion <- samples$info$multiregion
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    ## We assume this below
    n_chains <- length(unique(samples$chain))
    stopifnot(
      identical(samples$chain,
                rep(seq_len(n_chains),
                    each = length(samples$chain) / n_chains)))
  }
  cols <- rev(viridisLite::viridis(n_chains))

  i <- reorder_beta(colnames(samples$pars_full))
  if (multiregion) {
    pars <- samples$pars_full[, i, ]
  } else {
    pars <- samples$pars_full[, i]
  }
  nms <- colnames(pars)
  probs <- samples$probabilities_full

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  new_grid <- function(n, title) {
    par(mfrow = rep(ceiling(sqrt(n + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1 + as.integer(title), 1))
  }

  plot_traces1 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    ess <- coda::effectiveSize(coda::as.mcmc(traces))

    if (name == "log_likelihood") {
      main <- ""
    } else {
      main <- sprintf("ess = %s", round(sum(ess)))
    }
    matplot(traces, type = "l", lty = 1,
              xlab = "Iteration", bty = "n",
              ylab = name, col = cols,
              main = main,
              font.main = 1)
    rug(samples$iteration[samples$chain == 1], ticksize = 0.1)
  }

  if (multiregion) {
    nms_fixed <- samples$info$pars$fixed
    nms_varied <- samples$info$pars$varied
    region <- samples$info$region

    new_grid(length(nms_fixed), TRUE)
    for (nm in nms_fixed) {
      plot_traces1(pars[, nm, 1], nm)
    }
    plot_traces1(rowSums(probs[, "log_likelihood", ]), "log_likelihood")
    mtext("Fixed (shared across regions)", outer = TRUE)

    for (r in region) {
      new_grid(length(nms_varied), TRUE)
      for (nm in nms_varied) {
        plot_traces1(pars[, nm, r], nm)
      }
      plot_traces1(probs[, "log_likelihood", r], "log_likelihood")
      mtext(r, outer = TRUE)
    }
  } else {
    new_grid(length(nms), FALSE)
    for (nm in nms) {
      plot_traces1(samples$pars_full[, nm], nm)
    }
    plot_traces1(probs[, "log_likelihood"], "log_likelihood")
  }
}


##' @rdname spim_plot_fit
##' @export
spim_plot_fit_traces_separate <- function(samples) {
  multiregion <- samples$info$multiregion
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    ## We assume this below
    n_chains <- length(unique(samples$chain))
    stopifnot(
      identical(samples$chain,
                rep(seq_len(n_chains),
                    each = length(samples$chain) / n_chains)))
  }
  cols <- rev(viridisLite::viridis(n_chains))

  i <- reorder_beta(colnames(samples$pars_full))
  if (multiregion) {
    pars <- samples$pars_full[, i, ]
  } else {
    pars <- samples$pars_full[, i]
  }
  nms <- colnames(pars)
  probs <- samples$probabilities_full

  op <- par(no.readonly = TRUE)
  on.exit(par(op))


  plot_traces2 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    ess <- coda::effectiveSize(coda::as.mcmc(traces))

    if (name == "log_likelihood") {
      main <- ""
    } else {
      main <- sprintf("ess = %s", round(sum(ess)))
    }
    matplot(traces, type = "l", lty = 1,
            xlab = "Iteration", bty = "n",
            ylab = name, col = cols,
            main = main,
            font.main = 1)
    rug(samples$iteration[samples$chain == 1], ticksize = 0.1)
  }

  if (multiregion) {
    nms_fixed <- samples$info$pars$fixed
    nms_varied <- samples$info$pars$varied
    region <- samples$info$region

    for (nm in nms_fixed) {
      plot_traces2(pars[, nm, 1], nm)
    }
    plot_traces2(rowSums(probs[, "log_likelihood", ]), "log_likelihood")
    mtext("Fixed (shared across regions)", outer = TRUE)

    for (r in region) {
      for (nm in nms_varied) {
        plot_traces2(pars[, nm, r], nm)
      }
      plot_traces2(probs[, "log_likelihood", r], "log_likelihood")
      mtext(r, outer = TRUE)
    }
  } else {
    for (nm in nms) {
      plot_traces2(samples$pars_full[, nm], nm)
    }
    plot_traces2(probs[, "log_likelihood"], "log_likelihood")
  }
}


##' @rdname spim_plot_fit
##' @export
spim_plot_fit_posteriors <- function(samples, region = NULL) {
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    ## We assume this below
    n_chains <- length(unique(samples$chain))
    stopifnot(identical(
      samples$chain,
      rep(seq_len(n_chains), each = length(samples$chain) / n_chains)))
  }
  cols <- rev(viridisLite::viridis(n_chains))

  plot_posteriors1 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    if (diff(range(traces)) == 0) {
      ## We could pick up nicer things about the parameter range here,
      ## but this is only an issue when things have gone very wrong,
      ## or on a short run
      breaks <- seq(min(0, min(traces)), max(1, max(traces)), length.out = 20)
    } else {
      breaks <- seq(min(traces), max(traces), length.out = 20)
    }
    h <- apply(traces, 2, hist, plot = FALSE, breaks = breaks)

    y <- vapply(h, "[[", numeric(length(breaks) - 1), "density")
    y <- rbind(y, y[nrow(y), ])


    if (is.null(region)) {
      ess <- coda::effectiveSize(coda::as.mcmc(traces))
      main <- ifelse(name == "log_likelihood", "",
                   paste("ess =", round(sum(ess))))
    } else {
      ess <- coda::effectiveSize(coda::as.mcmc(traces))
      main <- ifelse(name == "log_likelihood", sprintf("Region %s", region),
                     sprintf("%s; ess = %d", region, round(sum(ess))))
    }

    matplot(breaks, y, type = "s", lty = 1,
              xlab = name, ylab = "density", col = cols, main = main,
              font.main = 1, bty = "n")
  }

  if (is.null(region)) {

    n_pars <- length(colnames(samples$pars))

    par(mfrow = rep(ceiling(sqrt(n_pars + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1, 1))

    for (nm in colnames(samples$pars)) {
      plot_posteriors1(samples$pars[, nm], nm)
    }
  } else {

    n_pars <- length(rownames(samples$pars))

    par(mfrow = rep(ceiling(sqrt(n_pars + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1, 1))

    for (nm in rownames(samples$pars)) {
      plot_posteriors1(samples$pars[nm, region, ], nm)
    }
  }

  legend("left", fill = cols, bty = "n",
         legend = paste("chain", seq_len(n_chains)))
  rhat <- gelman_diagnostic(samples, region)
  if (!is.null(rhat)) {
    mtext(side = 3, text = paste("Rhat =", round(rhat$mpsrf, 1)))
  }
}


gelman_diagnostic <- function(samples, region = NULL) {
  if (is.null(samples$chain)) {
    return(NULL)
  }

  if (is.null(region)) {
    chains <- lapply(unname(split(data.frame(samples$pars), samples$chain)),
                    coda::as.mcmc)
  } else {
    chains <- lapply(unname(split(data.frame(t(samples$pars[, region, ])),
                                  samples$chain)),
                     coda::as.mcmc)
  }

  tryCatch(
    coda::gelman.diag(chains),
    error = function(e) NULL)
}


reorder_beta <- function(nms) {
  i <- grep("^beta[0-9]+$", nms)
  j <- order(as.integer(sub("^beta", "", nms[i])))
  k <- seq_along(nms)
  k[i] <- i[j]
  k
}
