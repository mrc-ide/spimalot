##' Load restart state from restart date.  This takes output from a
##' combined task that has been processed with
##' [spimalot::spim_fit_process] to return the model state at a given
##' date.
##'
##' @title Load restart data
##'
##' @param restart The loaded restart date (processed with
##'   [spimalot::spim_fit_pocess] in a previous task). We require that
##'   `restart` has a `state` element which is list with elements
##'   `time` (a sircovid date) and `state` (a 3d array of model state
##'   by particle by `time`).
##'
##' @param date The date that the restart should come from.  Will be
##'   passed through [sircovid::as_sircovid_date] so can be either an
##'   integer sircovid date, an R Date or a string representing an ISO
##'   date.
##'
##' @return A list with the same names as `restart`, but with
##'   `restart$state` filtered down to a single date's data and
##'   `restart$state$restart_date` replaced by the sircovid date
##'   corresponding to the restart.
##'
##' @export
spim_restart_initial_state <- function(restart, date) {
  ## TODO: this could actually go into mcstate.
  date <- sircovid::as_sircovid_date(date)
  i <- match(date, restart$state$time)
  if (is.na(i)) {
    pos <- as.character(sircovid::sircovid_date_as_date(restart$state$time))
    stop(sprintf("Can't restart at date '%s', must be one of %s",
                 date, paste(pos, collapse = ", ")))
  }
  restart$state$state[, , i, drop = TRUE]
}


##' Join parent and restart fits together
##'
##' @title Join parent and restart fits
##'
##' @param fit The new fit. A list with elements `samples`, `pmcmc`,
##'   etc.
##'
##' @param parent The parent information from the restart data created
##'   by [spimalot::spim_fit_process]; this contains elements
##'   `trajectories`, `rt` and `data` and we use these to
##'   create consistent trajectories that cover both the parent and
##'   restart fits.
##'
##' @return A list with the same elements as `fit`, but with
##'   concatenated trajectories as appropriate.
##'
##' @export
spim_restart_join_parent <- function(fit, parent, data) {
  ## First, fix first step; see
  ## https://github.com/mrc-ide/mcstate/issues/55
  fit$samples$trajectories$step <- fit$samples$trajectories$step[-1L]
  fit$samples$trajectories$date <- fit$samples$trajectories$date[-1L]
  fit$samples$trajectories$predicted <- fit$samples$trajectories$predicted[-1L]
  fit$samples$trajectories$state <- fit$samples$trajectories$state[, , -1L]

  ## Then stitch together. Work out what we keep from the parent fit:
  ## This is the restart date. There are many other places in the
  ## object this is stored, but this one is ok at least?  We can also
  ## get it from the trajectories
  stopifnot(
    isTRUE(fit$samples$trajectories$date[[1]] == fit$data$fitted$date[[1]]))
  restart_date <- fit$data$fitted$date[[1]] - 1L
  i <- which(parent$trajectories$date <= restart_date)

  ## Filter trajectories:
  for (v in c("step", "date", "predicted")) {
    fit$samples$trajectories[[v]] <- c(
      parent$trajectories[[v]][i], fit$samples$trajectories[[v]])
  }
  fit$samples$trajectories$state <- mcstate::array_bind(
    parent$trajectories$state[, , i, drop = FALSE],
    fit$samples$trajectories$state)

  join_rt <- function(parent, new, i) {
    ret <- Map(function(a, b) rbind(a[i, , drop = FALSE], b[-1L, ]),
               parent, new)
    class(ret) <- class(new)
    ret
  }

  fit$rt <- join_rt(parent$rt, fit$rt, i)
  fit$ifr_t <- join_rt(parent$ifr_t, fit$ifr_t, i)

  ## Here we need things from the parent fit really
  i_data <- which(parent$data$full$date <= restart_date)
  fit$data$full <- rbind(parent$data$full[i_data, ], fit$data$full)
  fit$data$fitted <- rbind(parent$data$fitted[i_data, ], fit$data$fitted)

  ## Extra metadata for the fit:
  fit$samples$info$restart <- TRUE
  fit$samples$info$restart_date <- restart_date

  predicted <- fit$samples$trajectories$predicted
  if (any(predicted)) {
    stop("Predicted not supported at the moment")
  }

  fit$samples$info$time_index <- list(
    parent = i,
    restart = seq(max(i) + 1, length(fit$samples$trajectories$step)),
    predicted = integer(0))

  fit
}
