##' Load restart state from restart date.  This takes output from a
##' combined task that has been processed with
##' [spimalot::spim_fit_process] to return the model state at a given
##' date.
##'
##' @title Load restart data
##'
##' @param restart The loaded restart date (processed with
##'   [spimalot::spim_fit_process] in a previous task). We require that
##'   `restart` has a `state` element which is list with elements
##'   `time` (a sircovid date) and `state` (a 3d array of model state
##'   by particle by `time`).
##'
##' @param date The date that the restart should come from.  Will be
##'   passed through [sircovid::as_sircovid_date] so can be either an
##'   integer sircovid date, an R Date or a string representing an ISO
##'   date.
##'
##' @param round Logical parameter, indicates whether or not to round
##'   the main model states. Can be used to ensure that a restart initial
##'   state from a deterministic parent fit can be used in a stochastic
##'   restart fit. Default is FALSE
##'
##' @return A list with the same names as `restart`, but with
##'   `restart$state` filtered down to a single date's data and
##'   `restart$state$restart_date` replaced by the sircovid date
##'   corresponding to the restart.
##'
##' @export
spim_restart_initial_state <- function(restart, date, round = FALSE) {
  ## TODO: this could actually go into mcstate.
  date <- sircovid::as_sircovid_date(date)
  i <- match(date, restart$state$time)
  if (is.na(i)) {
    pos <- as.character(sircovid::sircovid_date_as_date(restart$state$time))
    stop(sprintf("Can't restart at date '%s', must be one of %s",
                 date, paste(pos, collapse = ", ")))
  }
  ret <- restart$state$state[, , i, drop = TRUE]

  if (round) {
    info <- restart$info$info[[length(restart$info$info)]]

    states_to_round <-
      c("D_hosp", "D_non_hosp", "D", "S", "vaccine_missed_doses",
        "T_sero_neg_1", "T_sero_neg_2", "R", "T_PCR_neg",
        "E", "I_A", "I_P", "I_C_1", "I_C_2", "G_D", "ICU_pre_unconf",
        "ICU_pre_conf", "H_R_unconf", "H_R_conf", "H_D_unconf", "H_D_conf",
        "ICU_W_R_unconf", "ICU_W_R_conf", "ICU_W_D_unconf", "ICU_W_D_conf",
        "ICU_D_unconf", "ICU_D_conf", "W_R_unconf", "W_R_conf", "W_D_unconf",
        "W_D_conf", "T_sero_pre_1", "T_sero_pos_1", "T_sero_pre_2",
        "T_sero_pos_2", "T_PCR_pre", "T_PCR_pos")

    idx_states_to_round <- c(unlist(info$index[states_to_round]))

    random_round <- function(x) {
      floor(x) + rbinom(prod(dim(x)), 1, x - floor(x))
    }

    ret[idx_states_to_round, ] <- random_round(ret[idx_states_to_round, ])
  }

  ret

}


##' Create restart parameters
##' @title Create restart parameters
##'
##' @param pars The full set of parameters.  This is a list with the
##'   baseline (`base`) along with `prior`, `info` and `transform`.
##'
##' @param pars_parent The parameters from the parent fit; this must
##'   be a `spim_pars_mcmc` object
##'
##' @param restart_date The date that restart happens, will be pased
##'   through `spimalot::as_sircovid_date`
##'
##' @return An object with the same structure as `pars` but with the
##'   `prior` element updated and an `mcmc` element added.
##'
##' @export
spim_restart_pars <- function(pars, pars_parent, restart_date) {
  assert_is(pars_parent, "spim_pars_pmcmc")
  assert_has_names(pars, c("base", "prior", "info", "transform"))

  ## These are parameters that we no longer will try and fit, and will
  ## eject from the full parameters object
  beta_date <- pars$base$beta_date
  restart_date <- sircovid::as_sircovid_date(restart_date)
  i <- max(which(beta_date <= restart_date))
  fixed <- c(sprintf("beta%d", seq_len(i - 1)), "start_date")

  ## These are the ones we should use the prior in the pars object,
  ## along with anything new. Note the -1 is because the restart data
  ## contains region (we might remove that)
  beta_restart <- sprintf("beta%d", seq(i, length(beta_date)))
  priors_propagate <- setdiff(pars_parent$info$name, beta_restart)
  pars$prior[match(priors_propagate, pars$prior$name), ] <-
    pars_parent$prior[match(priors_propagate, pars_parent$prior$name), -1]

  pars_full <- spim_pars_mcmc_single(pars$info, pars$prior, pars$proposal,
                                     pars$transform)

  pars$mcmc <- pars_full$fix(pars_full$initial()[fixed])
  pars
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
spim_restart_join_parent <- function(fit, parent) {
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
  ## NOTE: Previously we've saved ifr here too, but that needs
  ## reworking for the multistage fits still

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
