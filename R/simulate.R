## 1. Load the data

## read in the combined fits
## - drop unwanted regions
## - upgrade state, perhaps
## - convert to multistrain
## - largely the old prepare function
spim_simulate_prepare <- function(combined, n_par,
                                  regions = NULL, inflate_strain = FALSE,
                                  inflate_booster = FALSE) {
  if (is.null(regions)) {
    regions <- sircovid::regions("all")
  }

  combined <- simulate_prepare_drop_regions(combined, regions)
  combined <- simulate_prepare_upgrade(combined)

  info <- combined$info

  ## Take a random sample of our parameters without replacement.
  n_regions <- length(regions)
  n_par_combined <- nrow(combined$pars[[1]])
  if (n_par > n_par_combined) {
    message(sprintf(
      "Reducing n_par from %d to %d as too few available in combined",
      n_par, n_par_combined))
    n_par <- n_par_combined
  }
  i <- sort(sample(n_par_combined, n_par, replace = FALSE))
  pars_mcmc <- lapply(combined$pars, function(x) x[i, , drop = FALSE])
  state <- lapply(combined$state, function(x) x[, i, drop = FALSE])

  ## TODO: this is quite slow
  message("Creating odin parameters from sampled parameters")
  pars <- Map(function(pars_region, transform)
    apply(pars_region, 1, transform),
    pars_mcmc, combined$transform)

  if (inflate_strain) {
    tmp <- simulate_prepare_inflate_strain(pars, state, info)
    pars <- tmp$pars
    state <- tmp$state
    info <- tmp$info
  }

  if (inflate_booster) {
    tmp <- simulate_prepare_inflate_vacc_classes(pars, state, info)
    pars <- tmp$pars
    state <- tmp$state
    info <- tmp$info
  }

  ## For the final object we will use a list-matrix of parameters and
  ## a 3d array of state as these will feed more easily into dust.
  pars <- unlist(unname(pars), FALSE)
  dim(pars) <- c(n_par, n_regions)
  colnames(pars) <- regions

  state <- array(unlist(unname(state)),
                 c(nrow(state[[1]]), n_par, n_regions))

  ## Our final object that we will use in the simulations
  ret <- combined[c("step", "date", "dt", "steps_per_day")]
  ret$pars <- pars
  ret$state <- state
  ret$info <- info
  ret
}


spim_simulate_one <- function(args, dat) {

}


spim_simulate_expand_grid <- function(...) {

}


spim_simulate_local <- function(dat, grid) {
  lapply(grid, spim_siumulate_one, dat)
}


spim_simulate_rrq <- function(dat, grid, rrq) {
  rrq$lapply(grid, spim_siumulate_one, dat)
}



### prep
simulate_prepare_upgrade <- function(combined) {
  ours <- packageVersion("sircovid")
  for (i in seq_along(combined$info)) {
    if (combined$info[[i]]$version != ours) {
      message(sprintf("Upgrading state for '%s' (%s => %s)",
                      names(combined$state)[[i]],
                      combined$info[[i]]$version, ours))
      ## NOTE: this won't work well if we have to add new parameters
      ## because we then break the transform function.
      p <- combined$transform[[i]](combined$pars[[i]][1, ])
      info_new <- sircovid::carehomes$new(p, 0, 1)$info()
      cmp <- combined$state[[i]]
      combined$state[[i]] <- sircovid::upgrade_state(
        combined$state[[i]],
        combined$info[[i]]$info,
        info_new)
      combined$info[[i]]$info <- info_new
    }
  }

  combined
}


simulate_prepare_drop_regions <- function(combined, regions) {
  nms <- setdiff(names(combined),
                 c("date", "step", "dt", "steps_per_day", "simulate"))
  for (i in nms) {
    msg <- setdiff(regions, names(combined[[i]]))
    if (length(msg) > 0) {
      stop(sprintf("Missing regions from %s: %s", i,
                   paste(squote(msg), collapse = ", ")))
    }
    combined[[i]] <- combined[[i]][regions]
  }
  combined
}


simulate_prepare_inflate_strain <- function(pars, state, info) {
  ## First inflate the parameters, because we need these in order to
  ## inflate the state.

  ## Easy updates:
  update <- list(cross_immunity = c(1, 1),
                 n_strains = 4,
                 strain_transmission = rep(1, 4))

  inflate_pars <- function(p_i) {
    p_i[names(update)] <- update
    for (rel in grep("^rel_", names(p_i), value = TRUE)) {
      rel_old <- p_i[[rel]]
      ## rel_gamma_X
      if (is.null(dim(rel_old))) {
        p_i[[rel]] <- rep(rel_old, 4)
      } else {
        rel_new <- array(0, c(nrow(rel_old), 4, nlayers(rel_old)))
        rel_new[, , ] <- rel_old[, 1, , drop = FALSE]
        p_i[[rel]] <- rel_new
      }
    }
    p_i
  }

  pars_new <- lapply(pars, lapply, inflate_pars)

  ## Then update the states given that:
  info_old <- info[[1]]$info
  info_new <- sircovid::carehomes$new(pars_new[[1]][[1]], 0, 1)$info()
  state_new <- lapply(state, sircovid::inflate_state_strains,
                      info_old, info_new)

  info_new <- lapply(info, function(x) {
    x$info <- info_new
    x
  })

  list(pars = pars_new, state = state_new, info = info_new)
}


## Sometimes if a run is done with no booster we still want
## to use it with boosters; this will update the parameters and the
## state in order to allow this. This takes the (almost) final output
## of simulate_prepare_
simulate_prepare_inflate_vacc_classes <- function(pars, state, info) {
  ## First inflate the parameters, because we need these in order to
  ## inflate the state.

  ## How many strata are we adding
  n_new_vacc_classes <- 1L # hardcoded for now

  ## of the new strata, are any corresponding to a dose (with a schedule)
  ## and where are they in the additional strata.
  ## hardcoded for now # first of the additional vaccine classes
  idx_new_dose <- 1L
  n_new_doses <- length(idx_new_dose)

  old_n_vacc_classes <- pars[[1]][[1]]$n_vacc_classes
  new_n_vacc_classes <- old_n_vacc_classes + n_new_vacc_classes
  old_idx_dose <- pars[[1]][[1]]$index_dose
  idx_booster <- old_idx_dose[length(old_idx_dose)] + idx_new_dose
  new_idx_dose <- c(old_idx_dose, idx_booster)
  new_n_vacc_classes <- pars[[1]][[1]]$n_vacc_classes + n_new_vacc_classes

  ## Easy updates: ## TOO: check if need anything else here
  update <- list(n_vacc_classes = new_n_vacc_classes,
                 n_doses = pars[[1]][[1]]$n_doses + n_new_doses,
                 vaccine_index_booster = idx_booster,
                 index_dose = new_idx_dose,
                 index_dose_inverse =
                   sircovid:::create_index_dose_inverse(new_n_vacc_classes,
                                                       new_idx_dose))

  inflate_pars <- function(p_i) {
    p_i[names(update)] <- update
    for (rel in grep("^rel_", names(p_i), value = TRUE)) {
      rel_old <- p_i[[rel]]
      if (!is.null(dim(rel_old))) {
        rel_new <- array(0, c(nrow(rel_old), ncol(rel_old),
                              nlayers(rel_old) + n_new_vacc_classes))
        rel_new[, , seq_len(nlayers(rel_old))] <- rel_old
        p_i[[rel]] <- rel_new
      }
    }
    rel <- "vaccine_progression_rate_base"
    rel_old <- p_i[[rel]]
    rel_new <- matrix(0, nrow(rel_old), ncol(rel_old) + n_new_vacc_classes)
    rel_new[, seq_len(ncol(rel_old))] <- rel_old
    p_i[[rel]] <- rel_new
    p_i
  }

  pars_new <- lapply(pars, lapply, inflate_pars)

  ## Then update the states given that:
  info_old <- info[[1]]$info
  info_new <- sircovid::carehomes$new(pars_new[[1]][[1]], 0, 1)$info()
  state_new <- lapply(state, sircovid::inflate_state_vacc_classes,
                      info_old, info_new)

  info_new <- lapply(info, function(x) {
    x$info <- info_new
    x
  })

  list(pars = pars_new, state = state_new, info = info_new)
}
