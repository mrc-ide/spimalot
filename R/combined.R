##' Load combined fits from a directory of regional fits
##'
##' @title Load combined fits
##'
##' @param path Directory name. Within here we expect to see
##'   `<region>/fits.rds` for each region in
##'   `sircovid::regions("all")`
##'
##' @param regions Region type passed to [sircovid::regions] (default
##'   is `all`, otherwise try `england`)
##'
##' @param get_severity Logical, indicating whether to extract severity (e.g.
##'   IFH, IHR, HFR, etc.) trajectories. Must default to `FALSE`, as these
##'   are routine (e.g. MTPs) `sircovid` outputs.
##'
##' @param get_onward Logical, indicating whether to extract projection
##'   trajectories for downstream (e.g. simulations) use. Default is `TRUE`,
##'   as these are routinely outputted in the RTM pipeline.
##'
##' @return A combined fit object
##' @export
spim_combined_load <- function(path, regions = "all", get_severity = FALSE,
                               get_onward = TRUE) {
  regions <- sircovid::regions(regions)

  files <- file.path(path, regions, "fit.rds")
  msg <- !file.exists(files)
  if (any(msg)) {
    msg <- sprintf("  - %s", file.path(regions[msg], "fit.rds"))
    stop(sprintf("Missing files at '%s':\n%s",
                 path, paste(msg, collapse = "\n")),
         call. = FALSE)
  }

  read_rds <- function(filename, region) {
    message(sprintf("Reading fit for %s", region))
    readRDS(filename)
  }
  dat <- Map(read_rds, files, regions)
  names(dat) <- regions

  info <- lapply(dat, function(x)
    x$samples$info[c("date", "multistrain", "model_type", "beta_date",
                     "restart_date")])
  if (length(unique(info)) != 1L) {
    stop("Incompatible regional fits")
  }

  ret <- list_transpose(dat)
  ret$info <- info[[1]]
  ret$info$date <- as.Date(ret$info$date)
  ## TODO: parent fits output restart_date as string and restart fits as
  ## sircovid date
  if (is.numeric(ret$info$restart_date)) {
    ret$info$restart_date <-
      sircovid::sircovid_date_as_date(ret$info$restart_date)
  } else {
    ret$info$restart_date <- as.Date(ret$info$restart_date)
  }

  message("Reordering trajectories")
  ## reorder by increasing cumulative incidence:
  rank_cum_inc <- lapply(ret$samples, sircovid::get_sample_rank)
  ret$samples <- Map(sircovid::reorder_sample, ret$samples, rank_cum_inc)
  if (!is.null(ret$rt[[1]])) {
    ret$rt <- Map(reorder_variant_rt, ret$rt, rank_cum_inc, weighted = TRUE)
  }

  message("Aggregating England/UK")
  ## Aggregate some of these to get england/uk entries
  ## Note that we do not store the aggregated outputs in ret yet to avoid
  ## including in the onward object below. The exception is the Rt values,
  ## as aggregated Rt values are used in onwards simulations
  agg_samples <- combined_aggregate_samples(ret$samples)
  agg_data <- combined_aggregate_data(ret$data)
  if (!is.null(ret$rt[[1]])) {
    ret$rt <- combined_aggregate_variant_rt(ret$rt, agg_samples, TRUE)
  }
  ret$model_demography <- combined_aggregate_demography(ret$model_demography)


  ## Check if severity calculations were outputed and, if so, aggregate
  if (get_severity) {
    message("Aggregating severity outputs")
    ret$severity <- Map(reorder_rt_ifr, ret$severity, rank_cum_inc)
    ret$severity <- combined_aggregate_severity(ret$severity, agg_samples)

    message("Aggregating prop_protected")
    agg_samples <- combined_aggregate_prop_protected(agg_samples)

    message("Aggregating and summarising intrinsic severity")
    ret$intrinsic_severity_raw <-
      Map(reorder_intrinsic_severity, ret$intrinsic_severity, rank_cum_inc)
    ret$intrinsic_severity_raw <-
      combined_aggregate_intrinsic_severity(ret$intrinsic_severity_raw,
                                            agg_samples)
    ret$intrinsic_severity <-
      summarise_intrinsic_severity(ret$intrinsic_severity_raw)
  }

  ## We don't need projections for the severity paper, these will be NULL to
  ## save memory as we already process very heavy outputs!
  ## We might want these for other analysis; output if get_onward == TRUE
  if (get_onward) {
    message("Creating data for onward use")
    ret$onward <- spim_combined_onward(ret)
  } else {
    ret$onward <- NULL
  }

  ## There are 3 elements in the parameter list that we need to join
  ## together; info, prior and proposal, anything else we will leave
  ## as a nested list (this includes the baseline parameter set which
  ## will come through as 'base')
  ret$parameters <- list_transpose(ret$parameters)
  pars_combine <- c("info", "prior", "proposal")
  ret$parameters[pars_combine] <-
    lapply(ret$parameters[pars_combine], dplyr::bind_rows)

  ## Now the onward object has been created, we can safely store the
  ## other aggregated outputs in ret
  ret$samples <- agg_samples
  ret$data <- agg_data

  ret
}


##' Load multiregion fits as a combined fit object
##'
##' @title Load multiregion fits as a combined fit object
##'
##' @param path Directory name. Within here we expect to see
##'   `fits.rds`
##'
##' @param get_severity Logical, indicating whether to extract severity (e.g.
##'   IFH, IHR, HFR, etc.) trajectories. Must default to `FALSE`, as these
##'   are routine (e.g. MTPs) `sircovid` outputs.
##'
##'
##' @return A combined fit object
##' @export
spim_combined_load_multiregion <- function(path, get_severity = FALSE) {

  message("Reading multiregion fit")
  filename <- file.path(path, "fit.rds")
  ret <- readRDS(filename)

  regions <- ret$samples$info$region

  ret$info <- ret$samples$info[c("date", "multistrain", "model_type",
                                 "beta_date", "restart_date")]
  ret$info$date <- as.Date(ret$info$date)

  region_samples <- function(r) {
    samples <- ret$samples
    samples$pars <- samples$pars_full[, , r]
    samples$probabilities <- samples$probabilities_full[, , r]
    samples$state <- samples$state[, r, ]
    samples$trajectories$state <- samples$trajectories$state[, r, , ]
    samples$predict$transform <- samples$predict$transform[[r]]
    samples
  }

  ret$samples <- lapply(seq_along(regions), region_samples)
  names(ret$samples) <- regions

  region_data <- function(region) {
    data <- ret$data
    data$fitted <- data$fitted[data$fitted$region == region, ]
    data$full <- data$full[data$full$region == region, ]
    data
  }

  ret$data <- lapply(regions, region_data)
  names(ret$data) <- regions

  message("Aggregating England/UK")
  ## Aggregate some of these to get england/uk entries
  ## Note that we do not store the aggregated outputs in ret yet to avoid
  ## including in the onward object below. The exception is the Rt values,
  ## as aggregated Rt values are used in onwards simulations
  agg_samples <- combined_aggregate_samples(ret$samples)
  agg_data <- combined_aggregate_data(ret$data)
  if (!is.null(ret$rt)) {
    ret$rt <- combined_aggregate_variant_rt(ret$rt, agg_samples, TRUE)
  }

  ## Check if severity calculations were outputed and, if so, aggregate
  if (get_severity) {
    message("Aggregating severity outputs")
    ret$severity <- combined_aggregate_severity(ret$severity, agg_samples)

    message("Aggregating prop_protected")
    agg_samples <- combined_aggregate_prop_protected(agg_samples)

    message("Aggregating and summarising intrinsic severity")
    ret$intrinsic_severity_raw <-
      combined_aggregate_intrinsic_severity(ret$intrinsic_severity,
                                            agg_samples)
    ret$intrinsic_severity <-
      summarise_intrinsic_severity(ret$intrinsic_severity_raw)
  }

  ## Now the onward object has been created, we can safely store the
  ## other aggregated outputs in ret
  ret$samples <- agg_samples
  ret$data <- agg_data

  ret
}


spim_combined_onward <- function(dat) {
  date <- dat$info$date
  steps_per_day <- dat$samples[[1]]$info$data$steps_per_day
  list(date = date,
       step = sircovid::sircovid_date(date) * steps_per_day,
       steps_per_day = steps_per_day,
       dt = 1 / steps_per_day,
       pars = lapply(dat$samples, "[[", "pars"),
       base = lapply(dat$parameters, function(x) x$base),
       state = lapply(dat$samples, "[[", "state"),
       data = lapply(dat$samples, function(x) x$predict$filter$data),
       transform = lapply(dat$samples, function(x) x$predict$transform),
       info = lapply(dat$samples, "[[", "info"),
       vaccine = lapply(dat$samples, "[[", "vaccine"),
       simulate = spim_combined_onward_simulate(dat))
}


spim_combined_onward_simulate <- function(dat) {
  simulate <- list_transpose(dat$simulate)

  dates <- dat$simulate[[1]]$date
  idx_dates <- dat$samples[[1]]$trajectories$date %in% dates
  state_names <- unique(c(rownames(simulate$state[[1]]), "deaths_hosp",
                          "sero_pos_1", "sero_pos_2"))
  state_names <- intersect(state_names,
                           rownames(dat$samples[[1]]$trajectories$state))

  state <- lapply(dat$samples, function(x)
    x$trajectories$state[state_names, , idx_dates])
  state <- aperm(abind_quiet(state, along = 4), c(1, 2, 4, 3))

  state_by_age <- lapply(list_transpose(simulate$state_by_age),
                        abind_quiet, along = 3)

  ret <- list(date = dates,
              state = state,
              state_by_age = state_by_age)

  ## This is not terrible:
  rt <- list_transpose(dat$rt)[c("Rt_general", "eff_Rt_general")]
    ## Rt_general and eff_Rt_general will have dimensions:
  ## [n particles x n regions x n dates]
  rt_combined <- lapply(rt, function(x)
    aperm(abind_quiet(x, along = 3), c(2, 3, 1))[, , idx_dates])

  idx_variant_dates <- dat$variant_rt[[1]]$date[, 1] %in% dates

  variant_rt <-
    list_transpose(dat$variant_rt)[c("Rt_general", "eff_Rt_general")]
  variant_rt_combined <- lapply(variant_rt, function(x)
    aperm(abind_quiet(x, along = 4), c(3, 4, 1, 2))[, , idx_variant_dates, ])

  rt_combined <- Map(function(rt, weighted_rt) {
    x <- abind_quiet(rt, weighted_rt, along = 4)
    dimnames(x)[[4]] <- c("strain_1", "strain_2", "both")
    x}, variant_rt_combined, rt_combined)

  ret <- c(ret, rt_combined)

  idx_dates_mv_rt <- dat$variant_rt[[1]]$date[, 1] %in% dates

  ## multivariant_Rt_general and multivariant_eff_Rt_general will have
  ## dimensions: [n particles x n regions x n variants x n dates]
  ## TODO: we are putting "multivariant" in the name here for clarity,
  ## so we may want to rename variant_rt wherever it appears to
  ## multivariant_rt
  mv_rt <-
    list_transpose(dat$variant_rt)[c("Rt_general", "eff_Rt_general")]
  mv_rt_combined <- lapply(mv_rt, function(x)
    aperm(abind_quiet(x, along = 4), c(3, 4, 2, 1))[, , , idx_dates_mv_rt])
  names(mv_rt_combined) <- paste0("multivariant_", names(mv_rt_combined))

  ret <- c(ret, mv_rt_combined)

  ret
}


combined_aggregate_samples <- function(samples) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(samples))) {
    samples$england$trajectories <-
      sircovid::combine_trajectories(samples[england], rank = FALSE)
  }

  if (all(nations %in% names(samples))) {
    samples$uk$trajectories <-
      sircovid::combine_trajectories(samples[nations], rank = FALSE)
  }

  samples
}


combined_aggregate_data <- function(data) {
  aggregate1 <- function(region_names, type = "full") {
    d <- lapply(data[region_names], "[[", type)
    x <- d[[1]]
    date_names <- grep("date", names(x), value = TRUE)
    data_names <- setdiff(names(x), date_names)
    dim_t <- min(sapply(d, nrow)) # TODO
    d <- lapply(d, function(x) x[seq_len(dim_t), data_names])
    ret <- Reduce("+", d)
    death_nms <- c("deaths_hosp", "deaths_carehomes", "deaths_comm")
    ret$deaths <- pmax(ret$deaths,
                       rowSums(ret[, death_nms], na.rm = TRUE),
                       na.rm = TRUE)
    ret[-seq_len(nrow(ret) - 4), death_nms] <- NA
    data_frame(ret, x[, date_names])
  }

  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(data))) {
    data$england$full <- aggregate1(england, "full")
    data$england$fitted <- aggregate1(england, "fitted")
    data$england$fitted$deaths <- NA
  }

  if (all(nations %in% names(data))) {
    data$uk$full <- aggregate1(nations, "full")
    data$uk$fitted <- aggregate1(nations, "fitted")
    data$uk$fitted$deaths <- NA
  }

  data
}


combined_aggregate_severity <- function(severity, samples) {

  # We first need to sum hospitalisations/infections incidence by strain across
  # metastrains
  for (i in names(samples)) {
    hospitalisations <- paste0("hospitalisations_inc_strain_", seq_len(4))
    samples[[i]]$trajectories$state[hospitalisations[1], , ] <-
      apply(samples[[i]]$trajectories$state[hospitalisations[c(1, 4)], , ],
            c(2, 3), sum, na.rm = TRUE)
    samples[[i]]$trajectories$state[hospitalisations[2], , ] <-
      apply(samples[[i]]$trajectories$state[hospitalisations[c(2, 3)], , ],
            c(2, 3), sum, na.rm = TRUE)

    infections <- paste0("infections_inc_strain_", seq_len(4))
    samples[[i]]$trajectories$state[infections[1], , ] <-
      apply(samples[[i]]$trajectories$state[infections[c(1, 4)], , ],
            c(2, 3), sum, na.rm = TRUE)
    samples[[i]]$trajectories$state[infections[2], , ] <-
      apply(samples[[i]]$trajectories$state[infections[c(2, 3)], , ],
            c(2, 3), sum, na.rm = TRUE)
  }

  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(severity))) {
    severity$england <- combined_severity(severity[england], samples[england])
  }

  if (all(nations %in% names(severity))) {
    severity$uk <- combined_severity(severity[nations], samples[nations])
  }

  severity
}


combined_aggregate_intrinsic_severity <- function(intrinsic_severity, samples) {

  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(intrinsic_severity))) {
    intrinsic_severity$england <-
      combined_intrinsic_severity(intrinsic_severity[england], samples[england])
  }

  if (all(nations %in% names(intrinsic_severity))) {
    intrinsic_severity$uk <-
      combined_intrinsic_severity(intrinsic_severity[nations], samples[nations])
  }

  intrinsic_severity
}


combined_aggregate_prop_protected <- function(samples) {

  for (r in names(samples)) {
    traj <- samples[[r]]$trajectories$state

    recovered_names <- c("recovered_1", "recovered_2", "recovered_historic")
    SR <- apply(traj[c("susceptible", recovered_names), , ], c(2, 3), sum)

    protected_names <- c("protected_S_vaccinated_weighted",
                         "protected_R_unvaccinated_weighted",
                         "protected_R_vaccinated_weighted")
    prop_protected <- vapply(protected_names,
                             function(i) {
                               traj[i, , ] / SR},
                             array(0, dim(SR)))
    prop_protected <- aperm(prop_protected, c(3, 1, 2))
    rownames(prop_protected) <-
      c("prop_protected_S_vaccinated",
        "prop_protected_R_unvaccinated",
        "prop_protected_R_vaccinated")

    eff_S_names <- c("effective_susceptible_1", "effective_susceptible_2")
    prop_protected_strain <- vapply(eff_S_names,
                                    function(i) {
                                      (SR - traj[i, , ]) / SR},
                                    array(0, dim(SR)))
    prop_protected_strain <- aperm(prop_protected_strain, c(3, 1, 2))
    rownames(prop_protected_strain) <- c("prop_protected_1", "prop_protected_2")

    prop_protected <-
      abind1(prop_protected, prop_protected_strain)
    samples[[r]]$trajectories$state <-
      abind1(traj, prop_protected)

  }
  samples
}


combined_severity <- function(severity, samples) {

  ## Produce England/UK level severity estimates. We weight IFRs and IHRs by
  ## the relevant infection incidence, HFRs by the relevant hospitalisations
  ## incidence

  out <- combined_severity_1(severity, samples, c("ifr", "ihr"),
                             "infections_inc", time_date = TRUE)

  out <- c(out, combined_severity_1(severity, samples, c("hfr"),
                                    "hospitalisations_inc"))

  age <- seq(0, 80, by = 5)
  for (i in age) {
    out <- c(out, combined_severity_1(severity, samples,
                                      paste0(c("ifr", "ihr"), "_age_", i),
                                      paste0("infections_inc_age_", i)))

    out <-
      c(out, combined_severity_1(severity, samples, paste0("hfr_age_", i),
                                 paste0("hospitalisations_inc_age_", i)))
  }

  strain <- c(1, 2)
  for (i in strain) {
    out <- c(out, combined_severity_1(severity, samples,
                                      paste0(c("ifr", "ihr"), "_strain_", i),
                                      paste0("infections_inc_strain_", i)))

    out <-
      c(out, combined_severity_1(severity, samples, paste0("hfr_strain_", i),
                                 paste0("hospitalisations_inc_strain_", i)))
  }


  out
}

combined_severity_1 <- function(severity, samples, what, weight,
                                time_date = FALSE) {
  for (i in names(severity)) {
    severity[[i]] <- severity[[i]][c(what, "time", "date")]

    ## Where incidence is zero (and thus where the weight is zero), the severity
    ## metric will be NA. To avoid the weighted severity being NA, we set the
    ## corresponding severity values to 0. The weighted severity value should
    ## therefore only be NA where the incidence across all regions is 0
    zero_weight <-
      which(t(samples[[i]]$trajectories$state[weight, , ]) == 0)
    for (w in what) {
      severity[[i]][[w]][zero_weight] <- 0
    }
  }

  ret <- sircovid::combine_rt(severity, samples,
                              rank = FALSE, weight = weight)

  if (!time_date) {
    ret$time <- NULL
    ret$date <- NULL
  }
  ret
}

combined_intrinsic_severity <- function(intrinsic_severity, samples) {

  get_reg_pop <- function(r) {
    p <- samples[[r]]$predict$transform(samples[[r]]$pars[1, ])
    sum(p[[1]]$pars$population)
  }

  wts <- vapply(names(samples), get_reg_pop, numeric(1))

  out <- list(period = intrinsic_severity[[1]]$period,
              variant = intrinsic_severity[[1]]$variant)

  intrinsic_severity <- list_transpose(intrinsic_severity)

  weight1 <- function(what) {
    sev <- abind_quiet(intrinsic_severity[[what]], along = 4)

    apply(sev, c(1, 2, 3), weighted.mean, w = wts)
  }

  out$IFR <- weight1("IFR")
  out$IHR <- weight1("IHR")
  out$HFR <- weight1("HFR")

  out
}


combined_aggregate_rt <- function(rt, samples) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(rt))) {
    rt$england <- sircovid::combine_rt(rt[england], samples[england],
                                       rank = FALSE)
  }
  if (all(nations %in% names(rt))) {
    rt$uk <- sircovid::combine_rt(rt[nations], samples[nations],
                                  rank = FALSE)
  }
  rt
}


combined_aggregate_variant_rt <- function(variant_rt, samples,
                                          weighted = FALSE) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(variant_rt))) {
    variant_rt$england <- combine_variant_rt(variant_rt[england],
                                             samples[england],
                                             rank = FALSE,
                                             weighted = weighted)
  }
  if (all(nations %in% names(variant_rt))) {
    variant_rt$uk <- combine_variant_rt(variant_rt[nations], samples[nations],
                                        rank = FALSE, weighted = weighted)
  }
  variant_rt
}


combine_variant_rt <- function(variant_rt, samples, rank, weighted = FALSE) {

  ## the samples trajectories and the variant_rt do not necessarily have the
  ## same span of dates so we need to filter the trajectories
  dates <- variant_rt[[1]]$date
  idx_dates <- samples[[1]]$trajectories$date %in% dates
  for (r in names(samples)) {
    samples[[r]]$trajectories$state <-
      samples[[r]]$trajectories$state[, , idx_dates]
    samples[[r]]$trajectories$date <-
      samples[[r]]$trajectories$date[idx_dates]
  }

  what <- setdiff(names(variant_rt[[1]]), c("time", "date", "beta"))

  ## combine_rt in sircovid only works for one variant so we need to combine
  ## each variant separately
  combine_variant_rt_j <- function(j) {
    get_variant_j <- function(v) {
      for (w in what) {
        v[[w]] <- v[[w]][, j, ]
      }
      v
    }

    sircovid::combine_rt(lapply(variant_rt, get_variant_j),
                         samples, rank)

  }

  variant1 <- combine_variant_rt_j(1)
  variant2 <- combine_variant_rt_j(2)
  if (weighted) {
    variant_weighted <- combine_variant_rt_j(3)
  }

  ## finally we join the variants together
  ret <- variant1
  for (w in what) {
    if (weighted) {
      ret[[w]] <- aperm(abind_quiet(variant1[[w]], variant2[[w]],
                                    variant_weighted[[w]],
                                    along = 3), c(1, 3, 2))
      colnames(ret[[w]]) <- c("strain_1", "strain_2", "weighted")
    } else {
      ret[[w]] <- aperm(abind_quiet(variant1[[w]], variant2[[w]],
                                    along = 3), c(1, 3, 2))
    }

  }

  ret
}


combined_aggregate_demography <- function(demography) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(demography))) {
    demography$england <- combine_demography(demography[england])
  }
  if (all(nations %in% names(demography))) {
    demography$uk <- combine_demography(demography[nations])
  }
  demography
}


combine_demography <- function(demography) {
  demography <- combined_switch_levels(demography)

  agg_demography <- function(demo) {
    demo <- combined_switch_levels(demo)
    total <- Reduce(`+`, demo$total)
    mean_t <- Reduce(`+`, demo$mean_t)

    prop_total <- t(total) / colSums(total) * 100

    list(total = total,
         prop_total = prop_total,
         mean_prop_total = colMeans(prop_total),
         mean_t = mean_t,
         date = demo$date[[1]])
  }

  ret <- lapply(demography, agg_demography)
  ret
}


## Used to invert regions/variable in rt data
combined_switch_levels <- function(x) {
  nms <- names(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}

##' Get region and country population from a combined fits object
##' @title Get population from combined
##'
##' @param combined Combined fits object
##'
##' @param ignore_uk If `TRUE` population of UK is returned as NA
##'
##' @param by_age Logical, indicating if population should be returned by age
##'
##' @param group Vector of group indices (for the carehomes model on
##'   `1:19`)
##'
##' @return data.frame of region/country populations
##' @export
spim_population <- function(combined, ignore_uk = FALSE, by_age = TRUE,
                            group = 1:19) {
  pop <- vapply(names(combined$pars), function(r)
    as.integer(combined$transform[[r]](combined$pars[[r]][1, ])$N_tot[group]),
    group)
  pop_england <- rowSums(pop[, sircovid::regions("england")])
  if (!ignore_uk) {
    pop_uk <- rowSums(pop[, sircovid::regions("all")])
  } else {
    pop_uk <- rep(NA_real_, length(pop_england))
  }

  df <- data.frame(pop, england = pop_england, uk = pop_uk)
  if (!by_age) {
    df <- colSums(df)
  }

  df
}

##' Get proportion infected from country population from a combined fits object
##' @title Get proportion infected from combined
##' @param combined Combined fits object
##' @param population Population from [spim_population]
##' @param regions Regions for which to get proportion infected
##' @return data.frame of region/country proportion infected
##' @export
spim_prop_infected <- function(combined, population,
                              regions = sircovid::regions("england")) {
  n_epoch <- length(combined$info[[1]]$info)
  idx_cum_infections <-
    combined$info[[1]]$info[[n_epoch]]$index$cum_infections
  prop_infected <- sapply(regions, function(r) {
    mean_ci(combined$state[[r]][idx_cum_infections, ]) / sum(population[[r]])
  })
  t(prop_infected)
}

reorder_variant_rt <- function(x, rank, weighted = FALSE) {

  what <- setdiff(names(x), c("time", "date"))

  ## reorder_rt_ifr will only work for one variant so we will have to reorder
  ## each variant separately
  if (weighted) {
    dim_strains <- 3
  } else {
    dim_strains <- 2
  }
  for (j in seq_len(dim_strains)) {
    v <- x

    for (i in what[what != "beta"]) {
      v[[i]] <- v[[i]][, j, ]
    }

    v <- sircovid::reorder_rt_ifr(v, rank)

    for (i in what) {
      if (i == "beta") {
        x[[i]] <- v[[i]]
      } else {
        x[[i]][, j, ] <- v[[i]]
      }

    }
  }

  x
}

reorder_intrinsic_severity <- function(x, rank) {

  what <- setdiff(names(x), c("period", "variant"))

  dim_strains <- 2

  for (j in seq_len(dim_strains)) {
    v <- x[what]

    for (i in what) {
      v[[i]] <- v[[i]][, j, ]
    }

    class(v) <- "IFR_t_trajectories"
    v <- sircovid::reorder_rt_ifr(v, rank)

    for (i in what) {
      x[[i]][, j, ] <- v[[i]]
    }
  }

  x
}


summarise_intrinsic_severity <- function(intrinsic_severity) {

  regions <- names(intrinsic_severity)

  periods <- intrinsic_severity[[1]]$period
  variants <- intrinsic_severity[[1]]$variant

  what <- setdiff(names(intrinsic_severity[[1]]), c("period", "variant"))

  sev_vector <- function(x) {
    mean <- mean(x)
    lb <- quantile(x, 0.025)
    ub <- quantile(x, 0.975)

    out <- c(mean, lb, ub)
    names(out) <- c("mean", "lb", "ub")
    out
  }

  get_region <- function(r) {

    get_what <- function(w) {
      get_variant <- function(j) {
        y <- t(apply(intrinsic_severity[[r]][[w]][, j, ], 1, sev_vector))
        data.frame(period = periods, y) %>%
          pivot_longer(!period, names_to = "estimate")
      }

      v <- lapply(seq_len(length(variants)), get_variant)
      names(v) <- variants
      dplyr::bind_rows(v, .id = "name")
    }

    ret_what <- lapply(what, get_what)
    names(ret_what) <- what
    dplyr::bind_rows(ret_what, .id = "source")

  }

  ret <- lapply(regions, get_region)
  names(ret) <- regions
  ret <- dplyr::bind_rows(ret, .id = "region")
  ret$period <- factor(ret$period, levels = periods)

  ret
}

##' Add forecasts to combined fit object
##'
##' @title Add forecasts to combined fit object
##'
##' @param dat A combined fit object
##'
##' @param forecast_days Number of days of forecasts to add
##'
##' @return A combined fit object with forecasts added in
##' @export
spim_add_forecast <- function(dat, forecast_days) {

  message("Adding forecasts")
  dat$samples <- lapply(dat$samples[sircovid::regions("all")],
                        function(x) spim_add_forecast_region(x, forecast_days))
  dat$samples <- combined_aggregate_samples(dat$samples)

  dat
}


spim_add_forecast_region <- function(samples, forecast_days) {

  ## Run forecast
  time_predict <- seq(samples$predict$time,
                      length.out = forecast_days + 1L,
                      by = samples$predict$rate)
  forecast <- spim_pmcmc_predict(
    samples, time_predict,
    prepend_trajectories = FALSE)
  forecast$date <- forecast$time / forecast$rate

  ## Get forecast trajectories in a shape compatible with samples trajectories
  forecast_samples <- samples
  forecast_samples$trajectories <- forecast
  forecast_samples <- reduce_trajectories(forecast_samples, FALSE)

  ## Join trajectories from forecast and samples together
  ## We remove the first date from the forecast as this is the same as the last
  ## fitted date
  samples$trajectories$state <-
    abind_quiet(samples$trajectories$state,
                forecast_samples$trajectories$state[, , -1L],
                along = 3)
  samples$trajectories$time <- c(samples$trajectories$time,
                                 forecast_samples$trajectories$time[-1L])
  samples$trajectories$date <- c(samples$trajectories$date,
                                 forecast_samples$trajectories$date[-1L])
  samples$trajectories$predicted <-
    c(samples$trajectories$predicted,
      forecast_samples$trajectories$predicted[-1L])

  samples

}

## This is adapted from mcstate::pmcmc_predict
spim_pmcmc_predict <- function(object, time, prepend_trajectories = FALSE,
                               n_threads = NULL, seed = NULL) {

  if (is.null(object$predict)) {
    stop("mcmc was run with return_state = FALSE, can't predict")
  }
  if (length(time) < 2) {
    stop("At least two time steps required for predict")
  }
  if (time[[1]] != object$predict$time) {
    stop(sprintf("Expected time[1] to be %d", object$predict$time))
  }
  if (prepend_trajectories && is.null(object$trajectories)) {
    stop(paste("mcmc was run with return_trajectories = FALSE,",
               "can't prepend trajectories"))
  }

  state <- object$state
  info <- object$info$info[[length(object$info$info)]]
  index <- sircovid::lancelot_index(info)$state
  model <- object$predict$filter$model
  n_threads <- n_threads %||% object$predict$filter$n_threads

  pars <- lapply(seq_len(nrow(object$pars)), function(i)
    last(object$predict$transform(object$pars[i, ]))$pars)

  mod <- model$new(pars, time[[1]], NULL, n_threads = n_threads,
                   seed = seed, pars_multi = TRUE)

  mod$update_state(state = state)
  if (!is.null(index)) {
    mod$set_index(index)
  }
  y <- mod$simulate(time)

  res <-
    mcstate:::mcstate_trajectories_discrete(time, object$predict$rate, y, TRUE)

  if (prepend_trajectories) {
    res <-
      mcstate:::bind_mcstate_trajectories_discrete(object$trajectories, res)
  }

  res
}
