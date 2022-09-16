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
  ret$rt <- Map(sircovid::reorder_rt_ifr, ret$rt, rank_cum_inc)
  ret$variant_rt <- Map(reorder_variant_rt, ret$variant_rt, rank_cum_inc)

  message("Aggregating England/UK")
  ## Aggregate some of these to get england/uk entries
  ## Note that we do not store the aggregated outputs in ret yet to avoid
  ## including in the onward object below. The exception is the Rt values,
  ## as aggregated Rt values are used in onwards simulations
  agg_samples <- combined_aggregate_samples(ret$samples)
  agg_data <- combined_aggregate_data(ret$data)
  ret$rt <- combined_aggregate_rt(ret$rt, agg_samples)
  ret$variant_rt <- combined_aggregate_variant_rt(ret$variant_rt, agg_samples)
  ret$model_demography <- combined_aggregate_demography(ret$model_demography)


  ## Check if severity calculations were outputed and, if so, aggregate
  if (get_severity) {
    message("Aggregating severity outputs")
    ret$severity <- combined_aggregate_severity(ret$severity, agg_samples)
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
##' @return A combined fit object
##' @export
spim_combined_load_multiregion <- function(path) {

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
  ret$rt <- combined_aggregate_rt(ret$rt, agg_samples)
  ret$variant_rt <- combined_aggregate_variant_rt(ret$variant_rt, agg_samples)

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
    agg_positivity <-
      any(grepl("pillar2_positivity",
                rownames(samples$england$trajectories$state)))
    if (agg_positivity) {
      samples$england$trajectories <-
        combined_aggregate_positivity(samples$england$trajectories,
                                      samples[england])
    }
  }

  if (all(nations %in% names(samples))) {
    samples$uk$trajectories <-
      sircovid::combine_trajectories(samples[nations], rank = FALSE)
    agg_positivity <-
      any(grepl("pillar2_positivity",
                rownames(samples$uk$trajectories$state)))
    if (agg_positivity) {
      samples$uk$trajectories <-
        combined_aggregate_positivity(samples$uk$trajectories,
                                      samples[sircovid::regions("all")])
    }
  }

  samples
}

combined_aggregate_positivity <- function(agg_trajectories, samples) {

  p_NC_names <- c("p_NC_under15", "p_NC_15_24", "p_NC_25_49",
                  "p_NC_50_64", "p_NC_65_79", "p_NC_80_plus",
                  "p_NC_weekend_under15", "p_NC_weekend_15_24",
                  "p_NC_weekend_25_49", "p_NC_weekend_50_64",
                  "p_NC_weekend_65_79", "p_NC_weekend_80_plus")

  dim <- dim(samples[[1]]$trajectories$state)[c(2, 3)]

  pos_under15 <- array(0, dim = dim)
  pos_15_24 <- array(0, dim = dim)
  pos_25_49 <- array(0, dim = dim)
  pos_50_64 <- array(0, dim = dim)
  pos_65_79 <- array(0, dim = dim)
  pos_80_plus <- array(0, dim = dim)

  neg_under15 <- array(0, dim = dim)
  neg_15_24 <- array(0, dim = dim)
  neg_25_49 <- array(0, dim = dim)
  neg_50_64 <- array(0, dim = dim)
  neg_65_79 <- array(0, dim = dim)
  neg_80_plus <- array(0, dim = dim)

  for (r in names(samples)) {
    state <- samples[[r]]$trajectories$state
    date <- sircovid::sircovid_date_as_date(samples[[r]]$trajectories$date)

    pos_under15 <- pos_under15 + state[paste0("sympt_cases_under15_inc"), , ]
    pos_15_24 <- pos_15_24 + state[paste0("sympt_cases_15_24_inc"), , ]
    pos_25_49 <- pos_25_49 + state[paste0("sympt_cases_25_49_inc"), , ]
    pos_50_64 <- pos_50_64 + state[paste0("sympt_cases_50_64_inc"), , ]
    pos_65_79 <- pos_65_79 + state[paste0("sympt_cases_65_79_inc"), , ]
    pos_80_plus <- pos_80_plus + state[paste0("sympt_cases_80_plus_inc"), , ]

    pars_model <- lapply(seq_len(nrow(samples[[r]]$pars)), function(i)
      last(samples[[r]]$predict$transform(samples[[r]]$pars[i, ]))$pars)

    pars_base <- pars_model[[1]]

    pars <- t(vapply(pars_model, function(p) unlist(p[p_NC_names]),
                     numeric(length(p_NC_names))))

    calc_negs <- function(group) {
      neg <- pars_base[[paste0("N_tot_", group)]] -
        state[paste0("sympt_cases_", group, "_inc"), , ]

      neg[, grepl("^S", weekdays(date))] <-
        neg[, grepl("^S", weekdays(date))] *
        pars[, paste0("p_NC_weekend_", group)]
      neg[, !grepl("^S", weekdays(date))] <-
        neg[, !grepl("^S", weekdays(date))] *
        pars[, paste0("p_NC_", group)]

      neg
    }

    neg_under15 <- neg_under15 + calc_negs("under15")
    neg_15_24 <- neg_15_24 + calc_negs("15_24")
    neg_25_49 <- neg_25_49 + calc_negs("25_49")
    neg_50_64 <- neg_50_64 + calc_negs("50_64")
    neg_65_79 <- neg_65_79 + calc_negs("65_79")
    neg_80_plus <- neg_80_plus + calc_negs("80_plus")
  }

  pos_over25 <- pos_25_49 + pos_50_64 + pos_65_79 + pos_80_plus
  pos_all <- pos_under15 + pos_15_24 + pos_over25

  neg_over25 <- neg_25_49 + neg_50_64 + neg_65_79 + neg_80_plus
  neg_all <- neg_under15 + neg_15_24 + neg_over25

  pars_base <-
    last(samples[[1]]$predict$transform(samples[[1]]$pars[1, ]))$pars

  calc_pos <- function(pos, neg) {
    positivity1 <-
      (pos * pars_base$pillar2_sensitivity +
         neg * (1 - pars_base$pillar2_specificity)) / (pos + neg) * 100
    array(positivity1, c(1, dim(positivity1)))
  }

  positivity <- calc_pos(pos_all, neg_all)
  positivity <- abind1(positivity, calc_pos(pos_over25, neg_over25))
  positivity <- abind1(positivity, calc_pos(pos_under15, neg_under15))
  positivity <- abind1(positivity, calc_pos(pos_15_24, neg_15_24))
  positivity <- abind1(positivity, calc_pos(pos_25_49, neg_25_49))
  positivity <- abind1(positivity, calc_pos(pos_50_64, neg_50_64))
  positivity <- abind1(positivity, calc_pos(pos_65_79, neg_65_79))
  positivity <- abind1(positivity, calc_pos(pos_80_plus, neg_80_plus))

  row.names(positivity) <- paste0("pillar2_positivity",
                                  c("", "_over25", "_under15", "_15_24",
                                    "_25_49", "_50_64", "_65_79", "_80_plus"))

  agg_trajectories$state[rownames(positivity), , ] <- positivity
  agg_trajectories
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

  # We first need to sum admitted_inc and diagnoses_inc trajectories to get
  # all_admissions_inc, which we'll use as weight to aggregate hospital severity
  for (i in names(samples)) {
    admitted_inc <- samples[[i]]$trajectories$state["admitted_inc", , ,
                                                    drop = FALSE]
    diagnoses_inc <- samples[[i]]$trajectories$state["diagnoses_inc", , ,
                                                     drop = FALSE]
    all_admissions <- admitted_inc + diagnoses_inc
    rownames(all_admissions) <- "all_admissions_inc"
    samples[[i]]$trajectories$state <-
      abind_quiet(samples[[i]]$trajectories$state, all_admissions, along = 1)
  }

  # Some unfortunate double running of the same function outputting a similar
  # object but with different weighting
  admissions <- combined_aggregate_severity_1(severity, samples,
                                              weight = "all_admissions_inc")

  infections <- combined_aggregate_severity_1(severity, samples,
                                              weight = "infections_inc")

  # Severity trajectories weighted by admissions
  admission_weighted <- grep("^hfr", names(admissions[[1]]), value = TRUE)

  severity$england <- infections$england
  for (i in admission_weighted) {
    severity$england[[i]] <- admissions$england[[i]]
  }

  severity
}


combined_aggregate_severity_1 <- function(x, samples, weight) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(x))) {
    x$england <- sircovid::combine_rt(x[england], samples[england],
                                       rank = FALSE, weight = weight)
  }

  if (all(nations %in% names(x))) {
    x$uk <- sircovid::combine_rt(x[nations], samples[nations],
                                 rank = FALSE, weight = weight)
  }
  x
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


combined_aggregate_variant_rt <- function(variant_rt, samples) {
  england <- sircovid::regions("england")
  nations <- sircovid::regions("nations")

  if (all(england %in% names(variant_rt))) {
    variant_rt$england <- combine_variant_rt(variant_rt[england],
                                             samples[england],
                                             rank = FALSE)
  }
  if (all(nations %in% names(variant_rt))) {
    variant_rt$uk <- combine_variant_rt(variant_rt[nations], samples[nations],
                                        rank = FALSE)
  }
  variant_rt
}


combine_variant_rt <- function(variant_rt, samples, rank) {

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

  what <- setdiff(names(variant_rt[[1]]), c("step", "date", "beta"))

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

  ## finally we join the variants together
  ret <- variant1
  for (w in what) {
    ret[[w]] <- aperm(abind_quiet(variant1[[w]], variant2[[w]], along = 3),
                         c(1, 3, 2))
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

reorder_variant_rt <- function(x, rank) {

  what <- setdiff(names(x), c("step", "date"))

  ## reorder_rt_ifr will only work for one variant so we will have to reorder
  ## each variant separately
  for (j in seq_len(2)) {
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
  steps_predict <- seq(samples$predict$step,
                       length.out = forecast_days + 1L,
                       by = samples$predict$rate)
  forecast <- spim_pmcmc_predict(
    samples, steps_predict,
    prepend_trajectories = FALSE)
  forecast$date <- forecast$step / forecast$rate

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
  samples$trajectories$step <- c(samples$trajectories$step,
                                 forecast_samples$trajectories$step[-1L])
  samples$trajectories$date <- c(samples$trajectories$date,
                                 forecast_samples$trajectories$date[-1L])
  samples$trajectories$predicted <-
    c(samples$trajectories$predicted,
      forecast_samples$trajectories$predicted[-1L])

  samples

}

## This is adapted from mcstate::pmcmc_predict
spim_pmcmc_predict <- function(object, steps, prepend_trajectories = FALSE,
                               n_threads = NULL, seed = NULL) {

  if (is.null(object$predict)) {
    stop("mcmc was run with return_state = FALSE, can't predict")
  }
  if (length(steps) < 2) {
    stop("At least two steps required for predict")
  }
  if (steps[[1]] != object$predict$step) {
    stop(sprintf("Expected steps[1] to be %d", object$predict$step))
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

  mod <- model$new(pars, steps[[1]], NULL, n_threads = n_threads,
                   seed = seed, pars_multi = TRUE)

  mod$update_state(state = state)
  if (!is.null(index)) {
    mod$set_index(index)
  }
  y <- mod$simulate(steps)

  res <- mcstate:::mcstate_trajectories(steps, object$predict$rate, y, TRUE)

  if (prepend_trajectories) {
    res <- mcstate:::bind_mcstate_trajectories(object$trajectories, res)
  }

  res
}
