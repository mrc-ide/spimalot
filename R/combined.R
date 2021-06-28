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
##' @return A combined fit object
##' @export
spim_combined_load <- function(path, regions = "all") {
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
    x$samples$info[c("date", "multistrain", "model_type")])
  if (length(unique(info)) != 1L) {
    stop("Incompatible regional fits")
  }

  ret <- list_transpose(dat)
  ret$info <- info[[1]]
  ret$info$date <- as.Date(ret$info$date)

  message("Reordering trajectories")
  ## reorder by increasing cumulative incidence:
  rank_cum_inc <- lapply(ret$samples, sircovid::get_sample_rank)
  ret$samples <- Map(sircovid::reorder_sample, ret$samples, rank_cum_inc)
  ret$rt <- Map(sircovid::reorder_rt_ifr, ret$rt, rank_cum_inc)
  ret$ifr_t <- Map(sircovid::reorder_rt_ifr, ret$ifr_t, rank_cum_inc)
  if (ret$info$multistrain) {
    ret$variant_rt <- Map(reorder_variant_rt, ret$variant_rt, rank_cum_inc)
  }

  ## NOTE: have not ported the "randomise trajectory order" bit over,
  ## but I do not think that we need to.
  message("Creating data for onward use")
  ## Do this step *before* aggregation because otherwise we end up
  ## with the onward data including england/uk
  ret$onward <- spim_combined_onward(ret)
  ret$parameters <- lapply(list_transpose(ret$parameters), dplyr::bind_rows)

  message("Aggregating England/UK")
  ## Aggregate some of these to get england/uk entries
  ret$samples <- combined_aggregate_samples(ret$samples)

  ret$data <- combined_aggregate_data(ret$data)
  ret$rt <- combined_aggregate_rt(ret$rt, ret$samples)
  ret$ifr_t <- combined_aggregate_rt(ret$ifr_t, ret$samples)
  if (ret$info$multistrain) {
    ret$variant_rt <- combined_aggregate_variant_rt(ret$variant_rt, ret$samples)
  }

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

  state <- lapply(dat$samples, function(x)
    x$trajectories$state[rownames(simulate$state[[1]]), , idx_dates])
  state <- aperm(abind::abind(state, along = 4), c(1, 2, 4, 3))

  state_by_age <- lapply(
    list_transpose(simulate$state_by_age),
    abind::abind, along = 3)

  ret <- list(date = dates,
              state = state,
              state_by_age = state_by_age,
              n_protected = abind::abind(simulate$n_protected, along = 2),
              n_doses = abind::abind(simulate$n_doses, along = 3))

  ## This is not terrible:
  rt <- list_transpose(dat$rt)[c("Rt_general", "eff_Rt_general")]
  ## Rt_general and eff_Rt_general will have dimensions:
  ## [n particles x n regions x n dates]
  rt_combined <- lapply(rt, function(x)
    aperm(abind::abind(x, along = 3), c(2, 3, 1))[, , idx_dates])

  ret <- c(ret, rt_combined)

  if (dat$info$multistrain) {
    idx_dates_variant <- dat$variant_rt[[1]]$date[, 1] %in% dates

    ## variant_Rt_general and variant_eff_Rt_general will have dimensions:
    ## [n particles x n regions x n variants x n dates]
    variant_rt <-
      list_transpose(dat$variant_rt)[c("Rt_general", "eff_Rt_general")]
    variant_rt_combined <- lapply(variant_rt, function(x)
      aperm(abind::abind(x, along = 4), c(3, 4, 2, 1))[, , , idx_dates_variant])
    names(variant_rt_combined) <- paste0("variant_", names(variant_rt_combined))

    ret <- c(ret, variant_rt_combined)
  }

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
    ret[-seq_len(nrow(ret) - 5), ] <- NA
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
  idx_dates <- samples[[1]]$trajectories$date %in% dates[-1L]
  idx_dates[1] <- TRUE
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
    ret[[w]] <- aperm(abind::abind(variant1[[w]], variant2[[w]], along = 3),
                         c(1, 3, 2))
  }

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
##' @param combined Combined fits object
##' @param ignore_uk If `TRUE` population of UK is retured as NA
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
