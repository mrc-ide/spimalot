##' Plot Rt distribution over time for various scenarios
##' @title Plot Rt distribution
##'
##' @param npi_key data.frame type object with columns:
##'   nations, npi (scenario), Rt (mean), Rt_sd
##' @param xlim Passed to `plot`
##' @param ylim Passed to `plot`
##' @param cols Colours for plots, should be same for `npi_key2` if given
##' @param labels Passed to [legend] `legend` parameter
##' @param legend_ncol Passed to [legend] `ncol` parameter
##' @param npi_key2 Optional second npi_key that will be plotted with dashed
##'   lines
##' @param multiplier Multiplier on the R values to plot - the distribution
##' plotted will be that of multiplier * x where x is drawn from a lognormal
##' with parameters specified in npi_key
##'
##' @export
spim_plot_rt_dist <- function(npi_key, xlim, ylim, cols, labels = NULL,
                              legend_ncol = 1, npi_key2 = NULL,
                              multiplier = 1) {

  x <- seq(min(xlim), max(xlim), length.out = 1e3)

  labels <- labels %||% rownames(npi_key)
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n",
       xlab = expression(R[excl_immunity]),
       ylab = "Density",
       las = 1)

  for (i in seq_rows(npi_key)) {
    dist <- distr6::dstr("Lognormal", mean = npi_key$Rt[i],
                         sd = npi_key$Rt_sd[i])
    lines(x * multiplier, dist$pdf(x) / multiplier, col = cols[i], lwd = 2)
    legend_lty <- rep(1, nrow(npi_key2))
  }
  if (!is.null(npi_key2)) {
    for (i in seq_rows(npi_key)) {
      dist <- distr6::dstr("Lognormal", mean = npi_key2$Rt[i],
                           sd = npi_key2$Rt_sd[i])
      lines(x * multiplier, dist$pdf(x) / multiplier,
            col = cols[i], lwd = 2, lty = 2)
    }
    legend_lty <- c(legend_lty, 1, 2)
    cols <- c(cols, 1, 1)
  }

  legend("topleft", legend = labels, ncol = legend_ncol, bty = "n",
         lty = legend_lty, lwd = 2, col = cols)
}


##' Plot seasonality trends over time
##' @title Plot seasonality over time
##'
##' @param peak_date Date where seasonal multiplier is highest
##' @param seasonality Seasonal multiplier
##'
##' @export
spim_plot_seasonality <- function(peak_date = as.Date("2020-02-15"),
                                  seasonality = 0.1) {
  x <- seq_len(365)
  dx <- sircovid::sircovid_date_as_date(x)
  plot(dx,
       spim_calc_seasonality(x, sircovid::sircovid_date(peak_date),
                             seasonality),
       type = "l", lwd = 2, ylim = c(0.9, 1.1),
       ylab = "Seasonal multiplier", xlab = "", xaxt = "n")
  axis.Date(1, dx, at = seq.Date(from = dx[1], to = as.Date("2021-01-01"),
                                 by = "1 month"))
  abline(v = c(peak_date, round(peak_date + 365 / 2)), lty = 2,
         col = "grey30")
}


##' Plot relative VOC mean and sd
##' @title Plot relative strain transmissibility ranges
##'
##' @param R1 Lognormal mean for strain 1
##' @param R1_sd Lognormal standard deviation for strain 1
##' @param epsilon_range Relative range for strain 2
##' @param epsilon_central Relative mean for strain 2
##'
##' @export
spim_plot_voc_range <- function(R1, R1_sd, epsilon_range, epsilon_central) {
  NPI <- central <- variant <- NULL

    R1_range <- as.list(
        distr6::dstrs(
            "Lognormal",
            data.frame(mean = R1, sd = R1_sd)
        )$quantile(c(0.025, 0.975))
    )

    R2_range <- lapply(seq_along(R1_range), function(i) {
        c(min(R1_range[[i]]) * min(epsilon_range),
          max(R1_range[[i]]) * max(epsilon_range))
    })

    data_frame(do.call(rbind, Map(rbind, R1_range, R2_range))) %>%
       `colnames<-`(c("min", "max")) %>%
        dplyr::mutate(
            central = c(R1[1], (R1 * epsilon_central)[1], R1[2],
                        (R1 * epsilon_central)[2]),
            variant = rep(c("Alpha", "Delta"), 2),
            NPI = c(rep("central R after NPI lift", 2),
                    rep("high R after NPI lift", 2)),
            NPI = factor(NPI, levels = unique(NPI))
        ) %>%
        ggplot2::ggplot() +
        geom_pointrange(aes(x = NPI, y = central, ymin = min, ymax = max,
                            colour = variant),
            size = 1, position = position_dodge(width = 0.2)
        ) +
        scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 11)) +
        theme_classic() +
        xlab(NULL) +
        ylab("R excluding immunity") +
        theme(
            text = element_text(size = 14), legend.position = c(0.2, 0.9),
            legend.title = element_blank()
        )
}

##' Select colours for plotting simulation scenarios
##' @title Return accessible scenario colours
##'
##' @param scenarios Unique scenario names. Scenarios with `[High R]` will
##'   darken the colour for the corresponding central scenario and scenarios
##'   with `[Low R]` will brighten the colour for the corresponding
##'   central scenario.
##' @param weight `weight` passed to `mix_cols` to darken/brighten colours for
##'   high/low R scenarios
##' @param palette Colour palette, passed to [khroma::colour]
##' @param preview If `TRUE` then plots the final colour scheme with
##'   [khroma::plot_scheme]
##'
##' @examples
##' spim_scenario_cols(c("Step 4", "Step 4 [High R]", "Step 4 [Low R]"))
##'
##' spim_scenario_cols(c("Step 3", "Step 3 [High R]", "Step 4"))
##'
##' @export
spim_scenario_cols <- function(scenarios, weight = 0.3, palette = "bright",
                               preview = FALSE) {

  stopifnot(all(table(scenarios)) == 1)

  cen_scenarios <- scenarios[!grepl("(High|Low) R", scenarios)]
  n_scens <- length(cen_scenarios)

  cols <- khroma::colour(palette)(n_scens)
  names(cols) <- cen_scenarios
  dark_cols <- light_cols <- character()

  dark_scenarios <- grep("[High R]", scenarios, fixed = TRUE, value = TRUE)
  if (length(dark_scenarios) > 0) {
    dark_cols <- match(gsub(" [High R]", "", dark_scenarios, fixed = TRUE),
                       names(cols))

    if (any(is.na(dark_cols))) {
      stop("Unrecognised High R scenarios")
    }

    dark_cols <- mix_cols(cols[dark_cols],
                          rep("#000000", length(dark_cols)), weight)
    names(dark_cols) <- dark_scenarios
  }

  light_scenarios <- grep("[Low R]", scenarios, fixed = TRUE, value = TRUE)
  if (length(light_scenarios) > 0) {
    light_cols <- match(gsub(" [Low R]", "", light_scenarios, fixed = TRUE),
                       names(cols))

    if (any(is.na(light_cols))) {
      stop("Unrecognised Low R scenarios")
    }

    light_cols <- mix_cols(cols[light_cols],
                           rep("#FFFFFF", length(light_cols)), weight)
    names(light_cols) <- light_scenarios
  }

  cols <- c(cols, dark_cols, light_cols)

  if (preview) {
    khroma::plot_scheme(cols)
  }

  cols
}


##' Prepare aggregated real data for plotting
##' @title Prepare aggregated data for plotting
##'
##' @param path Path to aggregated rds object containing a named list where
##'  names correspond to regions and each element is a list with `full` data
##'  and `fitted` data
##'
##' @return Returns a list where elements correspond to regions and `fitted`
##'  data is removed. Data processing includes: adding `deaths`
##'  column as the sum over all death compartments; fixes when NAs converted to
##'  0s erroneously; and fitted data removed.
##'
##' @export
spim_prepare_aggregated_data <- function(path) {

  agg_data <-
    readRDS(path) %>%
    lapply(function(x) {
      x <- x$full

      deaths <- cbind(x$deaths_hosp, x$deaths_comm,
                      x$deaths_carehomes, x$deaths_non_hosp)

      x$deaths <- rowSums(deaths, na.rm = TRUE)
      x$deaths[apply(deaths, 1, function(x) all(is.na(x)))] <- NA
      x
    })

  nr <- nrow(agg_data[[1]])

  f <- function(what) {
    mat <- vapply(agg_data[sircovid::regions("england")],
                  function(x) x[[what]], integer(nr))
    which <- apply(mat, 1, function(x) all(is.na(x)))
    mat <- rowSums(mat, na.rm = TRUE)
    mat[which] <- NA
    mat
  }

  what <- c("icu", "general", "hosp",  "admitted", "diagnoses", "all_admission")
  agg_data$england[, what] <- vapply(what, f, numeric(nr))
  agg_data
}


##' Diagnostic plots for checking Rt from simulation
##' @title Check Rt from simulation
##' @param summary_state State from simulation summary object
##' @param dates Dates to plot vertical lines at
##' @param combined_state Optional state from combined object
##' @param Rt_state Rt states to plot, default is effective average Rt
##' @export
spim_plot_check_rt <- function(summary_state, dates, combined_state = NULL,
                              Rt_state = "eff_Rt_general_both") {
  region <- state <- scenario <- analysis <- `50%` <- NULL

  summary_state <- summary_state %>%
    dplyr::filter(region == "england",
                  state %in% Rt_state) %>%
    dplyr::mutate(scenario = as.factor(scenario),
                  analysis = as.factor(analysis)) %>%
    tidyr::pivot_wider(names_from = quantile)

  ddates <- rep(as.Date(Sys.time()), length(dates))
  ddates[dates != "today"] <- as.Date(dates[dates != "today"])

  ci_limits <- summary_quantile_range(summary_state)

  p <- summary_state %>%
    ggplot(aes(x = date, y = `50%`, colour = scenario)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = get(ci_limits$ymin),
                    ymax = get(ci_limits$ymax),
                    fill = "scenario")) +
    geom_line() +
    facet_grid(cols = vars(analysis), rows = vars(state),
              scales = "free_y",
              labeller = label_wrap_gen(width = 7)) +
    scale_x_date(date_breaks = "1 month") +
    geom_vline(xintercept = ddates, lty = 2, color = "gray")

  if (!is.null(combined_state)) {
      combined_state <- combined_state %>%
        dplyr::filter(
          region == "england",
          state %in% Rt_state) %>%
        tidyr::pivot_wider(names_from = quantile)

      p <- p +
        geom_line(data = combined_state, aes(x = date, y = `50%`),
                  inherit.aes = FALSE)
  }

  p
}


##' Diagnostic plots for checking states from simulation
##' @title Check states from simulation
##' @param summary_state State from simulation summary object
##' @param combined_state Optional state from combined object
##' @export
spim_plot_check_state <- function(summary_state, combined_state = NULL) {
  region <- scenario <- analysis <- state <- group <- `50%` <- NULL

  summary_state <- summary_state %>%
    dplyr::filter(region == "england") %>%
    dplyr::mutate(scenario = as.factor(scenario),
                  analysis = as.factor(analysis)) %>%
    dplyr::filter(state %in% c("diagnoses_admitted_inc",
                              "deaths_inc", "infections_inc",
                              "hosp"),
                  group == "all") %>%
    tidyr::pivot_wider(names_from = quantile)

  ci_limits <- summary_quantile_range(summary_state)

  # TODO:: for SPI-M we report 90% CI, hence 5% and 95% here but
  # default below; shall we change default?
  p <- summary_state %>%
    ggplot(aes(x = date)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = get(ci_limits$ymin),
                    ymax = get(ci_limits$ymax),
                    fill = scenario)) +
    geom_line(aes(y = `50%`, color = scenario)) +
    facet_grid(rows = vars(state), cols = vars(analysis), scales = "free",
              labeller = label_wrap_gen(width = 5)) +
    scale_x_date(date_breaks = "1 month")

    if (!is.null(combined_state)) {
      combined_state <- dplyr::filter(combined_state, region == "england") %>%
        dplyr::filter(state %in% c("diagnoses_admitted_inc",
                                  "deaths_inc", "infections_inc",
                                  "hosp"),
                      group == "all") %>%
        tidyr::pivot_wider(names_from = quantile)

        p <- p +
              geom_line(data = combined_state, aes(x = date, y = `50%`),
                        inherit.aes = FALSE)
    }

    p
}


##' Diagnostic plots for checking states by age from simulation
##' @title Check states by age from simulation
##' @param summary_agestate State by age from simulation summary object
##' @param ana Analysis to check, usually central
##' @param scen Scenario to check, usually central
##' @export
spim_plot_check_state_by_age <- function(summary_agestate, ana, scen) {
  region <- analysis <- scenario <- value <- vaccine_status <- NULL
  state <- group <- NULL

  summary_agestate %>%
    dplyr::filter(
      region == "england",
      analysis == ana,
      scenario == scen
    ) %>%
    ggplot(aes(x = date, y = value, fill = vaccine_status)) +
    geom_area() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(vars(state), vars(group),
      scales = "free",
      labeller = label_wrap_gen(width = 10)
    ) +
    scale_x_date(date_breaks = "1 month")
}


#' Calculate doses given out from simulation
#' @title Calculate doses given out from simulation
#' @param summary Output from [spim_simulate_tidy_states]
#' @param population_england Vector of England population by 19 age groups
#'  and rows are age groups
#' @param scen Scenario to check, usually central
#' @export
spim_calculate_doses <- function(summary, population_england, scen) {
  ## Declare variables in scope so that R CMD check does not warn us
  region <- scenario <- state <- group <- value <- pop <- NULL
  state_first_dose_inc <- NULL
  state_second_dose_inc <- NULL
  state_booster_dose_inc <- NULL

  doses_g <- summary$n_doses %>%
    dplyr::filter(region == "england",
                  scenario == scen) %>%
    tidyr::pivot_wider(names_from = state, names_prefix = "state_") %>%
    dplyr::mutate(state_total_dose_inc = state_first_dose_inc +
                    state_second_dose_inc +
                    state_booster_dose_inc) %>%
    tidyr::pivot_longer(dplyr::starts_with("state_"), names_to = "state") %>%
    tidyr::pivot_wider(names_from = group, names_prefix = "group_") %>%
    dplyr::mutate(group_total = rowSums(dplyr::across(starts_with("group")),
                                        na.rm = TRUE)) %>%
    tidyr::pivot_longer(starts_with("group_"), names_to = "group")

  pop_df <- data.frame(group = unique(doses_g$group),
                        pop = c(population_england,
                                total = sum(population_england)))

  doses_g %>%
    dplyr::left_join(pop_df) %>%
    dplyr::mutate(prop = value / pop)
}


#' Diagnostic plots for checking doses given out from simulation over
#'  all analyses and age groups
#' @title Check doses given out from simulation
#' @param doses Output from [spim_calculate_doses]
#' @export
#' @import ggplot2
spim_plot_check_doses <- function(doses) {
  value <- group <- NULL
  doses %>%
    ggplot(aes(x = date, y = value, colour = group)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_line() +
    facet_grid(rows = vars(state), cols = vars(analysis), scales = "free",
                labeller = label_wrap_gen(width = 10)) +
    scale_x_date(date_breaks = "1 month")
}


#' Diagnostic plots for checking total doses given out in simulation
#' @title Check total doses given out from simulation
#' @param doses Output from [spim_calculate_doses]
#' @export
spim_plot_check_total_doses <- function(doses) {
  group <- state <- value <- analysis <- NULL
  doses %>%
    dplyr::filter(group == "group_total",
                  state == "state_total_dose_inc") %>%
    ggplot(aes(x = date, y = value * 7, colour = group)) +
    theme_bw() +
    ylab("Weekly doses") +
    geom_line() +
    facet_wrap(vars(analysis)) +
    theme(legend.position = "n") +
    scale_x_date(date_breaks = "1 month")
}


#' Diagnostic plots for checking dose uptake from simulation
#' @title Check dose uptake from simulation
#' @param doses Output from [spim_calculate_doses]
#' @export
spim_plot_check_uptake <- function(doses) {
  state <- prop <- group <- analysis <- NULL
  doses %>%
    dplyr::filter(!grepl("_inc", state)) %>%
    ggplot(aes(x = date, y = prop, colour = group)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_line() +
    facet_grid(rows = vars(state), cols = vars(analysis), scales = "free",
              labeller = label_wrap_gen(width = 10)) +
    scale_x_date(date_breaks = "1 month")
}


#' Diagnostic plots for checking daily infections and proportion of strain 2
#' @title Check daily infections from simulation
#' @param summary_state State from simulation summary object
#' @param scen Scenario to plot
#' @param ana Analysis to plot
#' @export
spim_plot_daily_infections <- function(summary_state, ana, scen) {
  state <- region <- scenario <- analysis <- `50%` <- NULL

  summary_state %>%
    dplyr::filter(
      state %in% c("n_strain_1_inc", "n_strain_2_inc", "prop_strain_2"),
      region == "england",
      scenario == scen,
      analysis == ana
    ) %>%
    dplyr::mutate(
      state = factor(
        state,
        labels = c(
          "Daily infections strain 1", "Daily infections strain 2",
          "Proportion of strain 2 among daily infections"
        )
      )
    ) %>%
    tidyr::pivot_wider(names_from = quantile) %>%
    ggplot(aes(x = date, y = `50%`)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_blank()) +
    geom_line() +
    facet_grid(rows = vars(state), scales = "free",
              labeller =  label_wrap_gen(width = 30)) +
    scale_x_date(date_breaks = "1 month")
}

summary_quantile_range <- function(sum_df) {
  if ("2.5%" %in% colnames(sum_df)) {
    ymin <- "2.5%"
    ymax <- "97.5%"
  } else {
    ymin <- "5%"
    ymax <- "95%"
  }
  list(
    ymin = ymin,
    ymax = ymax
  )
}
