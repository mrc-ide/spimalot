##' Port MTP outcomes to template
##'
##' @title Port MTP outcomes to template
##'
##' @param summary_tidy A tibble containing MTP projections
##'
##' @param date The date to start output from (not necessarily the model date)
##'
##' @param run_grid The grid of scenarios ran
##'
##' @param combined Output from combined fitting task
##'   rtm_inference_pmcmc_spim_fits2_combined
##'
##' @param spim_state_names Named vector of states to be extracted for SPI-M
##'
##'
##' @export
spim_mtp_summary_to_template <- function(summary_tidy, date, run_grid,
                                         combined, spim_state_names) {
  pop <- spim_mtp_population(combined)

  model_type <- combined$info[[1]]$model_type
  if (model_type == "BB") {
    output_str <- "Stochastic Compartmental Positivity"
  } else {
    output_str <- "Stochastic Compartmental Cases"
  }

  regions <- vapply(c(sircovid::regions("all"), "england", "uk"),
                    spimalot::spim_region_name, "")

  summary_tidy$state$spim_name <- summary_tidy$state$scenario

  ## create common template columns
  lapply(names(run_grid), mtp_template_common,
         date = date, model_type = output_str) %>%
    dplyr::bind_rows() %>%
    ## join to results
    dplyr::left_join(summary_tidy$state, by = c(Scenario = "spim_name")) %>%
    dplyr::filter(state %in% names(spim_state_names),
                  date >= !!date,
                  group == "all") %>%
    dplyr::mutate("Day of Value" = lubridate::day(date),
                  "Month of Value" = lubridate::month(date),
                  "Year of Value" = lubridate::year(date),
                  AgeBand = "All",
                  Geography = regions[as.character(region)],
                  ValueType = spim_state_names[as.character(state)],
                  quantile = as.numeric(gsub("%", "", quantile)) / 100,
                  value = if_else(state == "react_pos",
                                  value / pop[as.character(region)] * 100,
                                  value),
                  .after = "Creation Year") %>%
    dplyr::select(-c(scenario, vaccine_daily_doses, analysis,
                     date, state, region, group, beta_step,
                     booster_daily_doses, vaccine_status)) %>%
    tidyr::pivot_wider(names_from = quantile, values_from = value,
                       names_prefix = "Quantile ") %>%
    dplyr::mutate(Value = `Quantile 0.5`, .after = ValueType)
}


## TODO: this should be harmonised with similar code in combined_spim.R
mtp_template_common <- function(scenario, date, model_type) {
  version <- packageVersion("sircovid")
  version <- substr(version, 3, nchar(version))

  data_frame(Group = "Imperial",
             Model = model_type,
             Scenario = scenario,
             ModelType = "Multiple",
             Version = version,
             "Creation Day" = lubridate::day(date),
             "Creation Month" = lubridate::month(date),
             "Creation Year" = lubridate::year(date))
}

##' Generate MTP population from combined
##'
##' @title Generate MTP population from combined
##'
##' @param combined Output from combined fitting task
##'   rtm_inference_pmcmc_spim_fits2_combined
##'
##' @export
spim_mtp_population <- function(combined) {
  pop <- vnapply(combined$pars[1, ], function(x) sum(x$N_tot[2:18]))
  c(pop,
    england = sum(unlist(pop[sircovid::regions("england")])),
    uk = sum(unlist(pop[sircovid::regions("all")])))
}


##' Prepare outputs by age and vaccination class for plotting
##'
##' @title MTP simulation outputs by age and vaccination class
##'
##' @param res Output from MTP simulation
##'
##' @param region A string, region for which outputs will be plotted
##'
##'
##' @export
spim_mtp_age_vaccine_outputs <- function(res, region = "england") {

  ## Objects for saving list of plots and matrices
  scenario_plots <- NULL
  scenario_plots_proportion <- NULL
  scenario_matrices <- NULL

  res <- dplyr::filter(res, region == !!region,
                       group != "all",
                       state %in% c("infections", "diagnoses_admitted",
                                    "deaths"))

  for (s in unique(res$scenario)) {

    tmp <- dplyr::filter(res, scenario == s)
    plots_age_vacc <- NULL
    plots_age_vacc_prop <- NULL


    for (w in unique(res$state)) {

      plot_matrix <- tmp %>%
        dplyr::filter(state == w)  %>%
        dplyr::mutate(age = factor(group,
                                   levels = unique(group),
                                   labels = c("Under 30s", "30 to 49",
                                              "50 to 74", "75+")))
      plots_age_vacc[[w]] <- ggplot2::ggplot(
        plot_matrix,
        ggplot2::aes(date, value, fill = age)) +
        ggplot2::ylab(paste(stringr::str_to_sentence(w))) + ggplot2::xlab("") +
        ggplot2::geom_area() + ggplot2::theme_bw() +
        ggplot2::facet_wrap(vars(vaccine_status)) +
        ggsci::scale_fill_lancet() +
        ggplot2::theme(axis.title.y = element_text(size = rel(0.9)),
                       legend.position = "none",
                       axis.title.x = element_text(size = rel(0.8)),
                       legend.title = element_blank(),
                       strip.text.x = element_text(size = rel(0.7)))

      plots_age_vacc_prop[[w]] <-
        ## Re-arrange matrix for proportional stacked chart
        plot_matrix %>%
        dplyr::group_by(date, vaccine_status, age) %>%
        dplyr::summarise(n = sum(value)) %>%
        dplyr::mutate(percentage = n / sum(n)) %>%
        ggplot2::ggplot(., ggplot2::aes(x = date, y = percentage, fill = age)) +
        ggplot2::geom_area(alpha = 0.6, size = 1, colour = "black") +
        ggplot2::ylab("") + ggplot2::xlab("") +
        ggplot2::geom_area() + ggplot2::theme_bw() +
        ggplot2::facet_wrap(vars(vaccine_status)) +
        ggsci::scale_fill_lancet() +
        ggplot2::theme(axis.title.y = element_text(size = rel(0.9)),
                       axis.title.x = element_text(size = rel(0.8)),
                       legend.title = element_blank(),
                       strip.text.x = element_text(size = rel(0.7)))

      if (w != "diagnoses_admitted_inc") {
        plots_age_vacc_prop[[w]] <- plots_age_vacc_prop[[w]] +
          ggplot2::theme(legend.position = "none")
        }

      scenario_matrices[[s]][[w]] <- plot_matrix
    }

    ## Save plot and matrix object into lists

    scenario_plots[[s]] <-
      (plots_age_vacc[["infections"]] +
         plots_age_vacc_prop[["infections"]]) /
      (plots_age_vacc[[ "diagnoses_admitted"]] +
         plots_age_vacc_prop[[ "diagnoses_admitted"]]) /
      (plots_age_vacc[["deaths"]] + plots_age_vacc_prop[["deaths"]]) +
      plot_annotation(
        caption = paste("Partial immunity includes those with either one or",
                        "two doses that have not yet achieved full vaccine",
                        "efficacy"),
        title = paste(stringr::str_to_sentence(region)),
        subtitle = paste(stringr::str_to_sentence(sub("_", " ", s))))

  }

  ## Main scenario always first element in list
  main_scenario <- grep("main", names(scenario_matrices))
  other_scenarios <- grep("main", names(scenario_matrices), invert = TRUE)
  scenario_plots <- scenario_plots[c(main_scenario, other_scenarios)]
  scenario_matrices <- scenario_matrices[c(main_scenario, other_scenarios)]

  out <- NULL
  out[["plots"]] <- scenario_plots
  out[["matrices"]] <- scenario_matrices
  out
}
