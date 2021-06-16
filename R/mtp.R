##' Prepare current Rt values for MTP simulations
##'
##' @title Prepare Rt for MTP scenarios
##'
##' @param rt Outputs for Rt by region from rtm_inference_spim_fits2_combined
##'
##' @param type A string or vector of Rt type to use, can be one of
##'   "eff_Rt_general", "Rt_general", "eff_Rt_all" or "Rt_all".
##'
##'
##' @export
spim_mtp_rt <- function(rt, date, type = c("eff_Rt_general", "Rt_general")){
  ret <- NULL
  regions <- c(sircovid::regions("all"), "england", "uk")

  for (t in type) {
    for (r in regions) {
      ret[[t]][[r]] <-
        rt[[r]][[t]][rt[[r]]$date == sircovid::sircovid_date(date)]
    }
  }
  ret
}

##' Prepare NPIs effect on Rt for MTP simulations
##'
##' @title Prepare NPIs for MTP scenarios
##' @param npi_key A CSV file containing NPI parameters
##'
##' @param Rt Rt values for all regions and nations obtained via
##'   [spimalot::spim_mtp_rt]
##'
##' @param rt_reduction A double specifying a reduction on Rt (excluding
##'   immunity) given school closures. A default value of 0.5 roughly yields
##'   a 0.3 reduction in eff_Rt.
##'
##'
##' @export
spim_mtp_npi_key <- function(npi_key, Rt, rt_reduction = 0.5){

  curr_Rt <- Rt$Rt_general
  curr_eff_Rt <- Rt$eff_Rt_general

  npi_key$Rt <- round(npi_key$Reff_t * mean(curr_Rt[["uk"]] /
                                              curr_eff_Rt[["uk"]]), 1)

  npi_key["schools_main", "Rt"] <- round(mean(curr_Rt[["uk"]]), 1)

  npi_key_holidays <- npi_key
  npi_key_holidays$Rt <- npi_key$Rt - rt_reduction
  rownames(npi_key_holidays) <- gsub("schools", "holidays", rownames(npi_key))

  npi_key <- rbind(npi_key, npi_key_holidays)
}

##' Prepare scenario simulations for MTP commission
##'
##' @title Prepare MTP scenarios
##' @param mtp_commission A CSV specifying MTP scenarios for SPI-M
##'
##' @param npi_key A CSV with corresponding effect of NPIs on Rt values
##'
##' @param end_date End date for the simulation, this should usually be
##'   current date + 9 weeks
##'
##' @param vaccine_parameters Output from upstream task rtm_vaccination_parameters
##'
##' @param combined Output from combined fitting task
##'   rtm_inference_pmcmc_spim_fits2_combined
##'
##' @param rt Rt values for all regions and nations obtained via
##'   [spimalot::spim_mtp_rt]
##'
##' @param n_par Number of particles for simulation
##'
##' @param end_date Date to end simulation
##'
##' @param n_threads Number of threads to run the simulation on
##'
##' @export
##'
spim_mtp_prepare <- function(mtp_commission, combined,
                             vaccination_parameters, npi_key, rt,
                             n_par, end_date, n_threads) {
  date <- combined$date

  ## get current Rt for all regions and nations
  curr_Rt <- spim_mtp_rt(rt, date)

  ## get NPIs key for simulation
  npi_key_expand <- spim_mtp_npi_key(npi_key, curr_Rt)

  ## expand data frame to include all regions
  f <- function(df) {
    msg <- setdiff(names(combined$pars), df$region)
    if (length(msg) > 0) {
      f <- function(r) {
        d <- df[df$region == "combined", ]
        d$region <- r
        d
      }
      df <- dplyr::bind_rows(c(list(df), lapply(msg, f)))
    }
    df[df$region != "combined", ]
  }

  mtp_commission <- split(mtp_commission, mtp_commission$scenario)
  ret <- lapply(mtp_commission, f) %>%
    dplyr::bind_rows()


  run_grid <- ret %>%
    dplyr::select(scenario, daily_doses, spim_name, type_rt) %>%
    dplyr::distinct()

  rt_future <- ret %>%
    dplyr::select(-c(daily_doses, spim_name)) %>%
    tidyr::pivot_longer(c(starts_with("npi"), starts_with("date")),
                        names_to = c(".value", "changepoint"),
                        names_pattern = "(.+)_(.+)") %>%
    dplyr::mutate(Rt_sd = npi_key_expand[npi, "Rt_sd"],
                  Rt = npi_key_expand[npi, "Rt"],
                  .after = npi)

  for (r in unique(rt_future$region)){
    pivot <- rt_future$pivot[rt_future$region == r]
    if (pivot == "scaling"){
      type_rt <- rt_future$type_rt[rt_future$region == r]
      base_Rt <- mean(curr_Rt[[type_rt]][[r]])
      rt_future$Rt[rt_future$region == r] = base_Rt +
        rt_future$Rt_1[rt_future$region == r]
    }
  }

  rt_future <- split(rt_future, f = rt_future$scenario)
  rt_future <- rt_future[run_grid$scenario]

  ## Things that we want to keep:
  inc_states <- c("deaths", "admitted", "diagnoses", "infections")
  prev_states <- c("icu", "general", "hosp", "react_pos")
  output <- list(incidence = inc_states,
                 prevalence = prev_states,
                 keep = c(inc_states, prev_states))

  ## Vaccination parameters we need
  future_daily_doses <- lapply(vaccination_parameters$daily_doses, "[[", "mtp")
  vaccine_efficacy <- vaccination_parameters$vacc_efficacy$central
  vaccine_uptake <- vaccination_parameters$uptake_by_age$central
  eligibility_by_age <- vaccination_parameters$eligibility_by_age$min_18

  list(combined = combined,
       curr_Rt = curr_Rt,
       npi_key = npi_key_expand,
       run_grid = run_grid,
       rt_future = rt_future,
       n_par = n_par,
       end_date = end_date,
       n_threads = n_threads,
       n_run = nrow(run_grid),
       future_daily_doses = future_daily_doses,
       vaccine_efficacy = vaccine_efficacy,
       vaccine_uptake = vaccine_uptake,
       eligibility_by_age = eligibility_by_age,
       output = output)
}


##' Prepare MTP outcomes template
##'
##' @title Prepare MTP outcomes template
##'
##' @param scenario A string for the name of the scenario requested by SPI-M
##'
##' @param date Date of the data
##'
##' @param model_type Model used
##'
##'
##' @export
spim_mtp_template_common <- function(scenario, date, model_type) {

  version <- packageVersion("sircovid")
  version <- substr(version, 3, nchar(version))

  data.frame(Group = "Imperial",
             Model = model_type,
             Scenario = scenario,
             ModelType = "Multiple",
             Version = version,
             "Creation Day" = lubridate::day(date),
             "Creation Month" = lubridate::month(date),
             "Creation Year" = lubridate::year(date),
             stringsAsFactors = FALSE,
             check.names = FALSE)
}


##' Port MTP outcomes to template
##'
##' @title Port MTP outcomes to template
##'
##' @param summary_tidy A tibble containing MTP projections
##'
##' @param run_grid The grid of scenarios ran
##'
##' @param date A string for date of simulation outcomes
##'
##' @param combined Output from combined fitting task
##'   rtm_inference_pmcmc_spim_fits2_combined
##'
##' @param model_type Model used
##'
##' @param spim_state_names Named vector of states to be extracted for SPI-M
##'
##'
##' @export
spim_mtp_summary_to_template <- function(summary_tidy, run_grid, date, regions,
                                         combined, model_type,
                                         spim_state_names) {

  pop <- spim_mtp_get_population(combined)

  if (model_type == "BB"){
    output_str <- "Stochastic Compartmental Positivity"
  } else {
    output_str <- "Stochastic Compartmental Cases"
  }

  ## create common template columns
  lapply(run_grid$spim_name, spim_mtp_template_common,
         date = date, model_type = model_type) %>%
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
                  Geography = regions[as.character(region), "name"],
                  ValueType = spim_state_names[as.character(state)],
                  quantile = as.numeric(gsub("%", "", quantile)) / 100,
                  value = if_else(state == "react_pos",
                                  value / pop[as.character(region)] * 100,
                                  value),
                  .after = "Creation Year") %>%
    dplyr::select(-c(scenario, daily_doses, type_rt, date, state, region, group,
                     vaccine_status)) %>%
    tidyr::pivot_wider(names_from = quantile, values_from = value,
                       names_prefix = "Quantile ") %>%
    dplyr::mutate(Value = `Quantile 0.5`, .after = ValueType)
}


##' Get population by region for MTP outcomes to template
##'
##' @title Get population by region
##'
##' @param combined Output from combined fitting task
##'   rtm_inference_pmcmc_spim_fits2_combined
##'
##'
##' @export
spim_mtp_get_population <- function(combined) {
  pop <- vnapply(names(combined$pars), function(r)
    sum(combined$transform[[r]](combined$pars[[r]][1, ])$N_tot[2:18]))
  names(pop) <- names(combined$pars)
  pop_england <- sum(unlist(pop[sircovid::regions("england")]))
  pop_uk <- sum(unlist(pop[sircovid::regions("all")]))
  c(pop, england = pop_england, uk = pop_uk)
}



##' Extract MTP simulation outputs into format required by SPI-M
##'
##' @title Save MTP outputs for SPI-M
##'
##' @param dat Output from MTP simulation
##'
##' @param path_template An xlsx template as shared by SPI-M
##'
##' @param path_save A string, name of file to save in xlsx format
##'
##' @param root Folder for outputs to be saved in
##'
##'
##' @export
spim_mtp_save_results <- function(dat,
                                  path_template = "template_combined.xlsx",
                                  path_save = "Imperial_MTP.xlsx",
                                  root = "outputs") {
  sheets <- readxl::excel_sheets(path_template)
  template <- setNames(
    lapply(sheets, readxl::read_xlsx, path = path_template,
           .name_repair = "minimal"),
    sheets)
  stopifnot(identical(names(template[[2]]), names(dat)))
  template[[2]] <- dat
  writexl::write_xlsx(template, sprintf("%s/%s", root, path_save))
}


##' Prepare outputs by age and vaccination class for plotting
##'
##' @title MTP simulation outputs by age and vaccination class
##'
##' @param dat Output from MTP simulation
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
                       state %in% c("infections_inc", "diagnoses_admitted_inc",
                                    "deaths_inc", "deaths"))

  for (s in unique(res$scenario)) {

    tmp <- dplyr::filter(res, scenario == s)
    plots_age_vacc <- NULL
    plots_age_vacc_prop <- NULL


    for (w in unique(res$state)) {

      plot_matrix <- tmp %>%
        dplyr::filter(state == w)  %>%
        dplyr::mutate(age = factor(group,
                                   levels = unique(group),
                                   labels = c("Under 30s", "30 to 49", "50 to 74",
                                              "75+")))
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
        ggplot2::ggplot(., ggplot2::aes(x=date, y=percentage, fill=age)) +
        ggplot2::geom_area(alpha=0.6 , size=1, colour="black") +
        ggplot2::ylab("") + ggplot2::xlab("") +
        ggplot2::geom_area() + ggplot2::theme_bw() +
        ggplot2::facet_wrap(vars(vaccine_status)) +
        ggsci::scale_fill_lancet() +
        ggplot2::theme(axis.title.y = element_text(size = rel(0.9)),
                       axis.title.x = element_text(size = rel(0.8)),
                       legend.title = element_blank(),
                       strip.text.x = element_text(size = rel(0.7)))

      if (w != "diagnoses_admitted_inc"){
        plots_age_vacc_prop[[w]] <- plots_age_vacc_prop[[w]] +
          ggplot2::theme(legend.position = "none")}

      scenario_matrices[[s]][[w]] <- plot_matrix
    }

    ## Save plot and matrix object into lists

    scenario_plots[[s]] <-
      (plots_age_vacc[["infections_inc"]] + plots_age_vacc_prop[["infections_inc"]]) /
      (plots_age_vacc[[ "diagnoses_admitted_inc"]] + plots_age_vacc_prop[[ "diagnoses_admitted_inc"]]) /
      (plots_age_vacc[["deaths_inc"]] + plots_age_vacc_prop[["deaths_inc"]]) /
      (plots_age_vacc[["deaths"]] + plots_age_vacc_prop[["deaths"]]) +
      plot_annotation(
        caption = "Partial immunity includes those with either one or two doses that have not yet achieved full vaccine efficacy",
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


##' Run MTP simulation
##'
##' @title Run MTP simulation
##'
##' @param obj Large object containing parameters for MTP simulation as per
##'   [spimalot::spim_mtp_prepare]
##'
##' @param simulate_schedule Global simulate function, this needs porting over
##'   to spimalot
##'
##'
##' @export
spim_mtp_simulate <- function(obj, simulate_schedule) {

  n_run <- nrow(obj$run_grid)
  res <- vector("list", n_run)
  i <- 1L
  for (i in seq_len(n_run)) {
    message(sprintf("Running scenario %d / %d", i, n_run))
    res[[i]] <- simulate_schedule(
      combined = obj$combined,
      n_par = obj$n_par,
      end_date = obj$end_date,
      n_threads = obj$n_threads,
      keep = obj$output$keep,
      rt_future = obj$rt_future[[i]],
      rt_type = obj$run_grid$type_rt[i],
      future_daily_doses = obj$future_daily_doses,
      vaccine_efficacy = obj$vaccine_efficacy,
      uptake_by_age = obj$vaccine_uptake,
      eligibility_by_age = obj$eligibility_by_age[i],
      calculate_rt = TRUE,
      n_strain = 1)
  }
  names(res) <- obj$run_grid$scenario

  res
}
