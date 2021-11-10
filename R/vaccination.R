##' Prepare vaccination data
##'
##' @title Prepare vaccination data
##'
##' @inheritParams spim_pars
##'
##' @param uptake A vector of proportion update by age. Must be length 19
##'
##' @param end_date The final date of the simulation; we will expand
##'   out the vaccination schedule to this date (if the simulation
##'   continues beyond that date, it's fine)
##'
##' @param mean_days_between_doses The number of days between doses
##'
##' @param efficacy Vaccine efficacy data. This will be a list with
##'   elements `rel_susceptibility`, `rel_p_sympt`,
##'   `rel_p_hosp_if_sympt`, `rel_p_death`, `rel_infectivity` (relative effect
##'   of susceptibility, symptomatic disease, hospitalisation conditional
##'   on symptomatic, death conditional on hospitalisation, and infectivity).
##'   Each of the list elements must
##'   be a `19 x n` matrix (where `n` is the number of vaccination
##'   classes and is shared across all elements). Every number in the
##'   matrix must be a proportion and decrease across columns.
##'
##' @param data Vaccination data (TODO: DESCRIBE CONTENTS)
##'
##' @param booster_proportion Proportion of the groups in
##'  `priority_population` to boost, default is all groups; ignored if
##'  `booster_daily_doses_value` is NULL.
##'
##' @return A list suitable for passing to `spim_pars` as
##'   `vaccination`, containing the new vaccination schedule
##'   (important bits are `efficacy` and `schedule`, but other
##'   elements may be used elsewhere)
##'
##' @export
spim_vaccination_data <- function(date, region, uptake, end_date,
                                  mean_days_between_doses, efficacy,
                                  data, booster_proportion = rep(1L, 19)) {

  browser()

  data$region <- tolower(data$region)
  data$region <- gsub(" ", "_", data$region)

  ## We have a 16-19 age band for Scotland
  data$age_band_min[data$age_band_min == 16] <- 15

  ## Summing over vaccine types
  data <- data %>%
    dplyr::group_by(.data$date, .data$region,
                    .data$age_band_min, .data$age_band_max) %>%
    dplyr::summarise(
      dose1 = sum(.data$first_dose, na.rm = TRUE),
      dose2 = sum(.data$second_dose, na.rm = TRUE),
      dose3 = sum(.data$third_dose, na.rm = TRUE))

  ## This will ensure all regions have data over the same date
  ## and age band combinations (including age_band_min = NA -
  ## which corresponds to any non-age specific doses)
  data <- data %>%
    tidyr::pivot_longer(cols = starts_with("dose"),
                        names_to = "dose",
                        values_to = "count") %>%
    tidyr::pivot_wider(names_from = region,
                       values_from = count,
                       values_fill = 0) %>%
    tidyr::pivot_longer(cols = !c(date, age_band_min, age_band_max, dose),
                        names_to = "region",
                        values_to = "count")

  ## Filter to region and make sure all dates include all age bands
  data <- data[data$region == region, ]
  data <- data %>%
    tidyr::pivot_wider(names_from = date,
                       values_from = count,
                       values_fill = 0) %>%
    tidyr::pivot_longer(cols = !c(age_band_min, age_band_max, dose, region),
                        names_to = "date",
                        values_to = "count") %>%
    tidyr::pivot_wider(names_from = dose,
                       values_from = count)

  ## Order columns and rows nicely
  data$date <- as.Date(data$date)
  data <- data[, c("date", "region", "age_band_min", "age_band_max",
                   "dose1", "dose2", "dose3")]
  data <- dplyr::arrange(data, date, age_band_min)

  ## DONE UP TO HERE

  priority_population <- sircovid::vaccine_priority_population(region, uptake)

  if (region == "northern_ireland") {
    ## No age disaggregated data, so we have to do it ourselves
    if (!all(is.na(data$age_band_min))) {
      stop("Found unexpected age-disaggregated data for northern_ireland")
    }

    data$first_dose[is.na(data$first_dose)] <- 0
    data$second_dose[is.na(data$second_dose)] <- 0
    data$doses <- data$first_dose + data$second_dose
    ## Trim trailing zeros
    data <- data[data$date <= max(data$date[data$doses > 0]), ]

    ## add in missing values - name only for now so that we can check each case
    ## individually
    missing_date <- as.Date("2020-12-06")
    if (!(missing_date %in% data$date) && ((missing_date + 1) %in% data$date)) {
      missing_data <- data[1, ]
      missing_data$date <- missing_date
      data <- rbind(data, missing_data)
      data <- data[order(data$date), ]
    }

    # There are two repeated dates in NI data, let's remove them
    data <- data[data$date >= missing_date, ]
    if (any(diff(data$date) != 1)) {
      stop("non-consecutive dates found in northern_ireland vaccination data")
    }

    date_start <- sircovid::sircovid_date(data$date[[1]])

    doses <- c(
      data$doses,
      rep(mean(tail(data$doses, 7)), end_date - (date_start + nrow(data))))
      ## no boosters data for NI currently
    boosters <- 0

    schedule <- sircovid::vaccine_schedule_future(
      date_start, doses, mean_days_between_doses,
      priority_population, booster_daily_doses_value = boosters,
      booster_proportion = booster_proportion)
  } else {
    # A number of vaccines have unexpectedly been allocoated to an NA age-group
    # in Scotland, let's filter them out for the time being
    nations <- c("scotland", "wales")
    if (region %in% nations) {
      data <- data %>% dplyr::filter(!is.na(.data$age_band_min))
    }
    i <- is.na(data$age_band_min)
    f <- function(x) {
      !is.na(x) & x > 0
    }
    err <- f(data$first_dose[i]) | f(data$second_dose[i])
    if (any(err)) {
      stop("Found unexpected non-age-disaggregated data for ", region)
    }

    ## Remove any data after the date parameter
    data <- data[data$date <= as.Date(date) & !is.na(data$age_band_min), ]

    if (region == "scotland") {
      ## no boosters data for Scotland currently
      boosters <- 0
      booster_start_date <- NULL
    } else {
      ## booster data is age specific, but we do not have support for
      ## this in sircovid yet
      data_booster <- data %>%
        dplyr::group_by(.data$date) %>%
        dplyr::summarise(
          dose3 = sum(.data$third_dose, na.rm = TRUE))
      booster_start_date <-
        data_booster$date[min(which(data_booster$dose3 > 0))]
      boosters <- data_booster$dose3[data_booster$date >= booster_start_date]
    }

    data <- data %>%
      dplyr::group_by(.data$date, .data$age_band_min, .data$age_band_max) %>%
      dplyr::summarise(
        dose1 = sum(.data$first_dose, na.rm = TRUE),
        dose2 = sum(.data$second_dose, na.rm = TRUE))

    last_day <- (data %>%
                 dplyr::group_by(.data$date) %>%
                 dplyr::summarise(
                   doses = sum(.data$dose1 + .data$dose2, na.rm = TRUE)) %>%
                 dplyr::filter(.data$doses > 0) %>%
                 tail(1))$date

    data <- data[data$date <= last_day, ]

    schedule <- sircovid::vaccine_schedule_data_future(data, region, uptake,
                                                       end_date,
                                                       mean_days_between_doses,
                                                       boosters,
                                                       booster_proportion)


    if (!is.null(booster_start_date)) {
      ## This is a bit of a hack - currently booster doses are added on after
      ## the last day in the data. So we shift them to the correct day here.
      ## Once the vaccine schedule code is refactored to deal with full
      ## booster data, we can remove this
      n_days <- sircovid::sircovid_date(end_date) - schedule$date + 1
      n_days_data <- sircovid::sircovid_date(last_day) - schedule$date + 1
      schedule_boosters <-
        schedule$doses[, 3, seq(n_days_data + 1, dim(schedule$doses)[3])]
      booster_dates <-
        seq(sircovid::sircovid_date(booster_start_date) -
              schedule$date + 1, by = 1, length.out = length(boosters))
      schedule$doses[, 3, booster_dates] <- schedule_boosters
      schedule$doses <- schedule$doses[, , seq_len(n_days)]
    }

  }

  i <- seq_len(sircovid::sircovid_date(date) - schedule$date + 1)
  schedule_real <- schedule
  schedule_real$doses <- schedule$doses[, , i, drop = FALSE]

  ret <- list(
    data = data,
    mean_days_between_doses = mean_days_between_doses,
    priority_population = priority_population,
    schedule = schedule,
    schedule_real = schedule_real,
    efficacy = efficacy)
  class(ret) <- "spim_vaccination_data" # soon
  ret
}
