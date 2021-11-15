##' Prepare vaccination data
##'
##' @title Prepare vaccination data
##'
##' @inheritParams spim_pars
##'
##' @param uptake A matrix of proportion uptake by age, with (i, j) th
##'   value being the uptake of dose j for group i. Must have 19 rows.
##'
##' @param days_to_effect A vector of the number of days each dose is
##'   assumed to take to reach effect
##'
##' @param end_date The final date of the simulation; we will expand
##'   out the vaccination schedule to this date (if the simulation
##'   continues beyond that date, it's fine)
##'
##' @param data Vaccination data (TODO: DESCRIBE CONTENTS)
##'
##' @return A list suitable for passing to `spim_pars` as
##'   `vaccination`, containing the new vaccination schedule
##'   (important bits are `efficacy` and `schedule`, but other
##'   elements may be used elsewhere)
##'
##' @export
spim_vaccination_data <- function(date, region, uptake, days_to_effect,
                                  end_date, data) {

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

  data <- data[data$region == region, ]

  ## we remove any trailing days with zero doses if they are within
  ## min(vaccine_days_to_effect) of the end_date
  agg_data <- data %>%
    dplyr::group_by(date) %>%
    summarise(doses = sum(.data$dose1) + sum(.data$dose2) + sum(.data$dose3))
  dates_remove <-
    agg_data$date[agg_data$date > max(agg_data$date[agg_data$doses > 0])]
  dates_remove <-
    dates_remove[dates_remove > end_date - min(days_to_effect)]
  data <- data[!(data$date %in% dates_remove), ]

  schedule <- sircovid::vaccine_schedule_from_data(data, region, uptake)


  ret <- list(
    data = data,
    schedule = schedule)
  class(ret) <- "spim_vaccination_data" # soon
  ret
}
