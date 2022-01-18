##' Prepare vaccination data
##'
##' @title Prepare vaccination data
##'
##' @inheritParams spim_data
##'
##' @param uptake A matrix of proportion uptake by age, with (i, j) th
##'   value being the uptake of dose j for group i. Must have 19 rows.
##'
##' @param days_to_effect A vector of the number of days each dose is
##'   assumed to take to reach effect
##'
##' @param data Vaccination data (TODO: DESCRIBE CONTENTS)
##'
##' @return A list suitable for passing to `spim_pars` as
##'   `vaccination`, containing the new vaccination schedule
##'   (important bits are `efficacy` and `schedule`, but other
##'   elements may be used elsewhere)
##'
##' @export
spim_vaccination_data <- function(date, region, uptake, days_to_effect, data) {

  ## Vaccination started 2020-12-08. There should be no doses before this
  ## so we remove dates before this
  data <- data[data$date >= "2020-12-08", ]

  ## Boosters started 2021-09-15, so ignore any third doses before this date
  data$third_dose[data$date < "2021-09-15"] <- 0

  ## Remove any data after date parameter
  data <- data[data$date <= date, ]

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

  ## Fill in any missing dates
  missing_dates <-
    setdiff(as.character(seq.Date(as.Date(min(data$date)), date, by = 1)),
                         unique(data$date))
  if (length(missing_dates) > 0) {
    data_missing <- data.frame(date = missing_dates,
                               age_band_min = NA,
                               age_band_max = NA,
                               dose1 = 0,
                               dose2 = 0,
                               dose3 = 0)
    data <- rbind(data, data_missing)
  }
  data <- data %>%
    dplyr::arrange(date, .data$age_band_min)

  ## we remove any trailing days with zero doses if they are within
  ## min(days_to_effect) of the date parameter
  agg_data <- data %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(
      doses = sum(.data$dose1) + sum(.data$dose2) + sum(.data$dose3))
  dates_remove <-
    agg_data$date[agg_data$date > max(agg_data$date[agg_data$doses > 0])]
  dates_remove <-
    dates_remove[dates_remove > date - min(days_to_effect)]
  data <- data[!(data$date %in% dates_remove), ]

  schedule <- sircovid::vaccine_schedule_from_data(data, region, uptake)


  ret <- list(
    data = data,
    schedule = schedule)
  class(ret) <- "spim_vaccination_data" # soon
  ret
}
