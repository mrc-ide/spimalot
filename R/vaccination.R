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
##' @return A list suitable for passing to `spim_pars` as
##'   `vaccination`, containing the new vaccination schedule
##'   (important bits are `efficacy` and `schedule`, but other
##'   elements may be used elsewhere)
##'
##' @export
spim_vaccination_data <- function(date, region, uptake, end_date, data) {

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

  schedule <- sircovid::vaccine_schedule_from_data(data, region, uptake)

  ret <- list(
    data = data,
    schedule = schedule)
  class(ret) <- "spim_vaccination_data" # soon
  ret
}
