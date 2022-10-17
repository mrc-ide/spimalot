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
##' @param n_doses Number of doses to include
##'
##' @param dose_start_dates A vector of the start dates for each dose. Any doses
##'   before the corresponding start date will be ignored
##'
##' @return A list suitable for passing to `spim_pars` as
##'   `vaccination`, containing the new vaccination schedule
##'   (important bits are `efficacy` and `schedule`, but other
##'   elements may be used elsewhere)
##'
##' @export
spim_vaccination_data <- function(date, region, uptake, days_to_effect, data,
                                  n_doses, dose_start_dates,
                                  carehomes = TRUE) {

  dose_cols <- paste0("dose", seq_len(n_doses))
  if (!all(dose_cols %in% names(data))) {
    stop(sprintf("n_doses = %s so expected dose column names: %s",
                 n_doses, paste(squote(dose_cols), collapse = ", ")))
  }
  if (length(dose_start_dates) != n_doses) {
    stop(sprintf("n_doses = %s so expected length of dose_start_dates to be %s",
                 n_doses, n_doses))
  }

  ## Remove dose columns above n_doses
  keep <- c("date", "region", "vaccine", "age_band_min", "age_band_max",
            dose_cols)
  data <- data[, keep]

  ## Remove all data before earliest dose start date
  data <- data[data$date >= min(dose_start_dates), ]

  ## Ignore dose data before each corresponding start date
  for (i in seq_len(n_doses)) {
    data[[paste0("dose", i)]][data$date < dose_start_dates[i]] <- 0
  }

  ## Remove any data after date parameter
  data <- data[data$date <= date, ]

  ## We have a 16-19 age band for Scotland
  data$age_band_min[data$age_band_min == 16] <- 15

  data$region <- tolower(data$region)
  data$region <- gsub(" ", "_", data$region)
  data <- data[data$region == region, ]

  ## Summing over vaccine types
  data <- data %>%
    tidyr::pivot_longer(dose_cols, names_to = "dose", values_to = "count") %>%
    dplyr::group_by(.data$date, .data$region,
                    .data$age_band_min, .data$age_band_max, .data$dose) %>%
    dplyr::summarise(count = sum(.data$count, na.rm = TRUE))

  ## Fill in any missing dates
  missing_dates <-
    setdiff(as.character(seq.Date(as.Date(min(data$date)), date, by = 1)),
                         unique(data$date))
  if (length(missing_dates) > 0) {
    data_missing <- data.frame(date = missing_dates,
                               region = region,
                               age_band_min = NA,
                               age_band_max = NA,
                               dose = "dose1",
                               count = 0)
    data <- rbind(data, data_missing)
  }
  data <- data %>%
    dplyr::arrange(date, .data$age_band_min)

  ## we remove any trailing days with zero doses if they are within
  ## min(days_to_effect) of the date parameter
  agg_data <- data %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(count = sum(count))
  dates_remove <-
    agg_data$date[agg_data$date > max(agg_data$date[agg_data$count > 0])]
  dates_remove <-
    dates_remove[dates_remove > date - min(days_to_effect)]
  data <- data[!(data$date %in% dates_remove), ]

  data <- data %>%
    tidyr::pivot_wider(names_from = dose, values_from = count)

  schedule <- sircovid::vaccine_schedule_from_data(data, region, uptake,
                                                   carehomes)


  ret <- list(
    data = data,
    schedule = schedule)
  class(ret) <- "spim_vaccination_data" # soon
  ret
}
