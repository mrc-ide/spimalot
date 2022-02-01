##' Organise data inputs that the model will be fit to. Note that
##' there are other data inputs (notably severity and vaccination
##' data) that are not fit to but are treated as model parameters.
##'
##' @title Organise input data
##'
##' @param date The date of the data. This will be the final date
##'
##' @param model_type Either "BB" (beta-binomial) or "NB" (negative
##'   binomial). We use and exclude different data types depending on
##'   this argument.
##'
##' @param region The region that we are fitting to. This will either
##'   be one of the 10 NHS regions (one of
##'   `sircovid::regions("all")`), or an argument to that function
##'   (e.g., `all`)
##'
##' @param rtm The "RTM collate" data feed. This contains information
##'   about deaths, positive test results, react serology, etc. TODO:
##'   DESCRIBE FORMAT; TODO: DESCRIBE ORIGIN
##'
##' @param serology The additional serology data. TODO: DESCRIBE FORMAT
##'
##' @param trim_deaths The number of days of deaths to trim to avoid
##'   back-fill issues. We typically use a value of 4 days.
##'
##' @param trim_pillar2 The number of days of pillar 2 data to trim to avoid
##'   back-fill issues. We typically use a value of 7 days.
##'
##' @param full_data Not sure yet, we'll find out
##'
##' @return A [data.frame()] TODO: describe columns
##'
##' @export
spim_data <- function(date, region, model_type, rtm, serology,
                      trim_deaths, trim_pillar2, full_data = FALSE) {
  spim_check_model_type(model_type)
  if (length(region) == 1) {
    check_region(region)
    spim_data_single(date, region, model_type, rtm, serology, trim_deaths,
                     trim_pillar2, full_data)
  } else {
    ## TODO: better error message here:
    stopifnot(all(region %in% sircovid::regions("all")))
    data <- lapply(region, function(r)
      cbind(
        region = r,
        spim_data_single(date, r, model_type, rtm, serology, trim_deaths,
                         trim_pillar2, full_data),
        stringsAsFactors = FALSE))
    if (length(unique(lapply(data, "[[", "date"))) != 1) {
      ## If this errors, we just need to compute the union set of
      ## dates, then fill in NA data for the missing entries.  I am
      ## not actually sure how the filter will behave with this
      ## though; https://github.com/mrc-ide/mcstate/issues/182
      stop("Align data dates across regions")
    }
    ret <- dplyr::bind_rows(data)
    ret$region <- factor(ret$region)
    ret
  }
}


spim_data_single <- function(date, region, model_type, rtm, serology,
                             trim_deaths, trim_pillar2, full_data) {
  ## TODO: verify that rtm has consecutive days
  rtm <- spim_lancelot_data_rtm(date, region, model_type, rtm, full_data)
  serology <- spim_data_serology(date, region, serology)

  ## Merge the two datasets on date
  stopifnot(all(serology$date %in% rtm$date))
  i <- match(rtm$date, serology$date)
  serology <- serology[i, ]
  rownames(serology) <- NULL
  data <- cbind(rtm, serology[setdiff(names(serology), "date")])

  ## At this point we'll save our "real" date into date_string and
  ## work with "sircovid" dates which are just integer dates into 2020
  data$date_string <- data$date
  data$date <- sircovid::sircovid_date(data$date)

  ## Set last 'trim_deaths' days with deaths reported to NA, as these
  ## are too likely to be back-filled to be reliable
  deaths_hosp_age <- paste0("deaths_hosp_", c(0, seq(50, 80, 5)),
                            "_", c(seq(49, 79, 5), 120))
  deaths_hosp_age <- gsub("120", "plus", deaths_hosp_age)
  i <- seq(to = nrow(data), length.out = trim_deaths)
  data[i, c("deaths", "deaths_hosp", "deaths_comm",
            "deaths_carehomes",  "deaths_non_hosp",
            deaths_hosp_age)] <- NA

  ## Set last 'trim_pillar2' days with pillar 2 reported to NA, as these
  ## are too likely to be back-filled to be reliable
  ## TODO: this works only for Lancelot - carehomes will be soon removed
  i <- seq(to = nrow(data), length.out = trim_pillar2)
  cols_pillar2 <- grep("pillar2", colnames(data), value = TRUE)
  data[i, cols_pillar2] <- NA_integer_

  data
}


##' @importFrom dplyr %>%
spim_lancelot_data_rtm <- function(date, region, model_type, data, full_data) {

  pillar2_over25_age_bands <- c("25_49", "50_64", "65_79", "80_plus")
  pillar2_age_bands <- c("under15", "15_24", pillar2_over25_age_bands)


  deaths_hosp_age <- paste0("death_", c(0, seq(50, 80, 5)),
                            "_", c(seq(49, 79, 5), 120))
  deaths_hosp_age <- gsub("120", "plus", deaths_hosp_age)

  vars <- c("phe_patients", "phe_occupied_mv_beds",  "icu", "general",
            "admitted", "new", "phe_admissions", "all_admission",
            deaths_hosp_age, "death2", "death3", "death_chr", "death_comm",
            "ons_death_carehome", "ons_death_noncarehome", "react_positive",
            "react_samples",
            # VAM data
            "n_symp_alpha_variant", "n_symp_non_alpha_variant",
            "n_symp_delta_variant", "n_symp_non_delta_variant",
            "n_symp_omicron_variant", "n_symp_non_omicron_variant",
            # Other VOC data
            "n_alpha_pred", "n_non_alpha_pred",
            "n_alpha_variant", "n_non_alpha_variant",
            "n_delta_variant", "n_non_delta_variant",
            "n_omicron_variant", "n_non_omicron_variant", "s_positive_adj1",
            "s_negative_adj1",
            # Pillar 2 positives
            "positives", "positives_over25", "pillar2_positives",
            "pillar2_positives_over25",
            paste0("pillar2_positives_", pillar2_age_bands),
            # Pillar 2 positives symptomatic PCR only
            "pillar2_positives_symp_pcr_only",
            "pillar2_positives_symp_pcr_only_over25",
            paste0("pillar2_positives_symp_pcr_only_", pillar2_age_bands),
            # Pillar 2 positive PRC all (includes LFT+PCR and PCR only)
            "pillar2_positives_pcr_all", "pillar2_positives_pcr_all_over25",
            paste0("pillar2_positives_pcr_all_", pillar2_age_bands),
            # Pillar 2 negatives
            "negatives", "pillar2_negatives", "pillar2_negatives_over25",
            paste0("pillar2_negatives_", pillar2_age_bands),
            # Pillar 2 negative PCR
            "pillar2_negatives_total_pcr_over25", "pillar2_negatives_total_pcr",
            paste0("pillar2_negatives_total_pcr_", pillar2_age_bands))
  data <- data[c("region", "date", vars)]

  ## Remove any data after the date parameter
  data <- data[as.Date(data$date) <= as.Date(date), ]

  ## TODO: de-deplyr this and/or make it function-safe
  ## Make sure the dates for each region match up
  rows_out <- data %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(rows = dplyr::n())
  all_regions <- unique(data$region)

  dates_incomplete <- which(rows_out$rows < length(all_regions))

  if (length(dates_incomplete > 0)) {
    for (d in rows_out$date[dates_incomplete]) {
      missing_regions <- all_regions[!all_regions %in%
                                       data[data$date == d, "region"]]

      tmp <-  data %>% filter(date == d)
      tmp_add <- tmp[seq_along(missing_regions), ]
      tmp_add[3:ncol(tmp_add)] <- NA
      tmp_add$region <- missing_regions
      tmp_add$date <- d
      data <- data %>% dplyr::bind_rows(tmp_add)
    }
  }

  data <- data %>% dplyr::arrange(date, region)

  # Set NA deaths to 0
  data[which(is.na(data$death2)), "death2"] <- 0
  data[which(is.na(data$death3)), "death3"] <- 0
  data[, deaths_hosp_age][is.na(data[, deaths_hosp_age])] <- 0
  data[which(is.na(data$death_chr)), "death_chr"] <- 0
  data[which(is.na(data$death_comm)), "death_comm"] <- 0
  data[which(is.na(data$ons_death_carehome)), "ons_death_carehome"] <- 0
  data[which(is.na(data$ons_death_noncarehome)), "ons_death_noncarehome"] <- 0

  if (region == "uk") {
    ## This might be better done in the upstream task
    nations <- c("scotland", "england", "northern_ireland", "wales")
    sub <- data[data$region %in% nations, ]
    data <- aggregate(sub[vars], sub["date"], sum)
  } else {
    data <- data[data$region == region, ]
  }

  if (region %in% c("northern_ireland", "scotland", "wales", "uk")) {
    data$deaths <- data$death2
    data[, deaths_hosp_age] <- NA_integer_
    data$deaths_hosp <- NA_integer_
    data$deaths_comm <- NA_integer_
    data$deaths_carehomes <- NA_integer_
    data$deaths_non_hosp <- NA_integer_
  } else {
    data$deaths <- NA_integer_
    data$deaths_hosp <- data$death3
    data$deaths_non_hosp <- NA_integer_

    ## due to ONS data being lagged, we will use death linelist data
    ## for recent care home and community deaths
    date_death_change <- as.Date(date) - 45
    data$deaths_carehomes <- dplyr::case_when(
      data$date < date_death_change ~ as.integer(data$ons_death_carehome),
      data$date >= date_death_change ~ as.integer(data$death_chr)
    )
    data$deaths_comm <- dplyr::case_when(
      data$date < date_death_change ~ as.integer(data$ons_death_noncarehome),
      data$date >= date_death_change ~ as.integer(data$death_comm)
    )

  }

  ## Fit to Wildtype/Alpha using sgtf for England, COG for S/W/NI
  if (region %in% c("scotland", "wales", "northern_ireland")) {
    data$strain_non_variant <- data$n_non_alpha_variant
    data$strain_tot <- data$n_alpha_variant + data$n_non_alpha_variant
  } else {
    data$strain_non_variant <- data$n_non_alpha_pred
    data$strain_tot <- data$n_alpha_pred + data$n_non_alpha_pred
  }

  # Only use Wildtype/Alpha data between 2020-09-17 and 2021-03-01
  na_strain_dates <-
    data$date < as.Date("2020-09-17") | data$date > as.Date("2021-03-01")
  data$strain_non_variant[na_strain_dates] <- NA_integer_
  data$strain_tot[na_strain_dates] <- NA_integer_

  # Fit to Alpha/Delta using sgtf data for England, COG data for S/W/NI
  if (region %in% c("scotland", "wales", "northern_ireland")) {
    data$strain_non_variant[data$date >= "2021-03-08"] <-
      data$n_non_delta_variant[data$date >= "2021-03-08"]
    data$strain_tot[data$date >= "2021-03-08"] <-
      data$n_delta_variant[data$date >= "2021-03-08"] +
      data$n_non_delta_variant[data$date >= "2021-03-08"]
  } else {
    data$strain_non_variant[data$date >= "2021-03-08"] <-
      data$n_symp_non_delta_variant[data$date >= "2021-03-08"]
    data$strain_tot[data$date >= "2021-03-08"] <-
      data$n_symp_delta_variant[data$date >= "2021-03-08"] +
      data$n_symp_non_delta_variant[data$date >= "2021-03-08"]
  }

  ## Only fit to Alpha/Delta up to 2021-07-31
  na_strain_dates <- data$date > as.Date("2021-07-31")
  data$strain_non_variant[na_strain_dates] <- NA_integer_
  data$strain_tot[na_strain_dates] <- NA_integer_

  ## Fit to Delta/Omicron using VAM data for England, COG data for S/W/NI
  if (region %in% c("scotland", "wales", "northern_ireland")) {
    data$strain_non_variant[data$date >= "2021-11-20"] <-
      data$n_non_omicron_variant[data$date >= "2021-11-20"]
    data$strain_tot[data$date >= "2021-11-20"] <-
      data$n_omicron_variant[data$date >= "2021-11-20"] +
      data$n_non_omicron_variant[data$date >= "2021-11-20"]
  } else {
    data$strain_non_variant[data$date >= "2021-11-20"] <-
      data$n_symp_non_omicron_variant[data$date >= "2021-11-20"]
    data$strain_tot[data$date >= "2021-11-20"] <-
      data$n_symp_non_omicron_variant[data$date >= "2021-11-20"] +
      data$n_symp_omicron_variant[data$date >= "2021-11-20"]
  }


  # Use positives/negatives as Pillar 2 for Scotland
  # Set data$phe_patients to NA between 2020-06-01 and 2020-09-09 (inclusive)
  if (region == "scotland") {
    data$pillar2_positives <- data$positives
    ## Scotland negatives are by report date (while positives are by specimen
    ## date). We assume a 2 day reporting delay.
    data$pillar2_negatives <- c(data$negatives[-c(1, 2)], rep(NA_integer_, 2))
    data$pillar2_positives_over25 <- data$positives_over25
    ## We do not have any age breakdown for negatives for Scotland
    data$pillar2_negatives_over25 <- NA_integer_
    for (i in pillar2_age_bands) {
      data[, paste0("pillar2_positives_", i)] <- NA_integer_
      data[, paste0("pillar2_negatives_", i)] <- NA_integer_
    }

    data$phe_patients[data$date >= as.Date("2020-06-01") &
                        data$date <= as.Date("2020-09-10")] <- NA_integer_
  }

  # Use positives/negatives as Pillar 2 for NI
  if (region == "northern_ireland") {
    data$pillar2_positives <- data$positives
    data$pillar2_negatives <- data$negatives
  }

  data$pillar2_cases <- data$pillar2_positives
  data$pillar2_cases_over25 <- data$pillar2_positives_over25
  data[, paste0("pillar2_cases_", pillar2_age_bands)] <-
    data[, paste0("pillar2_positives_", pillar2_age_bands)]

  ## Use symp PCR only for cases where available
  if (!all(is.na(data$pillar2_positives_symp_pcr_only))) {
    data$pillar2_cases <- data$pillar2_positives_symp_pcr_only
  }
  if (!all(is.na(data$pillar2_positives_symp_pcr_only_over25))) {
    data$pillar2_cases_over25 <- data$pillar2_positives_symp_pcr_only_over25
  }

  ## Use symp PCR only for cases by age where available
  pillar2_symp_PCR_only_by_age <-
    data[, paste0("pillar2_positives_symp_pcr_only_", pillar2_age_bands)]
  if (!all(is.na(pillar2_symp_PCR_only_by_age))) {
    if (!full_data) {
      data$pillar2_cases_over25 <- NA_integer_
    }
    data[, paste0("pillar2_cases_", pillar2_age_bands)] <-
      data[, paste0("pillar2_positives_symp_pcr_only_", pillar2_age_bands)]
  }

  ## Use PCR all for positives where available
  if (!all(is.na(data$pillar2_positives_pcr_all))) {
    data$pillar2_positives <- data$pillar2_positives_pcr_all
  }
  if (!all(is.na(data$pillar2_positives_pcr_all_over25))) {
    data$pillar2_positives_over25 <- data$pillar2_positives_pcr_all_over25
  }

  ## Use PCR all for positives by age where available
  pillar2_positives_pcr_all_by_age <-
    data[, paste0("pillar2_positives_pcr_all_", pillar2_age_bands)]
  if (!all(is.na(pillar2_positives_pcr_all_by_age))) {
    if (!full_data) {
      data$pillar2_positives_over25 <- NA_integer_
    }
    data[, paste0("pillar2_positives_", pillar2_age_bands)] <-
      data[, paste0("pillar2_positives_pcr_all_", pillar2_age_bands)]

  }

  ## Use total PCR for negatives where available
  if (!all(is.na(data$pillar2_negatives_total_pcr))) {
    data$pillar2_negatives <- data$pillar2_negatives_total_pcr
  }
  if (!all(is.na(data$pillar2_negatives_total_pcr_over25))) {
    data$pillar2_negatives_over25 <- data$pillar2_negatives_total_pcr_over25
  }

  ## Use total PCR for negatives by age where available
  pillar2_negatives_total_pcr_by_age <-
    data[, paste0("pillar2_negatives_total_pcr_", pillar2_age_bands)]
  if (!all(is.na(pillar2_negatives_total_pcr_by_age))) {
    if (!full_data) {
      data$pillar2_negatives_over25 <- NA_integer_
    }
    data[, paste0("pillar2_negatives_", pillar2_age_bands)] <-
      data[, paste0("pillar2_negatives_total_pcr_", pillar2_age_bands)]
  }

  # Use hospital data from dashboard for all except Wales (linelist)
  if (region == "wales") {
    data$final_admissions <- data$all_admission
    data$final_icu <- data$icu
    data$final_general <- data$general
    data$final_hosp <- data$icu + data$general

  } else {
    data$final_admissions <- data$phe_admissions
    data$final_icu <- data$phe_occupied_mv_beds
    data$final_general <- data$phe_patients - data$phe_occupied_mv_beds
    data$final_hosp <- data$phe_patients
  }

  cols_pillar2 <- c("pillar2_positives", "pillar2_negatives", "pillar2_cases",
                    paste0("pillar2_positives_",
                           c("over25", pillar2_age_bands)),
                    paste0("pillar2_negatives_",
                           c("over25", pillar2_age_bands)),
                    paste0("pillar2_cases_",
                           c("over25", pillar2_age_bands)))

  # Turn NAs to zeroes for pillar 2 columns where data is available
  for (i in cols_pillar2) {
    if (!all(is.na(data[, i]))) {
      data[which(is.na(data[, i])), i] <- 0
    }
  }

  # ignore pillar 2 testing before 2020-06-18
  data[which(data$date < "2020-06-18"), cols_pillar2] <- NA_integer_

  last_week <- seq(to = nrow(data), length.out = 7)
  ## Remove last week admissions for Wales (due to backfill)
  if (region == "wales") {
    data[last_week, "final_admissions"] <- NA_integer_
  }

  ## Remove implausible value for MV beds occupancy in east_of_england
  ## on 2020-09-11
  if (region == "east_of_england") {
    data[which(data$final_general < 0), "final_general"] <- NA_integer_
  }

  ## Remove implausible values for pillar2_negatives data
  pillar2_negatives_cols <-
    c("pillar2_negatives", paste0("pillar2_negatives_",
                                  c("over25", pillar2_age_bands)))

  for (i in pillar2_negatives_cols) {
    data[which(data[, i] < 0), i] <- NA_integer_
  }

  ## If we do not have negatives, set corresponding positives to 0
  for (i in c(paste0("_", pillar2_age_bands), "_over25", "")) {
    if (all(is.na(data[, paste0("pillar2_negatives", i)]))) {
      data[, paste0("pillar2_positives", i)] <- NA_integer_
    }
  }

  ## Check all pillar 2 data is greater than 0 or NA
  data_pillar2 <- data[, cols_pillar2]
  stopifnot(all(data_pillar2 >= 0, na.rm = TRUE))

  ## TODO: with a stripped down compare function wee could drop the NA
  ## columns here.
  ret <- data_frame(
    date = sircovid::as_date(data$date),
    deaths_hosp = data$deaths_hosp,
    deaths_hosp_0_49 = data$death_0_49,
    deaths_hosp_50_54 = data$death_50_54,
    deaths_hosp_55_59 = data$death_55_59,
    deaths_hosp_60_64 = data$death_60_64,
    deaths_hosp_65_69 = data$death_65_69,
    deaths_hosp_70_74 = data$death_70_74,
    deaths_hosp_75_79 = data$death_75_79,
    deaths_hosp_80_plus = data$death_80_plus,
    deaths_comm = data$deaths_comm + data$deaths_carehomes,
    deaths_carehomes = NA_integer_,
    deaths_non_hosp = data$deaths_non_hosp,
    icu = data$final_icu,
    general = data$final_general,
    hosp = data$final_hosp,
    deaths = data$deaths,
    admitted = data$admitted,
    diagnoses = data$new,
    all_admission = data$final_admissions,
    pillar2_tot = data$pillar2_positives + data$pillar2_negatives,
    pillar2_pos = data$pillar2_positives,
    pillar2_cases = data$pillar2_cases,
    pillar2_over25_tot = data$pillar2_positives_over25 +
      data$pillar2_negatives_over25,
    pillar2_over25_pos = data$pillar2_positives_over25,
    pillar2_over25_cases = data$pillar2_cases_over25,
    pillar2_under15_tot = data$pillar2_positives_under15 +
      data$pillar2_negatives_under15,
    pillar2_under15_pos = data$pillar2_positives_under15,
    pillar2_under15_cases = data$pillar2_cases_under15,
    pillar2_15_24_tot = data$pillar2_positives_15_24 +
      data$pillar2_negatives_15_24,
    pillar2_15_24_pos = data$pillar2_positives_15_24,
    pillar2_15_24_cases = data$pillar2_cases_15_24,
    pillar2_25_49_tot = data$pillar2_positives_25_49 +
      data$pillar2_negatives_25_49,
    pillar2_25_49_pos = data$pillar2_positives_25_49,
    pillar2_25_49_cases = data$pillar2_cases_25_49,
    pillar2_50_64_tot = data$pillar2_positives_50_64 +
      data$pillar2_negatives_50_64,
    pillar2_50_64_pos = data$pillar2_positives_50_64,
    pillar2_50_64_cases = data$pillar2_cases_50_64,
    pillar2_65_79_tot = data$pillar2_positives_65_79 +
      data$pillar2_negatives_65_79,
    pillar2_65_79_pos = data$pillar2_positives_65_79,
    pillar2_65_79_cases = data$pillar2_cases_65_79,
    pillar2_80_plus_tot = data$pillar2_positives_80_plus +
      data$pillar2_negatives_80_plus,
    pillar2_80_plus_pos = data$pillar2_positives_80_plus,
    pillar2_80_plus_cases = data$pillar2_cases_80_plus,
    react_pos = data$react_positive,
    react_tot = data$react_samples,
    strain_non_variant = data$strain_non_variant,
    strain_tot = data$strain_tot,
    strain_over25_non_variant = NA_integer_,
    strain_over25_tot = NA_integer_)

  if (!full_data) {
    ## Typically we do not fit to this
    ret$strain_over25_non_variant <- NA_integer_
    ret$strain_over25_tot <- NA_integer_

    ## Do not fit to under 15 pillar 2
    ret$pillar2_under15_tot <- NA_integer_
    ret$pillar2_under15_pos <- NA_integer_
    ret$pillar2_under15_cases <- NA_integer_

    ## Do not fit to aggregated hospital deaths
    deaths_hosp_age <- gsub("death", "deaths_hosp", deaths_hosp_age)
    if (any(!is.na(ret[, deaths_hosp_age]))) {
      ret$deaths_hosp <- NA_integer_
    }

    if (model_type == "BB") {
      omit <- c("hosp", "admitted", "diagnoses", "pillar2_cases",
                "pillar2_over25_cases",
                paste0("pillar2_", pillar2_age_bands, "_cases"))

      for (i in omit) {
        ret[[i]] <- NA_integer_
      }

      ## If we fit pillar 2 to any over 25 sub-age bands, do not fit to
      ## aggregated over 25
      fit_to_pillar2_age_bands_over25 <-
        !all(is.na(ret[, paste0("pillar2_", pillar2_over25_age_bands, "_tot")]))
      if (fit_to_pillar2_age_bands_over25) {
        ret$pillar2_over25_tot <- NA_integer_
        ret$pillar2_over25_pos <- NA_integer_
      }

      ## If we fit pillar 2 to any age bands (including over 25), do not fit to
      ## all ages aggregated
      fit_to_pillar2_age_bands <-
        !all(is.na(ret[, paste0("pillar2_", c("over25", pillar2_age_bands),
                                "_tot")]))
      if (fit_to_pillar2_age_bands) {
        ret$pillar2_tot <- NA_integer_
        ret$pillar2_pos <- NA_integer_
      }

    }
    if (model_type == "NB") {
      omit <- c("hosp", "admitted", "diagnoses", "pillar2_tot", "pillar2_pos",
                "pillar2_over25_tot", "pillar2_over25_pos",
                paste0("pillar2_", pillar2_age_bands, "_tot"),
                paste0("pillar2_", pillar2_age_bands, "_pos"))

      for (i in omit) {
        ret[[i]] <- NA_integer_
      }

      ## If we fit pillar 2 to any over 25 sub-age bands, do not fit to
      ## aggregated over 25
      fit_to_pillar2_age_bands_over25 <-
        !all(is.na(ret[, paste0("pillar2_", pillar2_over25_age_bands,
                                "_cases")]))
      if (fit_to_pillar2_age_bands_over25) {
        ret$pillar2_over25_cases <- NA_integer_
      }

      ## If we fit pillar 2 to any age bands (including over 25), do not fit to
      ## all ages aggregated
      fit_to_pillar2_age_bands <-
        !all(is.na(ret[, paste0("pillar2_", c("over25", pillar2_age_bands),
                                "_cases")]))
      if (fit_to_pillar2_age_bands) {
        ret$pillar2_cases <- NA_integer_
      }

    }
  }

  stopifnot(
    all(ret$deaths_carehomes >= 0, na.rm = TRUE),
    all(ret$general >= 0, na.rm = TRUE))

  ret
}


##' @importFrom dplyr .data
spim_data_serology <- function(date, region, data) {
  ## For R CMD check's scoping check
  assay <- n_positive <- total_samples <- NULL

  ## We might have serology data that is too recent; subset it here:

  if (region == "scotland") {
    data <- data %>%
      dplyr::mutate(age_group = ifelse(.data$age_group %in% c("0_39", "40_59"),
                                       "15_64", .data$age_group)) %>%
      dplyr::mutate(assay = ifelse(.data$assay == "abbott",
                                   "euro_immun", assay)) %>%
      dplyr::group_by(region, date, .data$age_group, assay) %>%
      dplyr::summarise(n_positive = sum(.data$n_positive),
                       total_samples = sum(.data$total_samples))

  }

  ## TODO: needs effort to fix
  data <- data %>%
    dplyr::filter(.data$assay %in% c("euro_immun", "roche_n")) %>%
    tidyr::pivot_wider(names_from = assay,
                       values_from = c(n_positive, total_samples),
                       values_fill = NA_integer_)

  data <- data[data$region == region & data$age_group == "15_64", ]

  ## Remove any data after the date parameter
  data <- data[as.Date(data$date) <= as.Date(date), ]

  ## Set EuroImmun data to NA after date_remove
  date_remove <- "2021-01-15"
  euro_immun <- c("n_positive_euro_immun", "total_samples_euro_immun")
  data[as.Date(data$date) >= as.Date(date_remove), euro_immun] <- NA_integer_

  data_frame(date = sircovid::as_date(data$date),
             sero_pos_15_64_1 = data$n_positive_euro_immun,
             sero_tot_15_64_1 = data$total_samples_euro_immun,
             sero_pos_15_64_2 = data$n_positive_roche_n,
             sero_tot_15_64_2 = data$total_samples_roche_n)
}
