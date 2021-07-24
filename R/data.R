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
##' @param serology The additional serology data. TODO: DESCRBE FORMAT
##'
##' @param trim_deaths The number of days of deaths to trim to avoid
##'   back-fill issues. We typically use a value of 4 days.
##'
##' @param trim_pillar2 The number of days of pillar 2 data to trim to avoid
##'   back-fill issues. We typically use a value of 7 days.
##'
##' @param full_data Not sure yet, we'll find out
##'
##' @param fit_to_variants Logical, whether to fit to variants data or not
##'
##' @return A [data.frame()] TODO: describe columns
##'
##' @export
spim_data <- function(date, region, model_type, rtm, serology,
                      trim_deaths, trim_pillar2, full_data = FALSE,
                      fit_to_variants = FALSE) {
  check_region(region)
  spim_check_model_type(model_type)

  if (length(region) > 1) {
    ## See the original task
    stop("Not yet supported")
  } else {
    spim_data_single(date, region, model_type, rtm, serology,
                     trim_deaths, trim_pillar2, full_data, fit_to_variants)
  }
}


spim_data_single <- function(date, region, model_type, rtm, serology,
                             trim_deaths, trim_pillar2, full_data,
                             fit_to_variants) {
  ## TODO: verify that rtm has consecutive days
  rtm <- spim_data_rtm(date, region, model_type, rtm, full_data,
                       fit_to_variants)
  serology <- spim_data_serology(date, region, serology)

  ## Merge the two datasets on date
  stopifnot(all(serology$date %in% rtm$date))
  i <- match(rtm$date, serology$date)
  serology <- serology[i, ]
  rownames(serology) <- NULL
  data <- cbind(rtm, serology[setdiff(names(serology), "date")])

  ## We need a dummy row here so that the sampler stops before the
  ## first day in order to get the differences in cumulative deaths
  ## correct.
  ## TODO: This is no longer needed!
  data <- rbind(data[1, ], data)
  data[1, ] <- NA
  data$date[[1]] <- data$date[[2]] - 1

  ## At this point we'll save our "real" date into date_string and
  ## work with "sircovid" dates which are just integer dates into 2020
  data$date_string <- data$date
  data$date <- sircovid::sircovid_date(data$date)

  ## Set last 'trim_deaths' days with deaths reported to NA, as these
  ## are too likely to be back-filled to be reliable
  i <- seq(to = nrow(data), length.out = trim_deaths)
  data[i, c("deaths", "deaths_hosp", "deaths_comm",
            "deaths_carehomes",  "deaths_non_hosp")] <- NA

  ## Set last 'trim_pillar2' days with pillar 2 reported to NA, as these
  ## are too likely to be back-filled to be reliable
  i <- seq(to = nrow(data), length.out = trim_pillar2)
  cols_pillar2 <- c("pillar2_tot", "pillar2_pos", "pillar2_cases",
                    "pillar2_over25_tot", "pillar2_over25_pos",
                    "pillar2_over25_cases")
  data[i, cols_pillar2] <- NA_integer_

  data
}


##' @importFrom dplyr %>%
spim_data_rtm <- function(date, region, model_type, data, full_data,
                          fit_to_variants) {

  vars <- c("phe_patients", "phe_occupied_mv_beds",  "icu", "general",
            "admitted", "new", "phe_admissions", "all_admission",
            "death2", "death3", "death_chr", "death_comm",
            "ons_death_carehome", "ons_death_noncarehome",
            "pillar2_positives", "pillar2_negatives",
            "positives", "negatives", "react_positive", "react_samples",
            "pillar2_negatives_total_pcr_over25", "pillar2_negatives_total_pcr",
            "pillar2_positives_over25", "pillar2_negatives_over25",
            "positives_over25", "pillar2_positives_symp_pcr_only",
            "pillar2_positives_symp_pcr_only_over25",
            "pillar2_positives_pcr_all", "pillar2_positives_pcr_all_over25",
            "n_delta_variant", "n_non_delta_variant",
            "n_symp_delta_variant", "n_symp_non_delta_variant")
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

  # Set NA deaths to 0
  data[which(is.na(data$death2)), "death2"] <- 0
  data[which(is.na(data$death3)), "death3"] <- 0
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
    data$deaths_hosp <- NA_integer_
    data$deaths_comm <- NA_integer_
    data$deaths_carehomes <- NA_integer_
    data$deaths_non_hosp <- NA_integer_
  } else {
    data$deaths <- NA_integer_
    data$deaths_hosp <- data$death3
    data$deaths_non_hosp <- NA_integer_

    if (!full_data) {
      date_death_change <- "2020-09-01"
      data$deaths_carehomes <- dplyr::case_when(
        data$date < date_death_change ~ as.integer(data$ons_death_carehome),
        data$date >= date_death_change ~ NA_integer_
      )
      data$deaths_comm <- dplyr::case_when(
        data$date < date_death_change ~ as.integer(data$ons_death_noncarehome),
        data$date >= date_death_change ~ NA_integer_
      )
    } else {
      ## due to ONS data being lagged, in the full_data version (not used in
      ## fitting) we will use death linelist data for recent care home and
      ## community deaths
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
  }

  # Use VAM data available
  if (region %in% c("scotland", "wales", "northern_ireland")){
    data$strain_non_variant <- data$n_non_delta_variant
    data$strain_tot <- data$n_delta_variant + data$n_non_delta_variant
  } else {
    data$strain_non_variant <- data$n_symp_non_delta_variant
    data$strain_tot <- data$n_symp_delta_variant + data$n_symp_non_delta_variant
  }

  # Use positives/negatives as Pillar 2 for Scotland
  # Set data$phe_patients to NA between 2020-06-01 and 2020-09-09 (inclusive)
  if (region == "scotland") {
    data$pillar2_positives <- data$positives
    data$pillar2_negatives <- data$negatives
    data$pillar2_positives_over25 <- data$positives_over25
    data$pillar2_negatives_over25 <- data$negatives_over25

    data$phe_patients[data$date >= as.Date("2020-06-01") &
                      data$date <= as.Date("2020-09-10")] <- NA_integer_
  }

  data$pillar2_cases <- data$pillar2_positives
  data$pillar2_cases_over25 <- data$pillar2_positives_over25

  ## Use symp PCR only for cases where available
  if (!all(is.na(data$pillar2_positives_symp_pcr_only))) {
    data$pillar2_cases <- data$pillar2_positives_symp_pcr_only
  }
  if (!all(is.na(data$pillar2_positives_symp_pcr_only_over25))) {
    data$pillar2_cases_over25 <- data$pillar2_positives_symp_pcr_only_over25
  }

  ## Use PCR all for positives where available
  if (!all(is.na(data$pillar2_positives_pcr_all))) {
    data$pillar2_positives <- data$pillar2_positives_pcr_all
  }
  if (!all(is.na(data$pillar2_positives_pcr_all_over25))) {
    data$pillar2_positives_over25 <- data$pillar2_positives_pcr_all_over25
  }

  ## Use total PCR for negatives where available
  if (!all(is.na(data$pillar2_negatives_total_pcr))) {
    data$pillar2_negatives <- data$pillar2_negatives_total_pcr
  }
  if (!all(is.na(data$pillar2_negatives_total_pcr_over25))) {
    data$pillar2_negatives_over25 <- data$pillar2_negatives_total_pcr_over25
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
                    "pillar2_positives_over25", "pillar2_negatives_over25",
                    "pillar2_cases_over25")

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
  data[which(data$pillar2_negatives < 0), "pillar2_negatives"] <- NA_integer_
  data[which(data$pillar2_negatives_over25 < 0), "pillar2_negatives_over25"] <-
    NA_integer_

  stopifnot(
    all(data$pillar2_negatives >= 0, na.rm = TRUE),
    all(data$pillar2_positives >= 0, na.rm = TRUE),
    all(data$pillar2_cases >= 0, na.rm = TRUE),
    all(data$pillar2_negatives_over25 >= 0, na.rm = TRUE),
    all(data$pillar2_positives_over25 >= 0, na.rm = TRUE),
    all(data$pillar2_cases_over25 >= 0, na.rm = TRUE))

  ## TODO: with a stripped down compare function wee could drop the NA
  ## columns here.
  ret <- data_frame(
    date = sircovid::as_date(data$date),
    deaths_hosp = data$deaths_hosp,
    deaths_comm = data$deaths_comm,
    deaths_carehomes = data$deaths_carehomes,
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
    react_pos = data$react_positive,
    react_tot = data$react_samples,
    strain_non_variant = data$strain_non_variant,
    strain_tot = data$strain_tot,
    strain_over25_non_variant = NA_integer_,
    strain_over25_tot = NA_integer_)

  if (!fit_to_variants) {
    ret$strain_non_variant <- NA_integer_
    ret$strain_tot <- NA_integer_
    ret$strain_over25_non_variant <- NA_integer_
    ret$strain_over25_tot <- NA_integer_
  }

  if (!full_data) {
    ## Typically we do not fit to this
    ret$strain_over25_non_variant <- NA_integer_
    ret$strain_over25_tot <- NA_integer_

    if (model_type == "BB") {
      omit <- c("hosp", "admitted", "diagnoses", "pillar2_tot", "pillar2_pos",
                "pillar2_cases", "pillar2_over25_cases")
      for (i in omit) {
        ret[[i]] <- NA_integer_
      }
      if (all(is.na(ret$pillar2_over25_tot))) {
        ret$pillar2_tot <- data$pillar2_positives + data$pillar2_negatives
        ret$pillar2_pos <- data$pillar2_positives
        ret$pillar2_over25_tot <- NA_integer_
        ret$pillar2_over25_pos <- NA_integer_
      }
    }
    if (model_type == "NB") {
      omit <- c("hosp", "admitted", "diagnoses", "pillar2_tot", "pillar2_pos",
                "pillar2_cases", "pillar2_over25_tot", "pillar2_over25_pos")
      for (i in omit) {
        ret[[i]] <- NA_integer_
      }
      if (all(is.na(ret$pillar2_over25_cases))) {
        ret$pillar2_cases <- data$pillar2_positives
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


##' Load admissions data
##'
##' @title Load admissions data
##' @param admissions A data.frame from the admissions by age sitrep
##'
##' @param region Name of the region
##' @export
##' @importFrom dplyr .data
spim_data_admissions <- function(admissions, region) {
  nations <- c("scotland", "wales", "northern_ireland")

  if (region %in% nations) {
    admissions <- data.frame(
      date = unique(admissions$date),
      region = region,
      adm_0 = NA_integer_,
      adm_25 = NA_integer_,
      adm_55 = NA_integer_,
      adm_65 = NA_integer_,
      adm_75 = NA_integer_)
  } else {
    vector_age_bands <- unique(admissions$age_from)
    age_bands <- c("date", "region", "adm_0", "adm_25", "adm_55", "adm_65",
                   "adm_75")

    admissions$region <- gsub(" ", "_", admissions$region)
    admissions <- admissions[admissions$region == region, ]

    admissions <- admissions %>%
      dplyr::group_by(date, .data$age_from) %>%
      dplyr::mutate(age_from = paste0("adm_", .data$age_from)) %>%
      dplyr::mutate(value = sum(admissions)) %>%
      dplyr::slice(1) %>%
      dplyr::select(date, region, .data$age_from, value)

    stopifnot(nrow(admissions) == length(unique(admissions$date)) *
              length(vector_age_bands))

    admissions <- admissions %>%
      tidyr::pivot_wider(id_cols = c(date, region),
                         names_from = .data$age_from)

    admissions$adm_0 <- rowSums(admissions[, 3:4])
    admissions$adm_25 <- rowSums(admissions[, 5:7])
    admissions$adm_75 <- rowSums(admissions[, 11:12])

    admissions <- admissions %>% dplyr::select(dplyr::all_of(age_bands))
  }

  admissions
}
