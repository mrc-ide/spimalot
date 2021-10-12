transform_carehomes <- function(region, model_type, multistrain, beta_date,
                                     vaccination, cross_immunity = NULL,
                                     waning_rate) {
  beta_date <- sircovid::sircovid_date(beta_date)
  assert_is(vaccination, "spim_vaccination_data")

  severity <- read_csv(spimalot_file("extdata/support_severity.csv"))
  progression <- read_csv(spimalot_file("extdata/support_progression.csv"))

  n_strain <- if (multistrain) 2 else 1
  if (is.null(cross_immunity)) {
    cross_immunity <- rep(1, n_strain)
  } else {
    stopifnot(length(cross_immunity) == n_strain)
  }

  function(pars) {
    start_date <- pars[["start_date"]]
    p_G_D <- pars[["p_G_D"]] # death in *community*
    p_G_D_CHR <- pars[["p_G_D_CHR"]] # death in *care home*
    eps <- pars[["eps"]]
    m_CHW <- pars[["m_CHW"]]
    m_CHR <- pars[["m_CHR"]]
    p_ICU <- pars[["p_ICU"]]
    p_ICU_2 <- pars[["p_ICU_2"]]
    p_ICU_D <- pars[["p_ICU_D"]]
    p_H_D <- pars[["p_H_D"]]
    p_W_D <- pars[["p_W_D"]]
    mu_D <- pars[["mu_D"]]
    mu_D_2 <- pars[["mu_D_2"]]
    mu_gamma_H <- pars[["mu_gamma_H"]]
    mu_gamma_H_2 <- pars[["mu_gamma_H_2"]]
    mu_gamma_H_3 <- pars[["mu_gamma_H_3"]]
    p_H <- pars[["p_H"]]
    p_H_2 <- pars[["p_H_2"]]
    p_H_3 <- pars[["p_H_3"]]
    alpha_D <- pars[["alpha_D"]]
    alpha_H <- pars[["alpha_H"]]

    if ("strain_transmission_2" %in% names(pars)) {
      strain_transmission_2 <- pars[["strain_transmission_2"]]
    } else {
      strain_transmission_2 <- NULL
    }

    if ("prop_strain_2" %in% names(pars)) {
      prop_strain_2 <- pars[["prop_strain_2"]]
    } else {
      prop_strain_2 <- 0
    }

    if ("strain_seed_date" %in% names(pars)) {
      strain_seed_date <- pars[["strain_seed_date"]] + seq(0, 6)
      strain_seed_pp <- 2 / 1e6
      regional_pop <- sum(sircovid:::sircovid_population(region))
      strain_seed_rate <- rep(strain_seed_pp * regional_pop, 7)
    } else {
      strain_seed_date <- NULL
      strain_seed_rate <- NULL
    }

    beta_value <- unname(pars[paste0("beta", seq_along(beta_date))])

    if (model_type == "BB") {
      p_NC <- pars[["p_NC"]]
      rho_pillar2_tests <- pars[["rho_pillar2_tests"]]
      ## Total: 40 fitted parameters

      ## Unused in BB fits so these are dummy values
      phi_pillar2_cases <- 0.5
      kappa_pillar2_cases <- 2
    }
    if (model_type == "NB") {
      phi_pillar2_cases <- pars[["phi_pillar2_cases"]]
      kappa_pillar2_cases <- 1 / pars[["alpha_pillar2_cases"]]
      ## Total: 40 fitted parameters

      ## Unused in NB fits so these are dummy values
      p_NC <- 0.002
      rho_pillar2_tests <- 0.01
    }

    if ("p_NC_weekend" %in% names(pars)) {
      p_NC_weekend <- pars[["p_NC_weekend"]]
    } else {
      p_NC_weekend <- p_NC
    }

    if ("phi_pillar2_cases_weekend" %in% names(pars)) {
      phi_pillar2_cases_weekend <- pars[["phi_pillar2_cases_weekend"]]
    } else {
      phi_pillar2_cases_weekend <- phi_pillar2_cases
    }

    ## Set severity parameters based on Bob's analysis and fitted parameters.
    severity[severity$Name == "p_sero_pos_1", 2:18] <- 0.85
    severity[severity$Name == "p_sero_pos_2", 2:18] <- 0.85
    severity[severity$Name == "p_G_D", 2:18] <- 0.05

    ## TODO: make these dates easier to work with?
    p_D_date <- sircovid::sircovid_date(c("2020-04-01", "2020-06-01",
                                          "2020-10-01", "2020-12-15",
                                          "2021-01-15", "2021-02-01"))
    mu_D_vec <- c(1, mu_D, mu_D, mu_D_2, mu_D_2, mu_D) /
      max(c(1, mu_D, mu_D_2))
    p_ICU_D_value <- p_ICU_D * mu_D_vec
    p_H_D_value <- p_H_D * mu_D_vec
    p_W_D_value <- p_W_D * mu_D_vec

    p_ICU_date <- sircovid::sircovid_date(c("2020-04-01", "2020-06-01"))
    p_ICU_value <- c(p_ICU, p_ICU_2)

    p_H_date <- sircovid::sircovid_date(c("2020-10-01", "2020-12-15",
                                          "2021-02-15"))
    p_H_value <- c(p_H, p_H_2, p_H_3)

    p_star_date <- sircovid::sircovid_date(c("2020-03-15", "2020-07-01",
                                             "2020-09-20", "2021-06-27"))
    p_star_value <- c(0.1, 0.42, 0.2, 0.45)

    severity <- sircovid::carehomes_parameters_severity(
      0.25,
      severity,
      p_H = list(value = p_H_value, date = p_H_date),
      p_ICU = list(value = p_ICU_value, date = p_ICU_date),
      p_ICU_D = list(value = p_ICU_D_value, date = p_D_date),
      p_H_D = list(value = p_H_D_value, date = p_D_date),
      p_W_D = list(value = p_W_D_value, date = p_D_date),
      p_G_D = list(value = p_G_D),
      p_G_D_CHR = list(value = p_G_D_CHR),
      p_star = list(value = p_star_value, date = p_star_date))

    ## Set progression parameters based on Bob's analysis and fitted parameters
    k_parameters <-
      progression[grep("^k_", progression$parameter), ]
    gammas <-
      progression[grep("^gamma_", progression$parameter),
                  "value"]
    names(gammas) <-
      progression[grep("^gamma_", progression$parameter),
                  "parameter"]
    gammas <- as.list(gammas)

    # Reduce length of stay; same dates apply
    mu_gamma_H_date <- sircovid::sircovid_date(c("2020-12-01",
                                                 "2021-01-01",
                                                 "2021-03-01",
                                                 "2021-06-15"))
    mu_gamma_H_value <- c(1, 1 / mu_gamma_H, 1 / mu_gamma_H_2, 1 / mu_gamma_H_3)

    gamma_ICU_pre <- gammas$gamma_ICU_pre
    gamma_ICU_D <- gammas$gamma_ICU_D
    gamma_ICU_W_D <- gammas$gamma_ICU_W_D
    gamma_ICU_W_R <- gammas$gamma_ICU_W_R
    gamma_H_R_value <- gammas$gamma_H_R * mu_gamma_H_value
    gamma_H_D_value <- gammas$gamma_H_D * mu_gamma_H_value
    gamma_W_R_value <- gammas$gamma_W_R * mu_gamma_H_value
    gamma_W_D_value <- gammas$gamma_W_D * mu_gamma_H_value

    progression <- sircovid::carehomes_parameters_progression(
      0.25,
      gamma_ICU_pre = list(value = gamma_ICU_pre),
      gamma_H_D = list(value = gamma_H_D_value, date = mu_gamma_H_date),
      gamma_H_R = list(value = gamma_H_R_value, date = mu_gamma_H_date),
      gamma_ICU_D = list(value = gamma_ICU_D),
      gamma_ICU_W_D = list(value = gamma_ICU_W_D),
      gamma_ICU_W_R = list(value = gamma_ICU_W_R),
      gamma_W_D = list(value = gamma_W_D_value, date = mu_gamma_H_date),
      gamma_W_R = list(value = gamma_W_R_value, date = mu_gamma_H_date))
    progression[k_parameters$parameter] <- k_parameters$value

    ## These could possibly be moved to the sircovid as defaults
    progression$k_sero_pre_1 <- 1
    progression$gamma_sero_pre_1 <- 1 / 13
    progression$k_sero_pre_2 <- 1
    progression$gamma_sero_pre_2 <- 1 / 13
    progression$k_PCR_pre <- 1
    progression$gamma_PCR_pre <- 0.1922243
    progression$k_PCR_pos <- 1
    progression$gamma_PCR_pos <- 0.083
    progression$k_sero_pos_1 <- 1
    progression$gamma_sero_pos_1 <- 1 / 200
    progression$k_sero_pos_2 <- 1
    progression$gamma_sero_pos_2 <- 1 / 400
    # Time to diagnosis if admitted without test
    progression$gamma_U <- 1 / 3

    observation <- sircovid::carehomes_parameters_observation()

    observation$rho_pillar2_tests <- rho_pillar2_tests
    observation$phi_pillar2_cases <- phi_pillar2_cases
    observation$phi_pillar2_cases_weekend <- phi_pillar2_cases_weekend
    observation$kappa_pillar2_cases <- kappa_pillar2_cases

    ## kappa for hospital data streams (not all will actually be used)
    observation$kappa_ICU <- 1 / alpha_H
    observation$kappa_general <- 1 / alpha_H
    observation$kappa_hosp <- 1 / alpha_H
    observation$kappa_admitted <- 1 / alpha_H
    observation$kappa_diagnoses <- 1 / alpha_H
    observation$kappa_all_admission <- 1 / alpha_H
    observation$kappa_death_hosp <- 1 / alpha_H

    ## kappa for death data streams (not all will actually be used)
    observation$kappa_death_carehomes <- 1 / alpha_D
    observation$kappa_death_comm <- 1 / alpha_D
    observation$kappa_death_non_hosp <- 1 / alpha_D
    observation$kappa_death <- 1 / alpha_D

    sens_and_spec <-
      sircovid::carehomes_parameters_sens_and_spec(sero_sensitivity_1 = 1,
                                                   sero_specificity_1 = 0.99,
                                                   sero_sensitivity_2 = 1,
                                                   sero_specificity_2 = 0.99,
                                                   pillar2_sensitivity = 1,
                                                   pillar2_specificity = 1,
                                                   react_sensitivity = 1,
                                                   react_specificity = 1)

    ## TODO: this breaks the current "fake" multistrain fitting and
    ## that needs working out. Previous version was fragile and
    ## confusing.
    if (multistrain) {
      rel_efficacy <- vaccination$efficacy
      strain_transmission <- c(1, strain_transmission_2)
    } else {
      ## We're adding a dimension for strain to all vaccine efficacy
      ## parameters this is in the format required for
      ## sircovid:::carehomes_parameters_vaccination
      f <- function(x) {
        ## TODO check this, and see if we can instead use mcstate's
        ## array functions.
        dim(x) <- c(19, 1, ncol(x))
        x
      }
      rel_efficacy <- lapply(vaccination$efficacy, f)
      strain_transmission <- 1
    }

    n_doses <- vaccination$schedule$n_doses
    if (n_doses == 2) {
      vaccine_progression_rate <- c(0, 1 / 21, 0, 0)
      vaccine_index_booster <- NULL
    } else if (n_doses == 3) {
      vaccine_progression_rate <- c(0, 1 / 21, 0, 0, 0)
      vaccine_index_booster <- 4L
    } else {
      stop("Expected 2 or 3 n_doses")
    }


    ret <- sircovid::carehomes_parameters(
      ## Core
      start_date = start_date,
      region = region,
      beta_date = beta_date,
      beta_value = beta_value,
      ## Severity
      severity = severity,
      ## Progression
      progression = progression,
      ## Observation
      observation = observation,
      ## sensitivity and specificity
      sens_and_spec = sens_and_spec,
      ## initial infective seed
      initial_I = 30,
      ## waning immunity rate
      waning_rate = waning_rate,
      ## Transmission
      eps = eps,
      m_CHW = m_CHW,
      m_CHR = m_CHR,
      p_NC = p_NC,
      p_NC_weekend = p_NC_weekend,
      rel_susceptibility = rel_efficacy$rel_susceptibility,
      rel_p_sympt = rel_efficacy$rel_p_sympt,
      rel_p_hosp_if_sympt = rel_efficacy$rel_p_hosp_if_sympt,
      rel_p_death = rel_efficacy$rel_p_death,
      rel_infectivity = rel_efficacy$rel_infectivity,
      vaccine_progression_rate = vaccine_progression_rate,
      vaccine_schedule = vaccination$schedule,
      vaccine_index_dose2 = 3L,
      vaccine_index_booster = vaccine_index_booster,
      n_doses = n_doses,
      ## Strains
      strain_transmission = strain_transmission,
      strain_seed_date = strain_seed_date,
      strain_seed_rate = strain_seed_rate,
      cross_immunity = cross_immunity)

    ## Could be moved to sircovid as a default
    ret$I_A_transmission <- 0.223

    ret$prop_strain_2 <- prop_strain_2

    ret
  }
}
