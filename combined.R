pkgload::load_all(export_all = FALSE)

dat <- spimalot::spim_combined_load("example/fits")

dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)

nowcast <- spimalot::spim_summary_nowcast(dat)
spimalot::spim_summary_write(nowcast, "example/template_combined.xlsx",
                             "outputs/spim_output.xlsx")

time_series <- spimalot::spim_summary_time_series(dat, "2020-10-01")
spimalot::spim_summary_write(time_series, "example/template_combined.xlsx",
                             "outputs/time_series.xlsx")

spimalot::spim_pars_pmcmc_save(dat$parameters, "outputs/parameters")
saveRDS(dat$onward, "outputs/combined.rds")
saveRDS(dat$ifr_t, "outputs/ifr_t.rds")

png("figs/forest_plot.png", width = 2400, height = 1600, res = 200)
spim_plot_forest(dat, FALSE)
dev.off()

png("figs/forest_plot_betas.png", width = 2400, height = 1600, res = 200)
spim_plot_forest(dat, TRUE)
dev.off()

pkgload::load_all(export_all = FALSE)
png("figs/data_fits_regional.png", width = 2400 / 5 * 7, height = 1800,
    res = 200)
spimalot::spim_plot_trajectories(
  dat, sircovid::regions("england"),
  c("deaths_hosp", "deaths_carehomes", "deaths_comm", "icu", "general",
    "hosp", "all_admission"),
  with_forecast = FALSE, add_betas = FALSE)
dev.off()

png("figs/data_fits_national.png", width = 2400, height = 1800, res = 200)
spimalot::spim_plot_trajectories(
  dat, c(sircovid::regions("nations"), "uk"),
  c("deaths_hosp", "deaths_carehomes", "deaths", "icu", "general",
    "hosp", "all_admission"),
  with_forecast = FALSE, add_betas = FALSE)
dev.off()

png("figs/projections_regional.png", width = 2400 / 5 * 7, height = 1800,
    res = 200)
spimalot::spim_plot_trajectories(
  dat, sircovid::regions("england"),
  c("deaths_hosp", "deaths_carehomes", "deaths_comm", "icu", "general",
    "hosp", "all_admission"),
  date_min = dat$info$date - 42, with_forecast = TRUE, add_betas = FALSE)
dev.off()

png("figs/projections_national.png", width = 2400, height = 1800, res = 200)
spimalot::spim_plot_trajectories(
  dat, c(sircovid::regions("nations"), "uk"),
  c("deaths_hosp", "deaths_carehomes", "deaths", "icu", "general",
    "hosp", "all_admission"),
  date_min = dat$info$date - 42, with_forecast = TRUE, add_betas = FALSE)
dev.off()

png("figs/Rt_eff_all.png", width = 2400, height = 1200, res = 200)
spimalot::spim_plot_Rt(dat, "eff_Rt_all")
dev.off()

png("figs/Rt_all.png", width = 2400, height = 1200, res = 200)
spimalot::spim_plot_Rt(dat, "Rt_all")
dev.off()

png("figs/Rt_eff_general.png", width = 2400, height = 1200, res = 200)
spimalot::spim_plot_Rt(dat, "eff_Rt_general")
dev.off()

png("figs/Rt_general.png", width = 2400, height = 1200, res = 200)
spimalot::spim_plot_Rt(dat, "Rt_general")
dev.off()

png("figs/serology_euroimmun.png", width = 2400, height = 1200, res = 200)
spimalot::spim_plot_serology(dat, 1, 40)
dev.off()

png("figs/serology_roche_n.png", width = 2400, height = 1200, res = 200)
spimalot::spim_plot_serology(dat, 2, 40)
dev.off()



## add (zoomed in) plots of SPI-M-relevant trajectories
dir.create("spim_view", FALSE, TRUE)

png("spim_view/regions.png", width = 2400 / 5 * 7, height = 1800, res = 200)
spimalot::spim_plot_trajectories(
  dat, sircovid::regions("england"),
  c("deaths", "hosp", "all_admission"), date_min = as.Date(dat$info$date) - 45,
  with_forecast = FALSE, add_betas = TRUE)
dev.off()

png("spim_view/nations.png", width = 2400, height = 1800, res = 200)
spimalot::spim_plot_trajectories(
  dat, c(sircovid::regions("nations"), "uk"),
  c("deaths_hosp", "deaths_carehomes", "deaths", "icu", "general",
    "hosp", "all_admission"),
  with_forecast = FALSE, add_betas = FALSE)
dev.off()
