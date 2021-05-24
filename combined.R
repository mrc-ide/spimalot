pkgload::load_all()
dat <- spim_combined_load("example/fits")

nowcast <- spimalot::spim_summary_nowcast(dat)
## save_results(nowcast,
##              "template_combined.xlsx",
##              "outputs/spim_output.xlsx")

time_series <- spimalot::spim_summary_time_series(dat, "2020-10-01")
## save_results(time_series,
##              "template_combined.xlsx",
##              "outputs/time-series.xlsx")
