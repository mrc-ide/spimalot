pkgload::load_all()

combined <- readRDS("outputs/combined.rds")
rt <- readRDS("outputs/rt_uk.rds")
vaccine_parameters <- readRDS("example/vaccine_parameters.rds")
mtp_commission <- read_csv("example/mtp_commission.csv")
npi_key <- read.csv("example/npi_key.csv", row.names = "npi")
end_date <- "2021-08-10"
n_par <- 10

options(error = recover)
obj <- spim_mtp_prepare(mtp_commission, npi_key, n_par, end_date,
                        vaccine_parameters, combined, rt)
res <- spim_mtp_simulate(obj)
