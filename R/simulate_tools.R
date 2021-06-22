spim_simulate <- function(i, args, combined) {
  el <- args[[i]]
  message(sprintf("-----\nRunning scenario %d / %d", i, length(args)))
  time <- system.time(
    ret <- spim_simulate_one(el, combined))
  message(sprintf("Finished scenario %d in %2.1f s", i, time[["elapsed"]]))
  ret
}
