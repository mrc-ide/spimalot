% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_simulation.R
\name{spim_simulate_control_output}
\alias{spim_simulate_control_output}
\title{Create output control for simulations}
\usage{
spim_simulate_control_output(
  keep,
  time_series = TRUE,
  rt = TRUE,
  rt_weighted = rt,
  rt_type = NULL,
  state_by_age = FALSE,
  vaccination = TRUE
)
}
\arguments{
\item{keep}{Character vector of outputs to keep. Must be names of
variables as included in the index of the model run.}

\item{time_series}{Logical, indicating if a time series should be
output (these are quite large)}

\item{rt}{Logical, indicating if Rt should be calculated. This is
quite slow.}

\item{rt_weighted}{Logical, indicating if weighted Rt should be
calculated. This is also quite slow. Only allowed if \code{rt} is
\code{TRUE}}

\item{rt_type}{Character vecrtor of Rt types to use, passed to
\link[sircovid:lancelot_Rt]{sircovid::lancelot_Rt}.  Can be any or all of \code{eff_Rt_all},
\code{eff_Rt_general}, \code{Rt_all} and \code{Rt_general}.}

\item{state_by_age}{Logical, indicating if state by age should be
output. These are really large and \emph{will} cause memory issues
for you if you try to run with long time series unless you have
a lot of RAM.}

\item{vaccination}{Logical, indicating if vaccination information
should be output.}
}
\value{
An object of \code{spim_simulate_control_output}. Do not modify
this object after creation.
}
\description{
Control output created as a result of running simulations. Many of
the outputs are large, or require additional computation, so
reducing this will make your life more pleasant.
}
