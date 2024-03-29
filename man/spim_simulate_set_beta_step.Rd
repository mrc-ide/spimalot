% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_simulation.R
\name{spim_simulate_set_beta_step}
\alias{spim_simulate_set_beta_step}
\title{Add beta_step into control}
\usage{
spim_simulate_set_beta_step(control, beta_step)
}
\arguments{
\item{control}{A control object from \link{spim_simulate_control}}

\item{beta_step}{Beta values; either a named list (if \code{beta_step}
appears in the run grid, which it probably will) or a single
set. These should be 3d arrays with dimensions corresponding to
\enumerate{
\item particle, 2. region, and 3. step (so time multiplied by
\code{steps_per_day})
}}
}
\description{
Add \code{beta_step} into a \code{spim_simulate_control} object.  Typically
this happens just before running simulations, after the Rt values
have been converted into beta values according to the assumptions
of the simulation.
}
