% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_plots.R
\name{spim_scenario_cols}
\alias{spim_scenario_cols}
\title{Return accessible scenario colours}
\usage{
spim_scenario_cols(
  scenarios,
  weight = 0.3,
  palette = "bright",
  preview = FALSE
)
}
\arguments{
\item{scenarios}{Unique scenario names. Scenarios with \verb{[High R]} will
darken the colour for the corresponding central scenario and scenarios
with \verb{[Low R]} will brighten the colour for the corresponding
central scenario.}

\item{weight}{\code{weight} passed to \code{mix_cols} to darken/brighten colours for
high/low R scenarios}

\item{palette}{Colour palette, passed to \link[khroma:colour]{khroma::colour}}

\item{preview}{If \code{TRUE} then plots the final colour scheme with
\link[khroma:plot_scheme]{khroma::plot_scheme}}
}
\description{
Select colours for plotting simulation scenarios
}
\examples{
spim_scenario_cols(c("Step 4", "Step 4 [High R]", "Step 4 [Low R]"))

spim_scenario_cols(c("Step 3", "Step 3 [High R]", "Step 4"))

}
