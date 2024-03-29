% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{spim_simulation_shaps}
\alias{spim_simulation_shaps}
\title{Calculate SHAPs over predicted states}
\usage{
spim_simulation_shaps(summary, feats = NULL)
}
\arguments{
\item{summary}{A tidy summary object such as that returned by
\code{create_summary}}

\item{feats}{Features to calculate SHAPS for. If NULL then uses default
selection returned by \code{spim_simulation_predictors}}
}
\description{
Calculate SHAPs from a tidy summary of predictions over various features.
SHAPs calculated as the expected difference in predicted states with and
without the given feature of interest (over all feature levels).
}
