% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pars.R
\name{spim_pars_pmcmc_load}
\alias{spim_pars_pmcmc_load}
\title{Load a set of parameters for the pmcmc}
\usage{
spim_pars_pmcmc_load(
  path,
  info = "info.csv",
  prior = "prior.csv",
  proposal = "proposal.csv"
)
}
\arguments{
\item{path}{Directory where the csv files are found.}

\item{info}{Filename for the parameter info, relative to \code{path}}

\item{prior}{Filename for the parameter priors, relative to \code{path}}

\item{proposal}{Filename for the parameter proposal, relative to \code{path}}
}
\value{
A spim_pars_pmcmc object
}
\description{
Load the pmcmc parameters from disk. We expect three files; one
for the overall parameters (\code{info}), one with details of the priors
(\code{prior}) and one describing the proposal kernel (\code{proposal}).
}
