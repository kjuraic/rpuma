% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/puma_run.R
\name{puma_run}
\alias{puma_run}
\title{run 3 iteration calculation with puma software}
\usage{
puma_run(puma_par, puma_path = "~/Job/R/Rpackages/rPuma/inst/puma",
  d_delta = 50)
}
\arguments{
\item{puma_par}{list with model parameters}

\item{puma_path}{puma program path}

\item{d_delta}{second iteration d interval}
}
\value{
rez = list with 3 iteration calculation results
}
\description{
run 3 iteration calculation with puma software
}
\examples{
    \dontrun{puma_run(puma_par)}
}
