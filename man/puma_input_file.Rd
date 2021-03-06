% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/puma_io.R
\name{puma_input_file}
\alias{puma_input_file}
\title{Construct PUMA input file}
\usage{
puma_input_file(f_name, f_path = getwd(), df, wd_min = -1,
  wd_max = -1, wd_step = 1)
}
\arguments{
\item{f_name}{puma input file name (extension -dat.txt will be added)}

\item{f_path}{directory where to save data}

\item{df}{data.frame(wavelength, intensity) with raw (measured) TR data}

\item{wd_min}{minimum wavelengt to use in calculations}

\item{wd_max}{maximum wavelength to use in calculation}

\item{wd_step}{wavelength step resolution}
}
\value{
data.frame with data saved to file
}
\description{
Write TR data in PUMA input file format
             file name should have extension '-dat.txt'
             first row contain points number
             data(wavelength, transmission/reflection)
             transmission in interval [0,1]
}
\examples{
       \dontrun{puma_input_file("ime", df)}
}
\author{
K. Juraic
}
