% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdm_tools.R
\name{save_log}
\alias{save_log}
\title{Save log object}
\usage{
save_log(log_object, file_path)
}
\arguments{
\item{log_object}{log object generated with \code{\link[=gen_log]{gen_log()}}}

\item{file_path}{the path to save the file in .json format}
}
\value{
saved file
}
\description{
Save log object
}
\examples{
\dontrun{
save_log(sdm_log, "species_log.json")
}
}
