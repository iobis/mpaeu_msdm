% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outqc_get_info.R
\name{outqc_dist_tords}
\alias{outqc_dist_tords}
\title{Convert QC distances to RDS format}
\usage{
outqc_dist_tords(
  distfolder,
  do_parallel = TRUE,
  parallel_cores = NULL,
  remove_distances = TRUE
)
}
\arguments{
\item{distfolder}{the folder where the distances, calculated using \code{\link[=outqc_get_distances]{outqc_get_distances()}}, were saved}

\item{do_parallel}{if \code{TRUE} run in parallel using \code{\link[furrr:future_map]{furrr::future_map()}}}

\item{parallel_cores}{number of cores to use for parallel processing. If \code{NULL} uses half of available cores}

\item{remove_distances}{if \code{TRUE}, after converting to RDS the calculated distances are deleted,
and only the RDS folder is kept}
}
\value{
saved RDS files
}
\description{
Convert QC distances to RDS format
}
\examples{
\dontrun{
outqc_dist_tords("distances")
}
}
