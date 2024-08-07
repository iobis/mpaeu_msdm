% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_matching.R
\name{name_to_aphia}
\alias{name_to_aphia}
\title{Retrieve AphiaID for one or more scientific names}
\usage{
name_to_aphia(
  scientificname,
  marine_only = TRUE,
  fuzzy = FALSE,
  try_valid = TRUE,
  only_exact = FALSE,
  progress = TRUE
)
}
\arguments{
\item{scientificname}{a vector with one or more scientific names}

\item{marine_only}{only look for marine species}

\item{fuzzy}{try fuzzy matching for the names. Recomended is \code{FALSE}}

\item{try_valid}{the function will always get the accepted name. If none is available
and \code{try_valid=TRUE}, then it will check if a "valid_AphiaID" value is
present, and then retrieve the record for that valid id.}

\item{only_exact}{if \code{TRUE}, it will get only exact matches}

\item{progress}{show progress bar (only used if more than 100 scientific names are supplied)}
}
\value{
dataframe with valid AphiaIDs
}
\description{
This function is a wrapper around \code{\link[worrms:wm_records_names]{worrms::wm_records_names()}}. It will search
for a valid AphiaID based on the supplie scientific name
}
\examples{
\dontrun{
  name_to_aphia("Acanthurus chirurgus")
}
}
