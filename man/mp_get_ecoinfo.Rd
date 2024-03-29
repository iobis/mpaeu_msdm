% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mp_get_ecoinfo.R
\name{mp_get_ecoinfo}
\alias{mp_get_ecoinfo}
\title{Get habitat information for species}
\usage{
mp_get_ecoinfo(
  species_list,
  outfile = "data/species_ecoinfo.csv",
  overwrite = FALSE,
  return_table = FALSE,
  try_higher = TRUE,
  try_remarks = TRUE,
  show_progress = TRUE
)
}
\arguments{
\item{species_list}{a vector with AphiaIDs}

\item{outfile}{a file to save the habitat information. Can be set to \code{NULL}
to not save a file}

\item{overwrite}{if this is set to \code{FALSE}, then the function will first read
the outfile (if it exists), then it will see which records already exists in the file,
and finally will get the records for those non existent}

\item{return_table}{if \code{TRUE}, return a \code{data.frame} with the information}

\item{try_higher}{if \code{TRUE} try to find the information from higher taxons}

\item{try_remarks}{if \code{TRUE} try to find the information on the remarks}

\item{show_progress}{if \code{TRUE} show a progress bar}
}
\value{
data.frame or saved file
}
\description{
Get habitat information for species
}
\details{
Species information is retrieved from SeaLifeBase, FishBase and WoRMS (in that order).
If the information exists in multiple sources, the information of the first is returned.
The function will by default search the information in higher taxonomic levels (on WoRMS)
if nothing was found, and also search on the remarks (additional information) of
SeaLifeBase and FishBase.

If nothing is found, the function returns a value of NOT_FOUND for that record.
\subsection{Behavior when list is long}{

If your species list have more than 50 records, the function will first create
a \code{\link[storr:storr]{storr::storr()}} database (using RDS structure) and then populate it with the
records. That means that, if for some reason the function fails at any point
(e.g. by internet connection problems) or if you need to abort it, you can just call
the function again and it will start from where it stopped.

The function will only remove the database (stored in a folder "_species_info")
when it completes the cycle.
}
}
\examples{
mp_get_ecoinfo(c(367850, 513377), outfile = NULL, return_table = TRUE)
}
