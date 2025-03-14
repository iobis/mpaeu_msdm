% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_processing.R
\name{post_prepare}
\alias{post_prepare}
\title{Prepare files for use in Zonation}
\usage{
post_prepare(
  species,
  source_folder,
  study_area,
  model_acro = "mpaeu",
  target_mask = "fit_region",
  target_threshold = "p10",
  target_model = c("maxent", "rf", "xgboost", "ensemble", "esm"),
  future_scenarios = c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
  future_periods = c("dec50", "dec100"),
  boot_weight = 0.5,
  outfolder = "proc-layers",
  verbose = TRUE
)
}
\arguments{
\item{species}{a vector of species IDs (AphiaID)}

\item{source_folder}{the path to the results folder, typically "results"}

\item{study_area}{the path to the study area shapefile (any format readable by \code{\link[terra:vect]{terra::vect()}})}

\item{model_acro}{the modelling realization acronym (e.g. "mpaeu")}

\item{target_mask}{which mask layer/type to use. Should be available on the mask file}

\item{target_threshold}{which threshold to use to remove areas with lower probability of occurrence}

\item{target_model}{which model (algorithm) to be used. Can be a single value (e.g. "maxent") or a vector.
In case a vector is supplied, the function will check which of those are available (in the order supplied)
and then select the first one available. So, if you supply c("maxent", "rf", "ensemble") and "maxent" is not available,
but the others are, the function will use Random Forest ("rf")}

\item{future_scenarios}{which scenarios you want to also process (in the format "sspXXX"). Current scenario is always included.}

\item{future_periods}{which periods (in the format "decXX") to use. Current period is always included.}

\item{boot_weight}{weight to be applied to the bootstrap layer}

\item{outfolder}{folder to save the processed layers. If not existing, it will be created}

\item{verbose}{if \code{TRUE} print messages}
}
\value{
processed files and a data.frame with the processing status
}
\description{
Prepare files for use in Zonation
}
\examples{
\dontrun{
post_prepare(species = c(124287, 137098), 
             source_folder = "results", 
             study_area = "data/shapefiles/mpaeu_studyarea_v2.shp", 
             model_acro = "mpaeu")
}
}
