% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_get_conf.R
\name{get_envofgroup}
\alias{get_envofgroup}
\title{Load SpatRaster environmental layers for a group defined on a configuration
file}
\usage{
get_envofgroup(
  group,
  depth = "depthsurf",
  scenario = "current",
  period = NULL,
  hypothesis = NULL,
  load_all = FALSE,
  conf_file = "sdm_conf.yml",
  env_folder = "data/env",
  future_name = "future",
  fixed_name = "terrain",
  keep_in_future = TRUE,
  surface_ifnot = TRUE,
  accepted_formats = c("tif", "nc"),
  verbose = TRUE
)
}
\arguments{
\item{group}{the group for which you want to retrieve the layers}

\item{depth}{the depth of the layer (usually one of "surf", "mean", "min" or
"max")}

\item{scenario}{the scenario for the layers. Current or one of the
future/past. See details}

\item{period}{optional character indicating the period}

\item{hypothesis}{if \code{NULL}, it will retrieve the first available
hypothesis list (or all if \code{load_all=TRUE, see below}).
Otherwise, it will look for the list equivalent to the
hypothesis. Should be a character (e.g. "hypothesis1")}

\item{load_all}{if \code{TRUE} it will load the unique files from all
available hypothesis for that group.}

\item{conf_file}{the path for the YAML configuration file}

\item{env_folder}{the path for the environmental layers folder}

\item{fixed_name}{the name of the folder containing the fixed or non variable
layers (usually terrain ones). Can also be NULL if you don't have fixed layers.}

\item{keep_in_future}{if \code{TRUE}, in case a future layer is not
available, it will get its current period version}

\item{surface_ifnot}{if \code{TRUE}, it will use the surface layer in case it
is not found for the required depth. If you named your surface layers in a way different than "depthsurf",
then you can instead supply here the name of the surface depth.}

\item{accepted_formats}{a character vector of accepted formats. Those should be supported by \code{\link[terra:rast]{terra::rast()}}}

\item{verbose}{print messages}
}
\value{
a SpatRaster with as many layers as present in the configuration file
}
\description{
Load SpatRaster environmental layers for a group defined on a configuration
file
}
\details{
Although you can indicate any path for the env folder, this folder
must be organized internally in the following way:
\itemize{
\item current (for the primary fitting scenario it SHOULD be named current)
\item \code{future_name} (the name of the folder for future layers)
\item {scenario1} (folder containing layers for the scenario1)
\item {scenario2} (folder containing layers for the scenario2)
\item \code{fixed_name} (the name of the folder for fixed layers, i.e., layers that are non variable across periods like bathymetry)
}

Note that you can use any folder name for your scenarios (e.g. ssp1, ssp2).
The function is not case sensitive, and will convert everything to lower
case.

The function will automatically detect if the layer is fixed or not, by
searching in the fixed folder.
}
\section{Using user supplied layers vs package layers}{
If you used the function \code{\link[=get_envlayers]{get_envlayers()}} to download from the Bio-ORACLE,
then the function will work without any problem. However, if you added
additional layers take into account that the function expect that layers names contains three pieces of information
separated by underscores or dashes: name of variable, depth,
and variant (e.g. max). One example: siconc_baseline_depthsurf_mean. Note that here you have
additional information, what is fine if you have the essential information. For fixed ones,
you can have only the name (e.g. bathymetry).

Because you can have "depthmean" and a variant "mean", it is ESSENTIAL that the variant comes preceded by _ or -. For example, all those should be valid:
\itemize{
\item siconc_baseline_depthmean_mean
\item siconc_mean_baseline_depthmean
\item siconc-baseline-mean-depthmean
}

But this one not work:
\itemize{
\item mean_siconc_baseline_depthmean
}

Layers can have an additional information regarding the period. So, for example, future layers for the decade of 2050 can be named:
siconc_ssp126_depthsurf_dec50_mean

The function will do strict string match for everything, so ensure that you wrote the names correctly in your configuration file.
}

\section{Note}{
Layer names are changed to the retrieved combination of name + variation (e.g. thetao_mean). The old names are stored as an attribute
called "old_names" and can be accessed using \code{attr(raster_layers_object, "old_names")} for verification.
}

\examples{
\dontrun{
raster_layers <- get_envofgroup(group = "fishes")
}
}
\seealso{
\code{\link[=get_listbygroup]{get_listbygroup()}}, \code{\link[=get_conf]{get_conf()}}
}
