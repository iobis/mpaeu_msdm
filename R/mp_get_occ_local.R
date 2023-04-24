#' Get occurrence data from a local dataset
#'
#' @description
#' Get occurrence data from local datasets. Those have to contain, at least, 
#' one collumn called "scientificName" and two columns with longitude/latitude 
#' data coded as lon/lat, x/y or decimalLongitude/decimalLatitude. 
#' Date columns should be named "date". Additional columns are ignored.
#'
#' @param loc_file the absolute or relative path to the local dataset file
#' @param scientificName the name of the species
#' @param startdate an optional start date to filter occurrences (YYYY-MM-DD)
#' @param enddate an optional end date to filter occurrences (YYYY-MM-DD)
#' @param geometry an optional geometry (WKT format) to filter occurrences
#'
#' @return the filterd locations
#'
#' @examples
#' 
.mp_get_occ_local <- function(loc_file,
                              scientificName,
                              startdate = NULL,
                              enddate = NULL,
                              geometry = NULL) {
  
  dat <- switch(tools::file_ext(loc_file),
                "csv" = read.csv(loc_file),
                "parquet" = read_parquet(loc_file))
  
  colnames(dat) <- gsub("lon|longitude|x|decimalLongitude", "decimalLongitude", colnames(dat))
  colnames(dat) <- gsub("lat|latitude|y|decimalLatitude", "decimalLatitude", colnames(dat))
  
  dat <- st_as_sf(dat,
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = "EPSG:4326")
  
  dat_filtered <- dat %>%
    filter(scientificName == scientificName)
  
  if (!is.null(geometry)) {
    dat_filtered <- st_filter(dat_filtered, geometry)
  }
  
  if (!is.null(startdate) | !is.null(enddate)) {
    if (is.null(startdate)) {
      startdate <- as.Date("1950-01-01")
    }
    if (is.null(enddate)) {
      enddate <- Sys.Date()
    }
    
    dat_filtered <- dat_filtered %>%
      mutate(date = lubridate::as_date(date))
      filter(date >= startdate & date <= enddate)
  }
  
  return(dat_filtered)
}