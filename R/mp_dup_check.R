#' Remove duplicate datasets
#'
#' @description
#' Removes duplicate records based on geohash and year.
#'
#' @param data_a A named `list` of datasets which are to be compared or a `data.frame`.
#' Should have at least the collumns ScientificName, decimalLongitude,
#' decimalLatitude and year.
#' @param data_b (data.frame) If `data_a` is a `data.frame`, then `data_b` should be supplied.
#' @param gh_precision precision parameter for the geohash.
#' (see details at [geohashTools::gh_encode()]).
#'
#' @return A `data.frame` with cleaned occurrence records.
#' @export
#'
#' @examples
#' 
#' 
mp_dup_check <- function(data_a,
                         data_b = NULL,
                         gh_precision = 6) {
  
  if (!is.list(data_a)) {
    if (is.null(data_b)) {
      stop("Second dataset was not supplied. data_a should be a list of datasets or data_b should be supplied.")
    } else {
      data_a <- list(data_a = data_a,
                     data_b = data_b)
    }
  }
  
  # Select essential variables that should be kept in the final dataset
  essent_vars <- c("occurrenceID", "catalogNumber", "recordNumber",
                   "fieldNumber", "materialSampleID", "institutionID",
                   "collectionID", "datasetID", "collectionCode",
                   "institutionCode", "datasetName")
  
  # Filter and prepare data for duplicate checking
  data_a <- lapply(1:length(data_a), function(x){
    data_a[[x]] %>%
      select(scientificName, decimalLongitude, decimalLatitude, year, any_of(essent_vars)) %>%
      mutate(datasetName = str_trunc(datasetName, 40)) %>%
      mutate(year = lubridate::as_date(ifelse(
        is.na(year), NA, glue::glue("{year}-01-01")
      ))) %>%
      mutate(year = lubridate::year(year)) %>%
      filter(!is.na(year)) %>%
      filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
      mutate(dt_proj_id = names(data_a)[x])
  })
  
  data_a <- bind_rows(data_a)
  
  
  # Get geohash and unique code by gh + year
  data_f <- data_a %>%
    mutate(geohash = gh_encode(decimalLatitude, decimalLongitude, gh_precision)) %>%
    mutate(cell = factor(paste(geohash,
                               year, sep = "_"))) %>%
    distinct(cell, .keep_all = T) %>%
    select(-cell, -geohash)
  
  return(data_f)

}
