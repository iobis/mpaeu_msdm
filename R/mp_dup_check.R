#' Remove duplicate datasets
#'
#' @description Removes duplicate records based on geohash and year.
#'
#' @param data_a A named `list` of datasets which are to be compared or a
#'   `data.frame`. Should have at least the collumns ScientificName,
#'   decimalLongitude, decimalLatitude and year.
#' @param data_b (data.frame) If `data_a` is a `data.frame`, then `data_b`
#'   should be supplied in order to produce a comparison between data providers. 
#'   Otherwise, the comparison will be only within the dataset.
#' @param gh_precision precision parameter for the geohash. (see details at
#'   [geohashTools::gh_encode()]).
#'
#' @return A `data.frame` with cleaned occurrence records.
#' @export
#'
#' @import dplyr
#' 
#' @author Pieter Provoost (adapted by Silas Principe)
#'
#' @examples
#' \dontrun{
#' mp_dup_check(dataset_a, dataset_b)
#' }
#' 
mp_dup_check <- function(data_a,
                         data_b = NULL,
                         gh_precision = 6) {
  
  if (is.data.frame(data_a) | any("tbl" %in% class(data_a))) {
    if (!is.null(data_b)) {
      data_a <- list(data_a = data_a,
                     data_b = data_b)
    } else {
      data_a <- eval(parse(text = paste0("list(", substitute(data_a), "= data_a)")))
    }
  } else {
    if (is.null(names(data_a))) {
      names(data_a) <- c("data_a", "data_b")
    } else {
      if (any(!nzchar(names(data_a)))) {
        names(data_a) <- c("data_a", "data_b")
      }
    }
  }
  
  # Select essential variables that should be kept in the final dataset
  essent_vars <- c("occurrenceID", "catalogNumber", "recordNumber",
                   "fieldNumber", "materialSampleID", "institutionID",
                   "collectionID", "datasetID", "collectionCode",
                   "institutionCode", "datasetName",
                   "gbifid", "datasetkey", "eventID")
  
  # Filter and prepare data for duplicate checking
  data_a <- lapply(1:length(data_a), function(x){
    data_a[[x]] %>%
      {if ("date_year" %in% colnames(data_a[[x]])) mutate(., year = date_year) else .} %>%
      select(species, decimalLongitude, decimalLatitude, year, any_of(essent_vars)) %>%
      {if ("datasetName" %in% colnames(data_a[[x]])) mutate(., datasetName = stringr::str_trunc(datasetName, 40)) else .} %>%
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
    mutate(geohash = geohashTools::gh_encode(decimalLatitude, decimalLongitude, gh_precision)) %>%
    mutate(cell = factor(paste(geohash,
                               year, sep = "_"))) %>%
    distinct(cell, .keep_all = T) %>%
    select(-cell, -geohash)
  
  return(data_f)

}
