#' Get occurrence data from OBIS, GBIF and local dataset
#'
#' @param sp The scientific name of the species or the AphiaID.
#' @param database `character` vector containing the name(s) of the databases from
#' which data will be downloaded (currently "obis" and "gbif"). For getting data
#' from local datasets, add "local_*", being \* the name of the dataset (e.g. "local_belg").
#' This should match the name of the vector of local files.
#' @param local_file  named vector with the relative or absolute paths for the local
#' dataset files.
#' @param study_area a WKT format shapefile delimiting the study area. This will be used
#' to download the data from only the delimited area, which can speed up the dowload.
#' @param start_date the starting date for the records.
#' @param end_date the end date for the records.
#' @param save_results wether to save or not the results. 
#' If FALSE, a data.frame is returned with the final dataset.
#' @param save_folder folder where the resulting datasets should be save.
#' @param save_format which format to save the final dataset. One of "csv" or "parquet".
#' @param print.summary should a summary of the data processing be printed?
#' @param verbose should messages be printed?
#' @param ... additional parameters for the quality control step (see [mp_qc_check])
#'
#' @return
#' @export
#'
#' @examples
mp_get_occ <- function(sp,
                       database = c("obis", "gbif"),
                       local_file = NULL,
                       study_area = NULL,
                       start_date = NULL,
                       end_date = NULL,
                       save_results = T,
                       save_folder = "data/species/",
                       save_format = "parquet",
                       print.summary = F,
                       verbose = F,
                       ...
){
  
  if (verbose) cli_inform(c("i" = "Downloading data for {sp}.",
                            "{symbol$pointer} Data will be downloaded from {database}."))
  
  # Verify if range date is correct
  if (!is.null(range_date)) {
    if (length(range_date) < 2) {
      range_date <- c(range_date, Sys.Date())
    }
    range_date <- as.Date(range_date)
  }
  
  # Get data from each provider
  
  # Start a list to hold results
  occs <- list()
  
  # Assign databse type for local sources
  dbtype <- database
  dbtype[grepl("local", dbtype)] <- "local"
  
  for (i in 1:length(dbtype)) {
    
    if (verbose) cli_progress_step("Downloading from {database[i]}")
    
    # Obtain occurrence data from each provider
    occs[[i]] <- switch (dbtype[i],
                         obis = robis::occurrence(
                           scientificname = sp,
                           startdate = range_date[1],
                           enddate = range_date[2],
                           geometry = study_area
                         ),
                         gbif = as.data.frame(rgbif::occ_data(
                           scientificName = sp,
                           geometry = study_area,
                           eventDate = range_date
                         )$data),
                         local = .mp_get_occ_local(
                           local_file[database[i]],
                           scientificName = sp,
                           startdate = range_date[1],
                           enddate = range_date[2],
                           geometry = study_area
                         )
    )
    
    if (verbose) cli_progress_done()
    
  }
  
  # Assign names
  names(occs) <- database
  
  # Check for duplicates
  occs_nodups <- mp_dup_check(occs, exclude = T, as_single = T)
  
  # Perform quality control
  occs_qc <- mp_qc_check(occs_nodups, ...)
  
  if (print.summary) {
    .mp_print_occ_get(sp,
                      dbs = database,
                      full = occs,
                      dups = occs_nodups,
                      qc = occs_qc)
  }
  
  # Save
  if (save_results) {
    id <- worrms::wm_name2id(sp)
    
    if (any(save_format %in% c("parquet", "csv")) & length(save_format) == 1) {
      switch (save_format,
              parquet = write_parquet(occs_qc, glue("{save_folder}spdat_{id}.{save_format}")),
              csv = write_csv_arrow(occs_qc, glue("{save_folder}spdat_{id}.{save_format}"))
      )
    }
    
    return(invisible(NULL))
  } else {
    return(occs_qc)
  }
}