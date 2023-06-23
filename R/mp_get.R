#' Retrieve occurrence records from OBIS and save in the standard format
#'
#' @param sci_names a single value or a vector containing scientific names or, ideally, AphiaID (\code{numeric})
#' @param start_year the initial year from which occurrences should be get
#' @param save_format the format to save the files. Either "parquet" (recommended) or "csv"
#' @param mode mode of export. Either "single" (for single species) or "full" for multiple species. When using "full",
#' a single file containing all the records will be saved in the "data/raw" folder.
#' @param save_acro in case of full export, an acronym to be added to the save file. If NULL, "full" is used.
#'
#' @return Files saved
#' @export
#' 
#' @import arrow
#'
#' @examples
#' \dontrun{
#' 
#' }
mp_get_obis <- function(
  sci_names,
  start_year = 1950,
  save_format = "parquet",
  mode = "single",
  save_acro = NULL
) {
  
  if (mode == "single" & length(sci_names) > 1) {
    stop("When downloading multiple species you should use mode = 'full'")
  }
  
  # Download occurrences from OBIS
  if (is.character(sci_names)) {
    obis <- occurrence(sci_names,
                       startdate = as.Date(paste(start_year, 01, 01, sep = "-")))
    
    # If scientific name supplied, convert to AphiaID
    sci_names <- worrms::wm_name2id(sci_names)
    
  } else {
    obis <- occurrence(taxonid = sci_names,
                       startdate = as.Date(paste(start_year, 01, 01, sep = "-")))
  }
  
  # Save results
  if (nrow(obis) == 0) {
    stop("No occurrences found for this species.")
  } else {
    # Full save, for multiple species
    if (mode == "full") {
      fs::dir_create("data/raw/")
      save_folder <- "data/raw/"
      
      if (is.null(save_acro)) {
        save_acro <- "full"
      }
      
      # Save according to the chosen format
      if (any(save_format %in% c("parquet", "csv")) & length(save_format) == 1) {
        spath <- glue::glue("{save_folder}obis_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.{save_format}")
        switch (save_format,
                parquet = write_parquet(obis, spath),
                csv = readr::write_csv(obis, spath)
        )
      } else {
        stop("Check save format.")
      }
    
      # Individual species save  
    } else {
      save_folder <- glue::glue("data/species/key={sci_names}/date={format(Sys.Date(), '%Y%m%d')}/")
      fs::dir_create(save_folder)
      
      if (any(save_format %in% c("parquet", "csv")) & length(save_format) == 1) {
        spath <- glue::glue("{save_folder}type=obis.{save_format}")
        switch (save_format,
                parquet = write_parquet(obis, spath),
                csv = readr::write_csv(obis, spath)
        )
      } else {
        stop("Check save format.")
      }
    }
  }
  
  return(invisible(NULL))
  
}


#' Retrieve occurrence records from GBIF and save in the standard format
#'
#' @param sci_names a single value or a vector containing scientific names or, ideally, the taxon key (\code{numeric}).
#' Note that if scientific names are supplied, the function will first fetch the names with the GBIF taxonomic backbone. In case
#' it fails to find a match for one name, the function will return a list of non matched names. It is __strongly__ advised to use the taxon key.
#' @param start_year the initial year from which occurrences should be get
#' @param download_format the format to download data from GBIF. For now, only "SIMPLE_PARQUET" is possible.
#' @param remove_zip \code{logical} indicating if the downloaded zip should be removed (\code{TRUE}).
#' @param savelog \code{logical} indicating if a log file containing the GBIF download DOI should be saved (\code{TRUE}).
#' @param mode mode of export. Either "single" (for single species) or "full" for multiple species. When using "full",
#' a single file containing all the records will be saved in the "data/raw" folder
#' @param area_only if \code{TRUE}, then the scientific names are ignored, and the download will return
#' all points within the \code{geom_wkt}. Note that this only works when \code{mode = "full"}
#' @param save_acro in case of full export, an acronym to be added to the save file. If NULL, "full" is used.
#' @param geom_wkt a WKT geometry to limit the download. Should be simple, otherwise will fail
#'
#' @return Files saved and, if not all species matched the taxonomic backbone, a \code{data.frame} containing the non matched names.
#' @export
#' 
#' @import arrow
#'
#' @examples
#' \dontrun{
#' 
#' }
mp_get_gbif <- function(
    sci_names = NULL,
    start_year = 1950,
    download_format = "SIMPLE_PARQUET",
    remove_zip = T,
    savelog = T,
    mode = "full",
    area_only = FALSE,
    save_acro = NULL,
    geom_wkt = NULL
  ) {
  
  if (area_only & is.null(geom_wkt)) {
    stop("For area only you SHOULD supply a WKT.")
  }
  if (area_only & mode != "full") {
    stop("For area only, mode 'full' is the only currently supported.")
  }
  if (!area_only & is.null(sci_names)) {
    stop("You should supply scientific names when not using area only.")
  }
  
  if (!area_only) {
    # Check if names were supplied or GBIF codes
    if (is.character(sci_names)) {
      gbif_taxon_keys <- name_backbone_checklist(sci_names, strict = T)
      
      if (sum(gbif_taxon_keys$matchType == "NONE") > 0) {
        warning("GBIF failed to match names for ",
                sum(gbif_taxon_keys$matchType == "NONE"),
                " species. Check the output carefully.")
        non_match <- sci_names[gbif_taxon_keys$matchType == "NONE"]
        
        gbif_taxon_keys <- gbif_taxon_keys[gbif_taxon_keys$matchType != "NONE",]
      } else {
        non_match <- NULL
      }
      
      sel_keys <- gbif_taxon_keys$usageKey
      
    } else {
      sel_keys <- sci_names
      
      non_match <- NULL
    }
  }
  
  # Set the download on GBIF server
  if (!is.null(geom_wkt)) {
    if (area_only) {
      gbif_download <- occ_download(
        pred_within(geom_wkt),
        pred("hasCoordinate", TRUE),
        pred_gte("year", start_year),
        pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
        format = download_format
      )
    } else {
      gbif_download <- occ_download(
        pred_in("taxonKey", sel_keys),
        pred("hasCoordinate", TRUE),
        pred_gte("year", start_year),
        pred_within(geom_wkt),
        pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
        format = download_format
      )
    }
  } else {
    gbif_download <- occ_download(
      pred_in("taxonKey", sel_keys),
      pred("hasCoordinate", TRUE),
      pred_gte("year", start_year),
      pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
      format = download_format
    )
  }
  
  # Check the download status
  status <- occ_download_meta(gbif_download)$status
  cli::cli_progress_bar("Waiting for download...")
  while (status == "RUNNING" | status == "PREPARING") {
    Sys.sleep(5)
    cli::cli_progress_update()
    status <- try(occ_download_meta(gbif_download)$status)
    if (class(status) == "try-error") {
      Sys.sleep(5)
      status <- try(occ_download_meta(gbif_download)$status)
    }
  }
  cli::cli_progress_update(force = T)
  
  if (status != "SUCCEEDED") {
    stop("Downloaded failed.")
  } else {
    if (mode == "full") {
      fs::dir_create("data/raw/")
      dl <- occ_download_get(gbif_download, path = "data/raw")
      
      if (is.null(save_acro)) {
        save_acro <- "full"
      }
      
      unzip(paste0("data/raw/", gbif_download, ".zip"),
            exdir = "data/raw/")
      
      switch (download_format,
        SIMPLE_CSV = file.rename("data/raw/occurrence.csv", glue::glue("data/raw/gbif_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.csv")),
        SIMPLE_PARQUET = file.rename("data/raw/occurrence.parquet", glue::glue("data/raw/gbif_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.parquet"))
      )
      
      if (remove_zip) {
        fs::file_delete(paste0("data/raw/", gbif_download, ".zip"))
      }
      
      if (savelog) {
        write.table(
          data.frame(date = attr(gbif_download, "created"),
                     download_code = as.character(gbif_download),
                     doi = attr(gbif_download, "doi"),
                     citation = attr(gbif_download, "citation")),
          glue::glue("data/gbif_{save_acro}_download_log.txt"), row.names = F
        )
      }
    } else {
      
      if (is.character(sci_names)) {
        spnam <- sci_names
      } else {
        spnam <- name_usage(sci_names)
        spnam <- spnam$data$scientificName
      }
      
      spnam <- unlist(stringr::str_split(spnam, pattern = " "))
      spnam <- paste(spnam[1:2], collapse = " ")
      spnam <- worrms::wm_records_name(spnam)
      spnam <- spnam[spnam$status == "accepted", "AphiaID"][1,]
      spnam <- spnam$AphiaID
      
      sppath <- glue::glue("data/species/key={spnam}/date={format(Sys.Date(), '%Y%m%d')}")
      
      fs::dir_create(sppath)
      dl <- occ_download_get(gbif_download, path = sppath)
      
      unzip(paste0(sppath, "/", gbif_download, ".zip"),
            exdir = sppath)
      
      if (download_format == "SIMPLE_PARQUET") {
        da <- open_dataset(paste0(sppath, "/occurrence.parquet"))
      } else {
        da <- data.table::fread(paste0(sppath, "/occurrence.csv"))
      }
      
      write_parquet(da, paste0(sppath, "/type=gbif.parquet"))
      rm(da)
      
      switch (download_format,
        SIMPLE_PARQUET = fs::dir_delete(paste0(sppath, "/occurrence.parquet")),
        SIMPLE_CSV = fs::dir_delete(paste0(sppath, "/occurrence.csv"))
      )
      
      if (remove_zip) {
        fs::file_delete(paste0(sppath, "/", gbif_download, ".zip"))
      }
      
      if (savelog) {
        write.table(
          data.frame(date = attr(gbif_download, "created"),
                     doi = attr(gbif_download, "doi")),
          paste0(sppath, "/type=log_gbif.txt"), row.names = F
        )
      }
    }
  }
  
  if (!is.null(non_match)) {
    return(non_match)
  } else {
    return(invisible(NULL))
  }
  
}


#' Retrieve occurrence records from local file for individual species and save in the standard format
#'
#' @param local_file the relative or absolute path to the file (or folder in case of parquet datasets)
#' @param database the name of the database source, will be used to save the files
#' @param sci_name the scientific name of the species. If \code{NULL}, then the aphia_id or gbif_key will be used instead (recommended)
#' @param aphia_id the AphiaID of the species
#' @param gbif_key the GBIF taxon key
#' @param save_format the format to save the file (either "parquet" or "csv"; recommended is "parquet")
#'
#' @return individual file saved
#' @export
#'
#' @import arrow 
#'
#' @examples
#' \dontrun{
#' 
#' }
mp_get_local <- function(local_file,
                         database,
                         sci_name = NULL,
                         aphia_id = NULL,
                         gbif_key = NULL,
                         save_format = "parquet"
                         ) {
  
  if (is.null(sci_names) & is.null(aphia_id) & is.null(gbif_key)) {
    stop("One of sci_names, aphia_id or gbif_key should be supplied.")
  }
  
  if (tools::file_ext(local_file) == "parquet") {
    database <- open_dataset(local_file)
  } else {
    database <- open_csv_dataset(local_file)
  }
  
  if (is.null(aphia_id) & is.null(gbif_key)) {
    # Filter by scientific name
    sp_data <- database %>%
      filter(scientificName == sci_names) %>%
      collect()
    
  } else {
    if (is.null(aphia_id)) {
      # Filter by taxon key (GBIF)
      sp_data <- database %>%
        filter(taxonkey == gbif_key) %>%
        collect()
      
    } else {
      # Filter by AphiaID (OBIS)
      sp_data <- database %>%
        filter(AphiaID == aphia_id) %>%
        collect()
    }
  }
  
  if (any(save_format %in% c("parquet", "csv")) & length(save_format) == 1) {
    
    if (is.null(aphia_id)) {
      if (is.null(gbif_key)) {
        spnam <- worrms::wm_name2id(sci_name)
      } else {
        spnam <- name_usage(gbif_key)
        spnam <- spnam$data$scientificName
        spnam <- unlist(stringr::str_split(spnam, pattern = " "))
        spnam <- paste(spnam[1:2], collapse = " ")
        spnam <- worrms::wm_records_name(spnam)
        spnam <- spnam[spnam$status == "accepted", "AphiaID"][1,]
        spnam <- spnam$AphiaID
      }
    } else {
      spnam <- aphia_id
    }
    
    save_folder <- glue::glue("data/species/key={spnam}/date={format(Sys.Date(), '%Y%m%d')}")
    fs::dir_create(save_folder)
    
    spath <- glue::glue("{save_folder}/type={database}.{save_format}")
    switch (save_format,
            parquet = write_parquet(obis, spath),
            csv = readr::write_csv(obis, spath)
    )
  } else {
    stop("Check save format.")
  }
  
  return(invisible(NULL))
  
}