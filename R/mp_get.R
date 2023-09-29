############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################

#' Retrieve occurrence records from OBIS and save in the standard format
#'
#' @param sci_names a single value or a vector containing scientific names or, ideally, AphiaID (\code{numeric})
#' @param start_year the initial year from which occurrences should be get.
#' @param save_format the format to save the files. Either "parquet" (recommended) or "csv"
#' @param mode mode of export. Either "single" (for single species) or "full" for multiple species. When using "full",
#' a single file containing all the records will be saved in the "data/raw" folder.
#' @param save_acro in case of full export, an acronym to be added to the save file. If NULL, "full" is used.
#' @param geom_wkt a WKT geometry. Can be supplied together with
#'   \code{sci_names}, in which case data will be downloaded only for this area,
#'   or alone, in which case all occurrence records in a certain area will be
#'   downloaded. To use the geometry alone, \code{mode} should be \code{"full"}.
#'   Should be simple, otherwise will fail.
#'
#' @return Files saved
#' @export
#' 
#' @import arrow
#'
#' @examples
#' \dontrun{
#' mp_get_obis("Leptuca thayeri")
#' }
mp_get_obis <- function(
  sci_names = NULL,
  start_year = 1950,
  save_format = "parquet",
  mode = "single",
  save_acro = NULL,
  geom_wkt = NULL
) {
  
  if (is.null(sci_names) & is.null(geom_wkt)) {
    stop("You should supply at least one of sci_names or geom_wkt")
  }
  
  if (is.null(sci_names) & mode != "full") {
    stop("When using geometry alone, you should use mode 'full'")
  }
  
  if (is.null(sci_names)) {
    
    obis <- robis::occurrence(geometry = geom_wkt,
                       startdate = as.Date(paste(start_year, 01, 01, sep = "-")))
    
  } else {
    
    if (mode == "single" & length(sci_names) > 1) {
      stop("When downloading multiple species you should use mode = 'full'")
    }
    
    # Download occurrences from OBIS
    if (is.character(sci_names)) {
      sci_names <- worrms::wm_name2id(sci_names)
    }
    
    obis <- robis::occurrence(taxonid = sci_names,
                       startdate = as.Date(paste(start_year, 01, 01, sep = "-")))
    
  }
  
  # Save results
  if (nrow(obis) == 0) {
    stop("No occurrences found for this species or area.")
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
#' @param sci_names a single value or a vector containing scientific names or,
#'   ideally, the taxon key (\code{numeric}). Note that if scientific names are
#'   supplied, the function will first fetch the names with the GBIF taxonomic
#'   backbone. In case it fails to find a match for one name, the function will
#'   return a list of non matched names. It is __strongly__ advised to use the
#'   taxon key.
#' @param start_year the initial year from which occurrences should be get
#' @param download_format the format to download data from GBIF. The default is
#'   "SIMPLE_PARQUET".
#' @param remove_zip \code{logical} indicating if the downloaded zip should be
#'   removed (\code{TRUE}).
#' @param savelog \code{logical} indicating if a log file containing the GBIF
#'   download DOI should be saved (\code{TRUE}).
#' @param mode mode of export. Either "single" (for single species) or "full"
#'   for multiple species. When using "full", a single file containing all the
#'   records will be saved in the "data/raw" folder
#' @param save_acro in case of full export, an acronym to be added to the saved
#'   file. If NULL, "full" is used.
#' @param geom_wkt a WKT geometry. Can be supplied together with
#'   \code{sci_names}, in which case data will be downloaded only for this area,
#'   or alone, in which case all occurrence records in a certain area will be
#'   downloaded. To use the geometry alone, \code{mode} should be \code{"full"},
#'   otherwise will fail.
#'
#' @return Files saved and, if not all species matched the taxonomic backbone, a
#'   \code{data.frame} containing the non matched names.
#' @export
#' 
#' @import arrow 
#' @import rgbif
#'
#' @examples
#' \dontrun{
#' mp_get_gbif("Leptuca thayeri", mode = "single")
#' }
mp_get_gbif <- function(
    sci_names = NULL,
    start_year = 1950,
    download_format = "SIMPLE_PARQUET",
    remove_zip = TRUE,
    savelog = TRUE,
    mode = "full",
    save_acro = NULL,
    geom_wkt = NULL
  ) {
  
  if (is.null(sci_names) & is.null(geom_wkt)) {
    stop("You should supply at least one of sci_names or geom_wkt")
  }
  
  if (is.null(sci_names) & mode != "full") {
    stop("When using geometry alone, you should use mode 'full'")
  }
  
  if (!is.null(sci_names)) {
    
    if (mode == "single" & length(sci_names) > 1) {
      stop("When downloading multiple species you should use mode = 'full'")
    }
    
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
  gbif_download <- occurrence_gbif(
    scientificname = sci_names,
    startdate = start_year,
    geometry = geom_wkt,
    absence = "include",
    exclude = c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"),
    format = download_format,
    wait = FALSE
  )
  
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


#' Retrieve occurrence records from local file for individual species and save
#' in the standard format
#'
#' @param local_file the relative or absolute path to the file (or folder in
#'   case of parquet datasets).
#' @param database_name the name of the database source, will be used to save the
#'   files.
#' @param sci_names the scientific name of the species. If \code{NULL}, then the
#'   aphia_id or gbif_key will be used instead (recommended). For all cases can
#'   also be a \code{vector} of values (see details).
#' @param aphia_id the AphiaID of the species. Can also be passed together with
#'   the GBIF key to ensure the file will be saved with the right AphiaID
#'   (otherwise it will be retrieved from WoRMS - see details).
#' @param gbif_key the GBIF taxon key.
#' @param save_format the format to save the file (either "parquet" or "csv";
#'   recommended is "parquet").
#' @param progress show progress bar (only if length of names > 1)
#' @param db_mode if \code{TRUE}, then a [duckdb::duckdb()] connection is used
#'  to query the file. This should be faster, specially when working with OBIS
#'  full export. However, for arrow datasets organized as a folder (e.g. GBIF
#'  downloads), using the standard mode is preferable.
#'  @param sel_columns enable to select just part of the available columns. This
#'  can considerably speed up the processing and reduce file size. Should be written
#'  as a single character string or a character vector.
#'
#' @return individual file saved
#' @export
#' 
#' @details
#' If you have a list of species, then it's possible to query through a vector of
#' names, AphiaID or taxonKey. In that case, for each species the function will
#' query the name from the database and save.
#' 
#' The individual files are saved in a folder according to the AphiaID. For the
#' cases in which a name is passed, the function will look in WoRMS for the
#' AphiaID of the species. In the case a GBIF taxonKey is passed, then the
#' function will first look for that name on the taxonomic backbone of GBIF
#' (through [rgbif::name_usage()]) and then search in WoRMS for the name. Thus, it's
#' advised that when using GBIF keys you also supply their AphiaID equivalent
#' (if you have it).
#' 
#' ## How Files are saved?
#' 
#' Files will be saved on "data/species" (if the folder does not exist on the
#' working directory it will be created) following a hive structure:
#' 
#' - key='AphiaID of the species'
#' - date='date being generated'
#' - ftype='type of file, i.e., which database generated'
#' 
#' So, for a species with AphiaID 12345, from both OBIS and GBIF,
#' downloaded on 2023-01-05, a folder would be structured as that:
#' 
#' data/species/key=12345/date=20230105
#' 
#' With the following folders inside:
#' ftype=obis/
#' ftype=gbif/
#' 
#' Each with a file called spdata.parquet (or csv).
#' 
#' This structure enables easy indexing and querying with multiple species, specially
#' when working with parquet datasets. Each of those keys will become a filtering feature.
#' For more details see [arrow::open_dataset()]
#'
#' @import arrow
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' mp_get_local("data/full_obis.parquet", "obis", "Leptuca thayeri")
#' }
mp_get_local <- function(local_file,
                         database_name,
                         sci_names = NULL,
                         aphia_id = NULL,
                         gbif_key = NULL,
                         save_format = "parquet",
                         progress = TRUE,
                         db_mode = FALSE,
                         sel_columns = NULL
                         ) {
  
  if (is.null(sci_names) & is.null(aphia_id) & is.null(gbif_key)) {
    stop("One of sci_names, aphia_id or gbif_key should be supplied.")
  }
  
  if (length(sci_names) < 2 & length(gbif_key) < 2 & length(aphia_id) < 2) {
    if (progress) {
      progress <- FALSE
    }
  }
  
  name_sup <- FALSE
  
  if (!is.null(aphia_id) & !is.null(gbif_key)) {
    cat("AphiaID and GBIF taxonKey supplied.
taxonKey will be used for search and AphiaID to save the file.")
    aphia_info <- aphia_id
    aphia_id <- NULL
    name_sup <- TRUE
  }
  
  if (!is.null(sci_names)) {
    nlen <- length(sci_names)
  }
  if (!is.null(aphia_id)) {
    nlen <- length(aphia_id)
  }
  if (!is.null(gbif_key)) {
    nlen <- length(gbif_key)
  }
  
  
  if (is.null(sel_columns)) {
    sel_columns <- "*"
    add_sel <- FALSE
  } else {
    sel_columns <- paste(sel_columns, collapse = ",")
    add_sel <- TRUE
  }
  
  
  if (db_mode) {
    if (fs::is_dir(local_file)) {
      ds_temp <- open_dataset(local_file)
      
      con <- DBI::dbConnect(duckdb::duckdb())
      
      duckdb::duckdb_register_arrow(con, "data_table", ds_temp)
      
      local_file <- "data_table"
      
    } else {
      con <- DBI::dbConnect(duckdb::duckdb())
    }
  } else {
    if (tools::file_ext(local_file) == "parquet") {
      database <- open_dataset(local_file) %>%
        {if(add_sel) select(., eval(strsplit(sel_columns, ",")[[1]])) else .}
    } else {
      database <- open_csv_dataset(local_file) %>%
        {if(add_sel) select(., eval(strsplit(sel_columns, ",")[[1]])) else .}
    }
  }
  
  if (progress) {
    cli::cli_progress_bar(
      format = paste0(
        "{cli::pb_spin} Getting species {z} of {nlen} ",
        "{cli::pb_bar} {cli::pb_percent} ETA:{cli::pb_eta}"
      ),
      format_done = paste0(
        "{cli::col_green(cli::symbol$tick)} Processed {nlen} species ",
        "in {cli::pb_elapsed}."
      ),
      total = nlen
    )
  }
  
  if (!any(save_format %in% c("parquet", "csv")) & length(save_format) != 1) {
    stop("Check save format.")
  }
  
  save_fun <- function(dat, save_format) {
    switch (save_format,
            parquet = write_parquet(dat, spath),
            csv = readr::write_csv(dat, spath)
    )
  }
  
  for (z in 1:nlen) {
    
    if (is.null(aphia_id)) {
      if (is.null(gbif_key)) {
        spnam <- worrms::wm_name2id(sci_name[z])
      } else {
        if (name_sup) {
          spnam <- aphia_info[z]
        } else {
          spnam <- name_usage(gbif_key[z])
          spnam <- spnam$data$scientificName
          spnam <- unlist(stringr::str_split(spnam, pattern = " "))
          spnam <- paste(spnam[1:2], collapse = " ")
          spnam <- worrms::wm_records_name(spnam)
          spnam <- spnam[spnam$status == "accepted", "AphiaID"][1,]
          spnam <- spnam$AphiaID
        }
      }
    } else {
      spnam <- aphia_id[z]
    }
    
    save_folder <- glue::glue("data/species/key={spnam}/date={format(Sys.Date(), '%Y%m%d')}/ftype={database_name}")
    fs::dir_create(save_folder)
    
    spath <- glue::glue("{save_folder}/spdata.{save_format}")
    
    if (is.null(aphia_id) & is.null(gbif_key)) {
      
      if (db_mode) {
        sp_data <- DBI::dbGetQuery(con,
                                   glue::glue(
                                     "SELECT {sel_columns}
            FROM '{local_file}'
            WHERE scientificName = '{sci_names[z]}';"
                                   ))
        save_fun(sp_data, save_format)
      } else {
        # Filter by scientific name
        database %>%
          filter(scientificName == sci_names[z]) %>%
          save_fun(save_format = save_format)
      }
    } else {
      if (is.null(aphia_id)) {
        
        if (db_mode) {
          sp_data <- DBI::dbGetQuery(con,
                                     glue::glue(
                                       "SELECT {sel_columns}
            FROM '{local_file}'
            WHERE taxonkey = {gbif_key[z]};"
                                     ))
          save_fun(sp_data, save_format)
        } else {
          # Filter by taxon key (GBIF)
          database %>%
            filter(taxonkey == gbif_key[z])  %>%
            save_fun(save_format = save_format)
        }
      } else {
        
        if (db_mode) {
          sp_data <- DBI::dbGetQuery(con,
                                     glue::glue(
                                       "SELECT {sel_columns}
            FROM '{local_file}'
            WHERE AphiaID = {aphia_id[z]};"
                                     ))
          save_fun(sp_data, save_format)
        } else {
          # Filter by AphiaID (OBIS)
          sp_data <- database %>%
            filter(AphiaID == aphia_id[z]) %>%
            save_fun(save_format = save_format)
        }
      }
    }
    
    if (progress) {
      cli::cli_progress_update()
    }
  }
  if (progress) {
    cli::cli_progress_done()
  }
  
  if (db_mode) {
    DBI::dbDisconnect(con, shutdown = T)
  }
  
  return(invisible(NULL))
  
}


# OLD VERSION
# mp_get_gbif <- function(
#     sci_names = NULL,
#     start_year = 1950,
#     download_format = "SIMPLE_PARQUET",
#     remove_zip = T,
#     savelog = T,
#     mode = "full",
#     area_only = FALSE,
#     save_acro = NULL,
#     geom_wkt = NULL
# ) {
#   
#   if (area_only & is.null(geom_wkt)) {
#     stop("For area only you SHOULD supply a WKT.")
#   }
#   if (area_only & mode != "full") {
#     stop("For area only, mode 'full' is the only currently supported.")
#   }
#   if (!area_only & is.null(sci_names)) {
#     stop("You should supply scientific names when not using area only.")
#   }
#   
#   if (!area_only) {
#     # Check if names were supplied or GBIF codes
#     if (is.character(sci_names)) {
#       gbif_taxon_keys <- name_backbone_checklist(sci_names, strict = T)
#       
#       if (sum(gbif_taxon_keys$matchType == "NONE") > 0) {
#         warning("GBIF failed to match names for ",
#                 sum(gbif_taxon_keys$matchType == "NONE"),
#                 " species. Check the output carefully.")
#         non_match <- sci_names[gbif_taxon_keys$matchType == "NONE"]
#         
#         gbif_taxon_keys <- gbif_taxon_keys[gbif_taxon_keys$matchType != "NONE",]
#       } else {
#         non_match <- NULL
#       }
#       
#       sel_keys <- gbif_taxon_keys$usageKey
#       
#     } else {
#       sel_keys <- sci_names
#       
#       non_match <- NULL
#     }
#   }
#   
#   # Set the download on GBIF server
#   if (!is.null(geom_wkt)) {
#     if (area_only) {
#       gbif_download <- occ_download(
#         pred_within(geom_wkt),
#         pred("hasCoordinate", TRUE),
#         pred_gte("year", start_year),
#         pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
#         format = download_format
#       )
#     } else {
#       gbif_download <- occ_download(
#         pred_in("taxonKey", sel_keys),
#         pred("hasCoordinate", TRUE),
#         pred_gte("year", start_year),
#         pred_within(geom_wkt),
#         pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
#         format = download_format
#       )
#     }
#   } else {
#     gbif_download <- occ_download(
#       pred_in("taxonKey", sel_keys),
#       pred("hasCoordinate", TRUE),
#       pred_gte("year", start_year),
#       pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
#       format = download_format
#     )
#   }
#   
#   # Check the download status
#   status <- occ_download_meta(gbif_download)$status
#   cli::cli_progress_bar("Waiting for download...")
#   while (status == "RUNNING" | status == "PREPARING") {
#     Sys.sleep(5)
#     cli::cli_progress_update()
#     status <- try(occ_download_meta(gbif_download)$status)
#     if (class(status) == "try-error") {
#       Sys.sleep(5)
#       status <- try(occ_download_meta(gbif_download)$status)
#     }
#   }
#   cli::cli_progress_update(force = T)
#   
#   if (status != "SUCCEEDED") {
#     stop("Downloaded failed.")
#   } else {
#     if (mode == "full") {
#       fs::dir_create("data/raw/")
#       dl <- occ_download_get(gbif_download, path = "data/raw")
#       
#       if (is.null(save_acro)) {
#         save_acro <- "full"
#       }
#       
#       unzip(paste0("data/raw/", gbif_download, ".zip"),
#             exdir = "data/raw/")
#       
#       switch (download_format,
#               SIMPLE_CSV = file.rename("data/raw/occurrence.csv", glue::glue("data/raw/gbif_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.csv")),
#               SIMPLE_PARQUET = file.rename("data/raw/occurrence.parquet", glue::glue("data/raw/gbif_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.parquet"))
#       )
#       
#       if (remove_zip) {
#         fs::file_delete(paste0("data/raw/", gbif_download, ".zip"))
#       }
#       
#       if (savelog) {
#         write.table(
#           data.frame(date = attr(gbif_download, "created"),
#                      download_code = as.character(gbif_download),
#                      doi = attr(gbif_download, "doi"),
#                      citation = attr(gbif_download, "citation")),
#           glue::glue("data/gbif_{save_acro}_download_log.txt"), row.names = F
#         )
#       }
#     } else {
#       
#       if (is.character(sci_names)) {
#         spnam <- sci_names
#       } else {
#         spnam <- name_usage(sci_names)
#         spnam <- spnam$data$scientificName
#       }
#       
#       spnam <- unlist(stringr::str_split(spnam, pattern = " "))
#       spnam <- paste(spnam[1:2], collapse = " ")
#       spnam <- worrms::wm_records_name(spnam)
#       spnam <- spnam[spnam$status == "accepted", "AphiaID"][1,]
#       spnam <- spnam$AphiaID
#       
#       sppath <- glue::glue("data/species/key={spnam}/date={format(Sys.Date(), '%Y%m%d')}")
#       
#       fs::dir_create(sppath)
#       dl <- occ_download_get(gbif_download, path = sppath)
#       
#       unzip(paste0(sppath, "/", gbif_download, ".zip"),
#             exdir = sppath)
#       
#       if (download_format == "SIMPLE_PARQUET") {
#         da <- open_dataset(paste0(sppath, "/occurrence.parquet"))
#       } else {
#         da <- data.table::fread(paste0(sppath, "/occurrence.csv"))
#       }
#       
#       write_parquet(da, paste0(sppath, "/type=gbif.parquet"))
#       rm(da)
#       
#       switch (download_format,
#               SIMPLE_PARQUET = fs::dir_delete(paste0(sppath, "/occurrence.parquet")),
#               SIMPLE_CSV = fs::dir_delete(paste0(sppath, "/occurrence.csv"))
#       )
#       
#       if (remove_zip) {
#         fs::file_delete(paste0(sppath, "/", gbif_download, ".zip"))
#       }
#       
#       if (savelog) {
#         write.table(
#           data.frame(date = attr(gbif_download, "created"),
#                      doi = attr(gbif_download, "doi")),
#           paste0(sppath, "/type=log_gbif.txt"), row.names = F
#         )
#       }
#     }
#   }
#   
#   if (!is.null(non_match)) {
#     return(non_match)
#   } else {
#     return(invisible(NULL))
#   }
#   
# }




### OLD VERSION MP GET LOCAL
# mp_get_local <- function(local_file,
#                          database_name,
#                          sci_names = NULL,
#                          aphia_id = NULL,
#                          gbif_key = NULL,
#                          save_format = "parquet",
#                          progress = TRUE,
#                          db_mode = FALSE
# ) {
#   
#   if (is.null(sci_names) & is.null(aphia_id) & is.null(gbif_key)) {
#     stop("One of sci_names, aphia_id or gbif_key should be supplied.")
#   }
#   
#   if (length(sci_names) < 2 & length(gbif_key) < 2 & length(aphia_id) < 2) {
#     if (progress) {
#       progress <- FALSE
#     }
#   }
#   
#   name_sup <- FALSE
#   
#   if (!is.null(aphia_id) & !is.null(gbif_key)) {
#     cat("AphiaID and GBIF taxonKey supplied.
# taxonKey will be used for search and AphiaID to save the file.")
#     aphia_info <- aphia_id
#     aphia_id <- NULL
#     name_sup <- TRUE
#   }
#   
#   if (!is.null(sci_names)) {
#     nlen <- length(sci_names)
#   }
#   if (!is.null(aphia_id)) {
#     nlen <- length(aphia_id)
#   }
#   if (!is.null(gbif_key)) {
#     nlen <- length(gbif_key)
#   }
#   
#   
#   if (db_mode) {
#     if (fs::is_dir(local_file)) {
#       ds_temp <- open_dataset(local_file)
#       
#       con <- DBI::dbConnect(duckdb::duckdb())
#       
#       duckdb::duckdb_register_arrow(con, "data_table", ds_temp)
#       
#       local_file <- "data_table"
#       
#     } else {
#       con <- DBI::dbConnect(duckdb::duckdb())
#     }
#   } else {
#     if (tools::file_ext(local_file) == "parquet") {
#       database <- open_dataset(local_file)
#     } else {
#       database <- open_csv_dataset(local_file)
#     }
#   }
#   
#   if (progress) {
#     cli::cli_progress_bar(
#       format = paste0(
#         "{cli::pb_spin} Getting species {z} of {nlen} ",
#         "{cli::pb_bar} {cli::pb_percent} ETA:{cli::pb_eta}"
#       ),
#       format_done = paste0(
#         "{cli::col_green(cli::symbol$tick)} Processed {nlen} species ",
#         "in {cli::pb_elapsed}."
#       ),
#       total = nlen
#     )
#   }
#   for (z in 1:nlen) {
#     if (is.null(aphia_id) & is.null(gbif_key)) {
#       
#       if (db_mode) {
#         sp_data <- DBI::dbGetQuery(con,
#                                    glue::glue(
#                                      "SELECT *
#             FROM '{local_file}'
#             WHERE scientificName = {sci_names[z]};"
#                                    ))
#       } else {
#         # Filter by scientific name
#         sp_data <- database %>%
#           filter(scientificName == sci_names[z]) %>%
#           collect()
#       }
#     } else {
#       if (is.null(aphia_id)) {
#         
#         if (db_mode) {
#           sp_data <- DBI::dbGetQuery(con,
#                                      glue::glue(
#                                        "SELECT *
#             FROM '{local_file}'
#             WHERE taxonkey = {gbif_key[z]};"
#                                      ))
#         } else {
#           # Filter by taxon key (GBIF)
#           sp_data <- database %>%
#             filter(taxonkey == gbif_key[z]) %>%
#             collect()
#         }
#       } else {
#         
#         if (db_mode) {
#           sp_data <- DBI::dbGetQuery(con,
#                                      glue::glue(
#                                        "SELECT *
#             FROM '{local_file}'
#             WHERE AphiaID = {aphia_id[z]};"
#                                      ))
#         } else {
#           # Filter by AphiaID (OBIS)
#           sp_data <- database %>%
#             filter(AphiaID == aphia_id[z]) %>%
#             collect()
#         }
#       }
#     }
#     
#     if (any(save_format %in% c("parquet", "csv")) & length(save_format) == 1) {
#       
#       if (is.null(aphia_id)) {
#         if (is.null(gbif_key)) {
#           spnam <- worrms::wm_name2id(sci_name[z])
#         } else {
#           if (name_sup) {
#             spnam <- aphia_info[z]
#           } else {
#             spnam <- name_usage(gbif_key[z])
#             spnam <- spnam$data$scientificName
#             spnam <- unlist(stringr::str_split(spnam, pattern = " "))
#             spnam <- paste(spnam[1:2], collapse = " ")
#             spnam <- worrms::wm_records_name(spnam)
#             spnam <- spnam[spnam$status == "accepted", "AphiaID"][1,]
#             spnam <- spnam$AphiaID
#           }
#         }
#       } else {
#         spnam <- aphia_id[z]
#       }
#       
#       save_folder <- glue::glue("data/species/key={spnam}/date={format(Sys.Date(), '%Y%m%d')}/ftype={database_name}")
#       fs::dir_create(save_folder)
#       
#       spath <- glue::glue("{save_folder}/spdata.{save_format}")
#       switch (save_format,
#               parquet = write_parquet(sp_data, spath),
#               csv = readr::write_csv(sp_data, spath)
#       )
#     } else {
#       stop("Check save format.")
#     }
#     
#     if (progress) {
#       cli::cli_progress_update()
#     }
#   }
#   if (progress) {
#     cli::cli_progress_done()
#   }
#   
#   if (db_mode) {
#     DBI::dbDisconnect(con, shutdown = T)
#   }
#   
#   return(invisible(NULL))
#   
# }
