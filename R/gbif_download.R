############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################

#' Download marine data from GBIF
#'
#' @param scientificname the scientific name of the species or the taxon key for
#'   the species obtained from GBIF. The second is preferred and you can obtain
#'   the taxon key using the function [get_gbif_keys()]. If a scientific name is
#'   supplied, the function will first look in the GBIF taxonomic backbone for a
#'   match and then use the key to download, but you will have no control on the
#'   matching.
#' @param taxonid the AphiaID of the species. If supplied, then the function
#'   will first look into WoRMS for the species name and then obtain the taxon
#'   key from GBIF. Ignored if \code{scientificname} is supplied.
#' @param startdate the earliest date on which occurrence took place.
#' @param enddate the latest date on which the occurrence took place.
#' @param startdepth 	the minimum depth below the sea surface.
#' @param enddepth the maximum depth below the sea surface.
#' @param geometry a WKT geometry string (should not be too complex, or the
#'   download will fail - see [rgbif::occ_download()]).
#' @param absence only include absence records (\code{TRUE}), exclude absence
#'   records (\code{NULL}) or include absence records (\code{include}).
#' @param exclude if \code{NULL} include all types of specimens. Alternatively,
#'   you can pass a vector of types of specimen to exclude (e.g. supplying
#'   \code{c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")} will exclude fossils and
#'   living specimens).
#' @param mode one of \code{search}, \code{cell} or \code{download}. Search will
#'   use the [rgbif::occ_search()] function, and have a limit of 200.000
#'   records. It is a good approach to fast results, but in general
#'   \code{download} should be used (see details). Mode \code{cell} uses the
#'   GBIF mapper API through the package [speedy] and returns cell records
#'   aggregated in a certain resolution. At this moment only \code{download} is
#'   available.
#' @param username the GBIF username, if not stored on the R environment (see
#'   details).
#' @param password the GBIF password, if not stored on the R environment (see
#'   details).
#' @param format format to download. One of \code{"SIMPLE_CSV"} or
#'   \code{"SIMPLE_PARQUET"}. If \code{NULL} (default), then a simple CSV is
#'   downloaded.
#' @param wait whether it should wait for the download to be completed. The
#'   download from GBIF occurs in a different way then from OBIS. First, a call
#'   is sent to the server, which will prepare the file for download. In
#'   general, this takes just a few minutes, but depending on the size it can
#'   take several minutes. If \code{wait} is \code{FALSE}, then the function
#'   will return the download call, which can later be used with
#'   [rgbif::occ_download_get())]. The default is \code{TRUE}, and the function
#'   will wait until the download is ready and download it on the supplied
#'   folder.
#' @param raw_path the path to the download folder. If \code{NULL}, the working
#'   folder is used. Ignored if \code{wait} is \code{FALSE}.
#' @param import if \code{TRUE}, the data is imported as a data frame.
#' @param cell_resolution if \code{mode = "cell"}, which resolution to use.
#' @param verbose if \code{TRUE} messages are printed.
#' 
#' @return a data frame containing occurrences or a GBIF download object.
#' @export
#'
#' @examples
#' \dontrun{
#' occurrence_gbif("Leptuca thayeri")
#' }
occurrence_gbif <- function(scientificname = NULL,
                            taxonid = NULL,
                            startdate = NULL,
                            enddate = NULL,
                            startdepth = NULL,
                            enddepth = NULL,
                            geometry = NULL,
                            absence = NULL,
                            exclude = NULL,
                            mode = "search",
                            username = NULL,
                            password = NULL,
                            format = NULL,
                            wait = TRUE,
                            raw_path = NULL,
                            import = TRUE,
                            cell_resolution = NULL,
                            verbose = FALSE) {
  
  require(rgbif)
  
  # Check arguments
  args <- as.list(environment())
  args <- args[!grepl("verbose|username|password|format|wait|mode|cell_resolution", names(args))]
  
  if (all(is.null(unlist(args)))) {
    stop("You should supply at least one predicate.")
  }
  
  if (is.null(format)) {
    format <- "SIMPLE_CSV"
  }
  
  if (mode == "download") {
    # Get GBIF user/password
    username <- .get_user_pass(username, "username")
    password <- .get_user_pass(password, "password")
  }
  
  # Get predicates
  pred_occ <- list(pred("hasCoordinate", TRUE))
  
  # Remove or include types of specimen
  if (!is.null(exclude)) {
    pred_occ <- c(pred_occ,
                  list(pred_not(pred_in("basisOfRecord", exclude))))
  }
  
  # Species names
  if (!is.null(scientificname)) {
    # Check if user supplied scientific names or the taxon keys for GBIF
    if (!is.numeric(scientificname)) {
      scientificname <- name_backbone_checklist(scientificname)
      if (nrow(scientificname) < 1) {
        stop("No match was found for the supplied name")
      }
      scientificname <- scientificname$usageKey
    }
    
    pred_occ <- c(pred_occ, list(pred_in("taxonKey", scientificname)))
  }
  
  # Species names from AphiaID
  if (!is.null(taxonid)) {
    if (!is.null(scientificname)) {
      if (verbose) cli::cli_warn("Scientific names and taxon ID supplied, ignoring taxon ID")
    } else {
      if (verbose) cli::cli_alert_info("Obtaining names from WoRMS")
      
      sp_names <- worrms::wm_id2name(taxonid)
      
      sp_names <- name_backbone_checklist(sp_names)
      
      sp_names <- sp_names$usageKey
      
      pred_occ <- c(pred_occ, list(pred_in("taxonKey", sp_names)))
    }
  }
  
  # Start date
  if (!is.null(startdate)) {
    if (nchar(startdate) > 4) {
      startdate <- format(as.Date(startdate), "%Y-%m-%d")
    }
    pred_occ <- c(pred_occ,
                  list(pred_gte("eventDate", startdate)))
  }
  
  # End date
  if (!is.null(enddate)) {
    if (nchar(enddate) > 4) {
      enddate <- format(as.Date(enddate), "%Y-%m-%d")
    }
    pred_occ <- c(pred_occ,
                  list(pred_lte("eventDate", enddate)))
  }
  
  # Start depth
  if (!is.null(startdepth)) {
    pred_occ <- c(pred_occ,
                  list(pred_gte("depth", startdepth)))
  }
  
  # End depth
  if (!is.null(enddepth)) {
    pred_occ <- c(pred_occ,
                  list(pred_lte("depth", enddepth)))
  }
  
  # WKT geometry
  if (!is.null(geometry)) {
    # Check if the geometry is too complex
    if (nchar(geometry) > 1500) {
      if (verbose) cli::cli_warn("Geometry seems overly complex. Download may fail.")
    }
    pred_occ <- c(pred_occ,
                  list(pred_within(geometry)))
  }
  
  # Absence status - following same pattern from robis
  if (!is.null(absence)) {
    if (isTRUE(absence)) {
      absence <- "ABSENT"
      
      pred_occ <- c(pred_occ,
                    list(pred("occurrenceStatus", absence)))
    }
  } else {
    absence <- "PRESENT"
    
    pred_occ <- c(pred_occ,
                  list(pred("occurrenceStatus", absence)))
  }
  
  # Generate download from GBIF
  gbif_download <- do.call(occ_download, c(format = format, pred_occ))
  
  # If wait, then wait until download is ready, download and open it.
  if (wait) {
    
    if (is.null(raw_path)) {
      raw_path <- getwd()
    }
    
    # Wait for download
    status <- occ_download_meta(gbif_download)$status
    if (verbose) cli::cli_progress_bar("Waiting for download...")
    while (status == "RUNNING" | status == "PREPARING") {
      Sys.sleep(5)
      if (verbose) cli::cli_progress_update()
      status <- try(occ_download_meta(gbif_download)$status)
      if (class(status) == "try-error") {
        Sys.sleep(5)
        status <- try(occ_download_meta(gbif_download)$status)
      }
    }
    if (verbose) cli::cli_progress_update(force = T)
    
    if (status != "SUCCEEDED") {
      stop("Downloaded failed.")
    } else {
      # Create folder to open the file
      fs::dir_create(raw_path)
      dl <- occ_download_get(gbif_download, path = raw_path)
      
      if (import) {
        call_res <- occ_download_import(dl)
      } else {
        call_res <- dl
      }
      
      return(call_res)
    }
    
  } else {
    return(gbif_download)
  }
  
}

#' @export
.get_user_pass <- function(value, what){
  
  if (is.null(value)) {
    if (what == "username") {
      what_b <- "GBIF_USER"
    } else {
      what_b <- "GBIF_PWD"
    }
    
    env_var <- Sys.getenv(what_b, "")
    
    if (env_var == "") {
      env_var <- getOption(tolower(what_b), NA)
      
      if (is.na(env_var)) {
        prompt <- paste("Enter your", what, " ")
        if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
          env_var <- rstudioapi::askForPassword(prompt)
        } else if (interactive()) {
          env_var <- readline(prompt)
        } else {
          stop(paste("You should supply a GBIF", what))
        }
      }
    }
    
    value <- env_var
  }
  
  return(value)
}


#' Download partial or full marine data from GBIF (using AWS export)
#'
#' @param full_mode if `TRUE` (default) then the full export of GBIF is downloaded.
#'  That will consume something around 250GB of space. Note that in this case,
#'  `scientificname`, `taxonid`, `startdate` or/and `enddate` are ignored. 
#'  If `FALSE`, then a partial download is executed based on the filters (one of
#'   `scientificname`, `taxonid`, `startdate` or/and `enddate` needs to be supplied)
#' @param export_path path to a folder where file will be saved
#' @param scientificname the scientific name of the species or the taxon key for
#'   the species obtained from GBIF. The second is preferred and you can obtain
#'   the taxon key using the function [get_gbif_keys()]. If a scientific name is
#'   supplied, the function will first look in the GBIF taxonomic backbone for a
#'   match and then use the key to download, but you will have no control on the
#'   matching.
#' @param taxonid the AphiaID of the species. If supplied, then the function
#'   will first look into WoRMS for the species name and then obtain the taxon
#'   key from GBIF. Ignored if \code{scientificname} is supplied.
#' @param startdate the earliest date on which occurrence took place.
#' @param enddate the latest date on which the occurrence took place.
#' @param format format to download in case of partial download. One of \code{"parquet"} or
#'   \code{"csv"}. If \code{NULL} (default), then \code{"parquet"} is
#'   downloaded.
#' @param backend the backend to be used in case of partial download. One of `arrow` or `duckdb` (the later is working better).
#' @param version the version of the full export to use (see: https://registry.opendata.aws/gbif/). 
#'   If `NULL`, the latest one is used
#' @param bucket which bucket to use. If `NULL` the default bucket is used (see [gbifdb::gbif_remote()])
#' @param verbose if \code{TRUE} messages are printed.
#' 
#' @return path to saved file (and files saved)
#' @export
#' 
#' @details
#' Downloading the full dataset is much faster if your query is really big (i.e.
#' the resulting data have many millions of records), but with the cost of a lot of storage
#' consumed. Of course, you can later process and reduce the file.
#' 
#' When doing partial download, files are saved partitioned by the taxonomic order.
#' 
#' DuckDB is performing better and is the recommended backend in case of not full download.
#' 
#'
#' @examples
#' \dontrun{
#' dir <- tempdir()
#' occurrence_gbif_db("Leptuca thayeri", dir)
#' }
occurrence_gbif_db <- function(full_mode = TRUE,
                               export_path,
                               scientificname = NULL,
                               taxonid = NULL,
                               startdate = NULL,
                               enddate = NULL,
                               format = NULL,
                               backend = "duckdb",
                               version = NULL,
                               bucket = NULL,
                               verbose = FALSE) {
  
  fs::dir_create(export_path)
  
  if (!full_mode) {
    
    if (verbose) cli::cli_alert_info("Partial download. Initializing.")
    
    # Check arguments
    args <- as.list(environment())
    args <- args[grepl("scientificname|taxonid|startdate|enddate", names(args))]
    
    if (all(is.null(unlist(args)))) {
      stop("You should supply at least one predicate.")
    }
    
    require(dplyr)
    
    if (is.null(format)) {
      format <- "parquet"
    } else if (format != "parquet" & format != "csv") {
      warning("`format` should be one of 'parquet' or 'csv'. Using 'parquet'")
      format <- "parquet"
    }
    
    if (is.null(version)) {
      version <- gbifdb::gbif_version()
    }
    if (is.null(bucket)) {
      bucket <- gbifdb:::gbif_default_bucket()
    }
    
    bucket_add <- paste0("s3://", bucket, "/occurrence/", version, "/occurrence.parquet/")
    
    if (!is.null(scientificname) & !is.null(taxonid)) {
      cli::cli_warn("Only one of `scientificname` or `taxonid` should be supplied. Using `taxonid` instead")
      scientificname <- NULL
    }
    
    # Species names
    if (!is.null(scientificname)) {
      # Check if user supplied scientific names or the taxon keys for GBIF
      if (!is.numeric(scientificname)) {
        scientificname <- rgbif::name_backbone_checklist(scientificname)
        if (nrow(scientificname) < 1) {
          stop("No match was found for the supplied name")
        }
        tkey <- scientificname$usageKey
      } else {
        tkey <- scientificname
      }
    }
    
    # Species names from AphiaID
    if (!is.null(taxonid)) {
      if (verbose) cli::cli_alert_info("Obtaining names from WoRMS")
      
      sp_names <- worrms::wm_id2name(taxonid)
      
      sp_names <- rgbif::name_backbone_checklist(sp_names)
      
      tkey <- sp_names$usageKey
    }
    
    # Start date
    if (!is.null(startdate)) {
      if (nchar(startdate) > 4) {
        startdate <- format(as.Date(startdate), "%Y-%m-%d")
      }
    }
    
    # End date
    if (!is.null(enddate)) {
      if (nchar(enddate) > 4) {
        enddate <- format(as.Date(enddate), "%Y-%m-%d")
      }
    }
    
    if (backend == "arrow") {
      
      if (verbose) cli::cli_alert_info("Downloading through arrow backend.")
      
      export_path <- paste0(export_path, "/occurrence_gbif_", format(Sys.Date(), "%Y%m%d"), ".", format)
      
      ds <- arrow::open_dataset(bucket_add)
      
      ds <- ds %>%
        select(-identifiedby, -recordedby, -typestatus, -mediatype, -issue)
      
      if (!is.null(scientificname) | !is.null(taxonid)) {
        ds <- ds %>%
          filter(taxonkey %in% tkey)
      } 
      if (!is.null(startdate)) {
        ds <- ds %>%
          filter(year >= lubridate::year(startdate))
      }
      if (!is.null(enddate)) {
        ds <- ds %>%
          filter(year <= lubridate::year(enddate))
      }
      
      ds %>%
        arrow::write_parquet(sink = export_path)
      
    } else if (backend == "duckdb") {
      
      if (verbose) cli::cli_alert_info("Downloading through DuckDB backend.")
      
      require(DBI)
      require(duckdb)
      
      con <- dbConnect(duckdb())
      dbSendQuery(con, "install httpfs; load httpfs;")
      
      query_call <- "where "
      
      if (!is.null(scientificname) | !is.null(taxonid)) {
        query_call <- paste0(query_call, "taxonkey in (", paste0(tkey, collapse = ", "), ")")
        if (!is.null(startdate) | !is.null(enddate)) {
          query_call <- paste0(query_call, " and ")
        }
      } 
      if (!is.null(startdate) & !is.null(enddate)) {
        query_call <- paste0(query_call, "year >= ", lubridate::year(startdate), " and year <= ", lubridate::year(enddate))
      } else if (!is.null(startdate)) {
        query_call <- paste0(query_call, "year >= ", lubridate::year(startdate))
      } else if (!is.null(enddate)) {
        query_call <- paste0(query_call, "year <= ", lubridate::year(enddate))
      }
      
      dbSendQuery(con, glue::glue(
        "
        copy (
          select * exclude (mediatype, issue, identifiedby, recordedby, typestatus)
          from read_parquet('{bucket_add}*')
          {query_call}
        ) to '{export_path}' (format {format}, partition_by ('order'), overwrite_or_ignore)
       "
      ))
      
      DBI::dbDisconnect(con)
      
    } else {
      stop("backend not recognized - should be one of 'parquet' or 'duckdb'")
    }
    
    if (verbose) cli::cli_alert_success("File downloaded and available at {.path {export_path}}")
    
    return(export_path)
    
  } else {
    
    if (verbose) cli::cli_alert_info("Downloading full GBIF export (aws).")
    
    fs::dir_create(paste0(export_path, "/occurrence/"))
    
    if (is.null(version)) {
      version <- gbifdb::gbif_version()
    }
    if (is.null(bucket)) {
      bucket <- gbifdb:::gbif_default_bucket()
    }
    
    download_data <- function() {
      gbifdb::gbif_download(
        version = version, 
        dir = export_path,
        bucket = bucket
      )
    }
    
    base_url <- Sys.getenv("AWS_S3_ENDPOINT", "s3.amazonaws.com")
    if (getOption("gbif_unset_aws", TRUE)) gbifdb:::unset_aws_env()
    minioclient::mc_alias_set("aws", base_url, "", "")
    estimated_total <- minioclient::mc_du(paste0("aws/", bucket, "/occurrence/", version, "/occurrence.parquet/"))
    estimated_total <- as.numeric(gsub("GiB*.*", "", estimated_total$stdout))
    
    if (verbose) cli::cli_alert_warning("Estimated size of the file: {estimated_total}G. Abort if you think there will be no space.")
    
    check_folder_size <- function(folder_path) {
      folder_size <- fs::dir_info(folder_path, recurse = T)
      folder_size <- dplyr::filter(folder_size, type == "file")
      folder_size <- dplyr::summarize(folder_size, size = sum(size))
      folder_size$size
    }
    
    future::plan(future::multisession)
    download_future <- future::future(download_data())
    
    tstart <- Sys.time()
    if (verbose) cli::cli_progress_bar(format = 
        "{cli::pb_spin} Downloading GBIF full export. {as.character(folder_size)} out of {estimated_total}G. ETA: {eta}")
    while (!future::resolved(download_future)) {
      folder_size <- check_folder_size(paste0(export_path, "/occurrence/"))
      
      remaining <- fs::as_fs_bytes(paste0(estimated_total, "G")) - folder_size
      elapsed_time <- as.numeric(difftime(Sys.time(), tstart, units = "mins"))
      download_speed_mb_per_min <- folder_size / elapsed_time
      
      eta <- remaining / download_speed_mb_per_min
      
      tdif <- tstart - Sys.time()
      
      eta <- ifelse(eta <= 60, paste(round(eta, 1), "minutes"), paste(round(eta/60, 1), "hours"))
      
      Sys.sleep(2)
      if (verbose) cli::cli_progress_update()
    }
    
    if (verbose) cli::cli_progress_update(force = T)
    if (verbose) cli::cli_progress_done()
    
    if (verbose) cli::cli_alert_success("File downloaded and available at {.path {paste0(export_path, '/occurrence/')}}")
    
    return(paste0(export_path, "/occurrence/", version, "/occurrence.parquet/"))
  }
  
}
