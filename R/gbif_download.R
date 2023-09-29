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
#' @param username the GBIF username, if not stored on the R environment (see
#'   details).
#' @param password the GBIF password, if not stored on the R environment (see
#'   details).
#' @param exclude if \code{NULL} include all types of specimens. Alternatively,
#'   you can pass a vector of types of specimen to exclude (e.g. supplying
#'   \code{c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")} will exclude fossils and
#'   living specimens).
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
                            username = NULL,
                            password = NULL,
                            exclude = NULL,
                            format = NULL,
                            wait = TRUE,
                            raw_path = NULL,
                            import = TRUE,
                            verbose = FALSE) {
  
  # Check arguments
  args <- as.list(environment())
  args <- args[!grepl("verbose|username|password|format|wait", names(args))]
  
  if (all(is.null(unlist(args)))) {
    stop("You should supply at least one predicate.")
  }
  
  if (is.null(format)) {
    format <- "SIMPLE_CSV"
  }
  
  # Get GBIF user/password
  username <- .get_user_pass(username, "username")
  password <- .get_user_pass(password, "password")
  
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
