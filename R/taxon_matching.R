#' Retrieve AphiaID for one or more scientific names
#'
#' This function is a wrapper around [worrms::wm_records_names()]. It will search
#' for a valid AphiaID based on the supplie scientific name
#'
#' @param scientificname a vector with one or more scientific names
#' @param marine_only only look for marine species
#' @param fuzzy try fuzzy matching for the names. Recomended is \code{FALSE}
#' @param try_valid the function will always get the accepted name. If none is available
#' and \code{try_valid=TRUE}, then it will check if a "valid_AphiaID" value is
#' present, and then retrieve the record for that valid id.
#' @param only_exact if \code{TRUE}, it will get only exact matches
#' @param progress show progress bar (only used if more than 100 scientific names are supplied)
#'
#' @return dataframe with valid AphiaIDs
#' @export
#'
#' @examples
#' \dontrun{
#'   name_to_aphia("Acanthurus chirurgus")
#' }
name_to_aphia <- function(scientificname,
                          marine_only = TRUE,
                          fuzzy = FALSE,
                          try_valid = TRUE,
                          only_exact = FALSE,
                          progress = TRUE) {
  
  if (length(scientificname) <= 100) {
    progress <- FALSE
  }
  
  blocks <- seq(1, length(scientificname), by = 100)
  if (blocks[length(blocks)] != length(scientificname)) {
    blocks <- c(blocks, length(scientificname))
  }
  
  if (length(blocks) > 1) {
    blocks <- blocks[-1]
  }
  
  st <- 1
  
  if (progress) cli::cli_progress_bar(total = length(blocks))
  
  for (z in 1:length(blocks)) {
    
    block_list <- scientificname[st:blocks[z]]
    
    worms_res <- try(worrms::wm_records_names(block_list,
                                              fuzzy = fuzzy,
                                              marine_only = marine_only),
                     silent = TRUE)
    
    if (inherits(worms_res, "try-error")) {
      marine <- FALSE
    } else {
      marine <- TRUE
    }
    
    if (marine) {
      
      worms_res <- lapply(1:length(worms_res), function(x){
        
        wmr <- worms_res[[x]]
        
        if (nrow(wmr) > 0 & only_exact) {
          wmr <- wmr[wmr$match_type == "exact"]
        }
        
        if (nrow(wmr) > 0) {
          if ("accepted" %in% wmr$status) {
            wmr <- wmr[wmr$status == "accepted",]
          } else if (try_valid) {
            if (!is.na(wmr$valid_AphiaID[1])) {
              wmr <- wm_record(wmr$valid_AphiaID[1])[1,]
            } else {
              return(wmr[0,])
            }
          } else {
            return(wmr[0,])
          }
          wmr$original_scientificname <- block_list[x]
        }
        return(wmr)
      })
      
      worms_res <- bind_rows(worms_res)
      
      if (!exists("marine_list")) {
        marine_list <- worms_res
      } else {
        marine_list <- bind_rows(marine_list, worms_res)
      }
    }
    
    st <- blocks[z] + 1
    if (progress) cli::cli_progress_update()
  }
  
  if (progress) cli::cli_progress_done()
  
  if (!exists("marine_list")) {
    cat("No AphiaID match found for any species.")
    return(invisible(NULL))
  }
  
  return(marine_list)
}




#' Retrieve AphiaID from GBIF taxonKey
#' 
#' Search WoRMS AphiaID for a GBIF taxonKey. The function first gets the scientific
#' name for the provided key(s) and then search for the AphiaID(s) using
#' [name_to_aphia()].
#'
#' @param taxonkey a vector of one or more GBIF taxonKeys
#' @param ... additional parameters passed to [name_to_aphia()]
#'
#' @return dataframe with names
#' @export
#'
#' @examples
#' \dontrun{
#'   key_to_aphia(2347567)
#' }
key_to_aphia <- function(taxonkey, ...) {
  
  key_name <- lapply(taxonkey, function(x){
    name_use <- try(rgbif::name_usage(x), silent = T)
    
    if (inherits(name_use, "try-error")) {
      return(NULL)
    } else {
      return(name_use$data)
    }
  })
  
  key_name <- bind_rows(key_name)
  
  if (nrow(key_name) < 1) {
    stop("No match found for any of the keys. Check and try again.")
  } else {
    aphia_ids <- name_to_aphia(key_name$species, ...)
    
    if (!is.null(aphia_ids)) {
      if (nrow(aphia_ids) > 0) {
        key_name <- key_name[,c("key", "species")]
        colnames(key_name) <- c("original_key", "original_scientificname")
        
        aphia_ids <- left_join(aphia_ids, key_name, by = "original_scientificname")
      }
    }
  }
  
  if (is.null(aphia_ids)) {
    return(invisible(NULL))
  } else {
    return(aphia_ids)
  }
}