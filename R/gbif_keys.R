############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################

#' Retrieve GBIF Taxon keys for species
#'
#' Get taxon keys used by GBIF based on species name or AphiaID.
#'
#' @param species a \code{character} vector containing species names or a
#'   \code{numeric} vector containing AphiaIDs (which will be converted to
#'   names) or a \code{data.frame} containing a column with species name (it
#'   should be named 'scientificName') and alternatively other columns with
#'   additional information (see details).
#' @param strict if \code{TRUE}, it matches only the given name and never an
#'   upper classification.
#' @param try_fuzzy if \code{TRUE}, it will also include names that the match
#'   was of the type "fuzzy" (see details on [rgbif::name_backbone()]).
#'   Otherwise (default), only exact matches are accepted.
#' @param try_names if \code{TRUE}, in case a name is not matched by including
#'   the full taxonomic name (i.e. species name + authorship), a second trial is
#'   done on the unmatched names using only species name.
#' @param verbose if \code{TRUE}, messages are printed.
#'
#' @return matched names with usage keys or a list with matched and unmatched
#'   names
#' @export
#'
#' @details The best way to ensure that the correct name is matched is to also
#' supply the authorship. In the case you supply only a vector of names, the
#' authorship can simply be pasted together.
#'
#' The best option, however, is to supply a data frame that contains a column
#' called 'scientificName' and other columns (named as in parenthesis) with
#' additional taxa information:
#' - authorship (scientificNameAuthorship) - strongly advised
#' - taxa rank (taxonRank)
#' - taxa phylum (phylum)
#' - taxa kingdom (kingdom)
#' - taxa class (class)
#' - taxa order (order)
#' - taxa family (family)
#'
#' Such a table is easily produced if you first get the names from WoRMS. So, if
#' you have a list of AphiaIDs, you can simply supply it to the
#' [worrms::wm_record()] and use the returned data frame in this function.
#' 
#' @examples
#' \dontrun{
#' get_gbif_keys("Leptuca thayeri")
#' }
get_gbif_keys <- function(species,
                          strict = TRUE,
                          try_fuzzy = FALSE,
                          try_names = FALSE,
                          verbose = TRUE) {
  
  if (try_fuzzy) {
    mtype <- c("EXACT", "FUZZY")
  } else {
    mtype <- "EXACT"
  }
  
  if (is.vector(species)) {
    
    if (!is.character(species)) {
      species <- worrms::wm_id2name(species)
    }
    
    taxa_name <- species
    
    taxa_rank <- taxa_kingdom <- taxa_phylum <- 
      taxa_class <- taxa_order <- taxa_family <- taxa_author <- NULL
    
  } else {
    if (!any(grepl("scientificname", colnames(species), ignore.case = T))) {
      stop("species should have a column named scientificName")
    }
    
    colnames(species) <- gsub("scientificname", "scientificName", colnames(species))
    
    taxa_name <- species[["scientificName"]]
    taxa_author <- species[["scientificNameAuthorship"]]
    taxa_rank <- species[["taxonRank"]]
    taxa_kingdom <- species[["kingdom"]]
    taxa_phylum <- species[["phylum"]]
    taxa_class <- species[["class"]]
    taxa_order <- species[["order"]]
    taxa_family <- species[["family"]]
    
  }
  
  if (is.null(taxa_author)) {
    taxa_full_name <- taxa_name
  } else {
    taxa_author <- ifelse(is.na(taxa_author), "", taxa_author)
    taxa_full_name <- paste(taxa_name, taxa_author)
    taxa_full_name <- trimws(taxa_full_name, "right")
  }
  
  # Checklist matching was returning json error
  # changing for name_backbone implementation
  # temporalily until the problem is solved
  # gbif_keys <- rgbif::name_backbone_checklist(
  #   name_data = taxa_full_name,
  #   rank = taxa_rank,
  #   kingdom = taxa_kingdom,
  #   phylum = taxa_phylum,
  #   class = taxa_class,
  #   order = taxa_order,
  #   family = taxa_family,
  #   strict = strict
  # )
  
  gbif_keys <- lapply(1:length(taxa_full_name), function(x){
    rgbif::name_backbone(
        name = taxa_full_name[x],
        rank = taxa_rank[x],
        kingdom = taxa_kingdom[x],
        phylum = taxa_phylum[x],
        class = taxa_class[x],
        order = taxa_order[x],
        family = taxa_family[x],
        strict = strict
    )
  })
  
  gbif_keys <- dplyr::bind_rows(gbif_keys)
  
  if (is.null(taxa_author)) {
    gbif_keys <- dplyr::bind_cols(obis_names = taxa_name,
                                  gbif_keys)
  } else {
    gbif_keys <- dplyr::bind_cols(obis_names = taxa_name,
                                  obis_author = taxa_author,
                                  gbif_keys)
  }
  
  gbif_keys$usageKey[!gbif_keys$matchType %in% mtype] <- NA
  
  non_match <- gbif_keys[is.na(gbif_keys$usageKey),]
  matched <- gbif_keys[!is.na(gbif_keys$usageKey),]
  
  if (nrow(non_match) > 0 & try_names) {
    
    if (verbose) cat("Missing match for some names. Trying with name only (without authorship).\n")
    
    matched$match_way <- "FULL"
    
    nm_gbif_keys <- lapply(1:nrow(non_match), function(x){
      rgbif::name_backbone(name = non_match$obis_names[x],
                           strict = strict)
    })
    
    nm_gbif_keys <- dplyr::bind_rows(nm_gbif_keys)
    
    if (is.null(taxa_author)) {
      nm_gbif_keys <- dplyr::bind_cols(obis_names = non_match$obis_names,
                                      nm_gbif_keys)
    } else {
      nm_gbif_keys <- dplyr::bind_cols(obis_names = non_match$obis_names,
                                      obis_author = non_match$obis_author,
                                      nm_gbif_keys)
    }
    
    if (!"usageKey" %in% colnames(nm_gbif_keys)) {
      nm_gbif_keys$usageKey <- NA
    }
    
    nm_gbif_keys$usageKey[!nm_gbif_keys$matchType %in% mtype] <- NA
    
    nm_gbif_keys$match_way <- "NAME_ONLY"
    
    non_match <- nm_gbif_keys[is.na(nm_gbif_keys$usageKey),]
    
    matched <- dplyr::bind_rows(
      matched,
      nm_gbif_keys[!is.na(nm_gbif_keys$usageKey),]
    )
    
  }
  
  if (nrow(non_match) == 0) {
    if (verbose) cat("All taxa matched")
    return(matched)
    
  } else {
    if (verbose) cat(paste("Match missing for", nrow(non_match),
                           "taxa, out of", nrow(matched), "- check output.\n"))
    return(list(
      match = matched,
      non_match = non_match
    ))
  }
  
}
# 
# teste <- get_gbif_keys(sp, try_fuzzy = T, try_names = T)
