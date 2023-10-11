#' Standardize occurrence data for SDMs
#'
#' @param species the AphiaID for the species
#' @param species_folder the folder where the species data is located (default
#'   is the package standard)
#' @param geo_out if \code{TRUE} assess geographical outliers
#' @param geo_out_mccore for computing geographic distances in parallel, number
#'   of cores
#' @param env_out if \code{TRUE} assess environmental outliers
#' @param prepare_sdm if \code{TRUE} prepare the data for SDM (see details)
#' @param skip_first_out if \code{TRUE} (default), it will perform the outlier checking
#' only for the SDM dataset, after converting to 1 per cell (see \code{out_from_cell})
#' @param sdm_base a base \code{SpatRaster} in the target resolution. You can
#'   supply one of the environmental layers that you will use (same extent and
#'   resolution planned for the SDM)
#' @param remove_land if \code{TRUE} remove points from land. If \code{sdm_base}
#'   is supplied, will use it as reference. Otherwise, will use the
#'   [obistools::check_onland()] function
#' @param out_from_cell only relevant when prepare_sdm is set \code{TRUE}. When this
#' option is \code{TRUE}, the function will first aggregate each dataset to 1 record per
#' cell. Then it will run the outlier procedure again based on the new rarefied set (see details)
#' @param narm_out if \code{TRUE} and \code{prepare_sdm = TRUE}, then the function
#' will not only remove the outliers, but also any NA value that there is (cells
#' for which no distance or environmental information was found)
#'
#' @return saved files (see details)
#' @export
#'
#' @details This function aims to standardize (i.e. prepare) the occurrence data
#' for use on the SDM. It does not have much flexibility in terms of arguments -
#' and this is in purpose, as the aim is to produce one particular configuration
#' for the SDMs. However, all the underlying functions are available for the user
#' in the package, so it can assembly its own 'standardize' function.
#' 
#' The function will perform the following steps:
#' 1. Open the most recent files from both OBIS and GBIF (those that are available)
#' 2. It will then verify if there are duplicated records, using GeoHash (resolution of 6)
#' and the year of the record.
#' 3. If \code{remove_land} is \code{TRUE}, will remove points from land
#' 4. If outlier removal methods were set as \code{TRUE}, it will tag (not remove)
#' the outliers (using limit of removal as 0.05%, see [outqc_geo()] or [outqc_env()] for more details).
#' Then, this file is saved, following this pattern:
#' 
#' species_folder/key=`SPECIES KEY`/date=`CURRENT DATE IN YMD FORMAT`/ftype=dupclean/spdata0.parquet
#' 
#' If you set prepare_sdm (there is no reason to not leave this as \code{TRUE}), then:
#' 
#' 1. The function will remove the outliers based on MAD (for environmental will
#' consider only bathymetry). If \code{out_from_cell = TRUE}, this will be done
#' based on a per-cell version of the dataset. This is recommended, because if you have
#' a very dense region, in terms of points, this can influence the outlier removal procedure.
#' Each dataset is converted to 1 per cell, then the equivalent cells in the aggregated
#' distance layer are retrieved (for the whole set of data), duplicated points are removed and outliers are detected. 
#' Any duplicated point (falling on the same point then other dataset) receives the same outlier status than its
#' duplicate.
#' 2. Then it will try to set one dataset to be used
#' as independent evaluation data. First it will calculate the standard distance
#' to assess how spread points are (spreadness in the rest of the doc)
#' and then the number of points in each dataset
#' (this considering the converted data to 1 per cell). After that it will
#' filter and remove those datasets with less than 10 records. If more than 1 dataset is
#' still available, it will order the datasets by their spreadness (high to low)
#' and remove the 50% with lowest spreadness. Then it will get the dataset with
#' the lowest impact in the total number of records to be used as an evaluation
#' dataset. 
#' 3. It will convert both the evaluation and the fitting dataset to 1 per cell and save.
#' 
#' It will save one single file following this pattern:
#' 
#' species_folder/key=`SPECIES KEY`/date=`CURRENT DATE IN YMD FORMAT`/ftype=stdpts/spdata0.parquet
#'
#' This file will have the following columns:
#' - decimalLongitude
#' - decimalLatitude
#' - data_type: one of "fit_points" and "eval_points". The last will only be available if the 
#'  function was able to retrieve an evaluation dataset. To use the data for SDMs
#'  you should then filter it and get a dataset for fiting the model (fit_points) and
#'  to evaluate the final model (eval_points).
#'  - dataset_sel: only for eval_points, the dataset that was used to retrieve the evaluation points
#'  - taxonID: the AphiaID, for control
#'  - species: scientificName, for control
#'  
#'  Aggregate resolution used in the outlier removal process is 0.8 (see [outqc_geo()])
#'  
#'  # Important note
#'  
#'  Even if the function succeed in setting one dataset apart for evaluation, note that the selected
#'  dataset may be important to model fitting! This is specially true if there was a small number of datasets.
#'  As explained, if after removing those datasets with less than 10 records there are at least 2 datasets,
#'  the function will proceed and select the 50% with highest spreadness for the next step. Thus, in this case
#'  the one dataset with high spreadness will be chosen and the one with lower spreadness will be used to
#'  model fitting.
#'  
#'  We recommend that you do one additional step of checking before using the dataset. For example, in this work
#'  we checked, before the SDM:
#'  1. If the number of records in the evaluation dataset was higher than in the fitting one, and;
#'  2. The inclusion of the evaluation dataset increased our coverage (number of TOTAL per cell records and
#'  spreadness).
#'  
#'  If this was the case, the evaluation dataset was incorporated in fitting (points were converted again
#'  to 1 per cell, as some points could overlap).
#'  
#' @examples
#' \dontrun{
#' mp_standardize(21510)
#' }
mp_standardize <- function(species,
                           species_folder = "data/species/",
                           geo_out = TRUE,
                           geo_out_mccore = NULL,
                           env_out = TRUE,
                           prepare_sdm = TRUE,
                           skip_first_out = TRUE,
                           sdm_base = NULL,
                           remove_land = TRUE,
                           out_from_cell = TRUE,
                           narm_out = TRUE
) {
  
  if (remove_land & is.null(sdm_base)) {
    warning("Points from land being removed using obistools.")
  }
  
  if (prepare_sdm & is.null(sdm_base)) {
    stop("To transform in 1 per cell, sdm_base should be supplied.")
  }
  
  lonlatdim <- c("decimalLongitude", "decimalLatitude")
  
  get_recent_file <- function(files) {
    dates <- gsub("\\/.*.", "", gsub("^(.*?)date=", "", files))
    dates <- as.Date(dates, "%Y%m%d")
    files[which.max(dates)]
  }
  
  gbif_file <- list.files(paste0(species_folder, "key=", species, "/"),
                          recursive = T)
  gbif_file <- get_recent_file(gbif_file[grepl("gbif", gbif_file)])
  
  obis_file <- list.files(paste0(species_folder, "key=", species, "/"),
                          recursive = T)
  obis_file <- get_recent_file(obis_file[grepl("obis", obis_file)])
  
  if (length(gbif_file) > 0) {
    gbif_data <- arrow::read_parquet(paste0(species_folder, "key=", species, "/", gbif_file))
    gbif_data <- dplyr::rename(gbif_data,
                               occurrenceID = occurrenceid,
                               catalogNumber = catalognumber)
  } else {
    gbif_data <- NULL
  }
  if (length(obis_file) > 0) {
    obis_data <- arrow::read_parquet(paste0(species_folder, "key=", species, "/", obis_file))
  } else {
    obis_data <- NULL
  }
  
  # Remove duplicates based on year and geohash
  if (!is.null(gbif_data) & !is.null(obis_data)) {
    non_dup_recs <- mp_dup_check(list(obis = obis_data, gbif = gbif_data))
  } else if (!is.null(gbif_data)) {
    non_dup_recs <- mp_dup_check(gbif_data)
  } else {
    non_dup_recs <- mp_dup_check(obis_data)
  }
  
  # Remove from land
  if (remove_land) {
    if (!is.null(sdm_base)) {
      pt_valid <- terra::extract(sdm_base, non_dup_recs[,lonlatdim])
      non_dup_recs <- non_dup_recs[!is.na(pt_valid[,2]),]
    } else {
      shoredistances <- obistools::lookup_xy(non_dup_recs[,lonlatdim], 
                                             shoredistance = TRUE, 
                            grids = FALSE, areas = FALSE, asdataframe = TRUE)
      pt_valid <- which(as.vector(shoredistances$shoredistance) >= 0)
      non_dup_recs <- non_dup_recs[pt_valid,]
    }
  }
  
  if (nrow(non_dup_recs) > 0) {
    # Perform geographical and environmental outlier removal
    if (!skip_first_out) {
      if (geo_out) {
        if (!is.null(sdm_base)) {
          base <- terra::mask(sdm_base, sdm_base, inverse = T, updatevalue = 1)
        } else {
          base <- outqc_get_base(0.05)
        }
        
        if (is.null(geo_out_mccore)) {
          geo_out_mccore <- 1
        }
        
        outqc_get_distances(base, target = non_dup_recs[, lonlatdim],
                            agg_res = 0.8, outfolder = "data/distances", skip_exists = T,
                            mc_cores = geo_out_mccore, verbose = F)
        
        non_dup_geotag <- outqc_geo(non_dup_recs[, lonlatdim],
                                    dist_folder = "data/distances",
                                    mc_cores = geo_out_mccore,
                                    limit_rem = 0.005)
        
        colnames(non_dup_geotag) <- paste0("GeoTag_", colnames(non_dup_geotag))
        
        # If you want to sum other methods and get an 'ensemble' like metric
        non_dup_geotag <- cbind(non_dup_geotag,
                                GeoTag_all = apply(non_dup_geotag, 1, function(x) sum(x)))
        
        non_dup_recs <- cbind(non_dup_recs, non_dup_geotag)
      }
      
      if (env_out) {
        non_dup_envtag <- outqc_env(non_dup_recs[, lonlatdim],
                                    limit_rem = 0.005)
        
        colnames(non_dup_envtag) <- paste0("EnvTag_", colnames(non_dup_envtag))
        
        non_dup_envtag <- cbind(non_dup_envtag,
                                EnvTag_all = apply(non_dup_envtag, 1, function(x) sum(x)))
        
        non_dup_recs <- cbind(non_dup_recs, non_dup_envtag)
      }
    }
    
    # Save duplicate records removal
    fs::dir_create(paste0(species_folder, "key=", species,
                          "/date=", format(Sys.Date(), "%Y%m%d"),
                          "/ftype=dupclean/"))
    
    arrow::write_parquet(non_dup_recs,
                         paste0(species_folder, "key=", species,
                                "/date=", format(Sys.Date(), "%Y%m%d"),
                                "/ftype=dupclean/spdata0.parquet"))
    
    if (prepare_sdm) {
      
      # Get species name
      species_name <- non_dup_recs$species[1]
      
      # See if there is any potential dataset to be left out for validation
      if (all(c("datasetkey", "datasetName", "datasetID") %in% colnames(non_dup_recs))) {
        non_dup_recs_filt <- non_dup_recs %>%
          mutate(unID = ifelse(is.na(datasetkey),
                               ifelse(is.na(datasetName), datasetID, datasetName),
                               datasetkey))
      } else if ("datasetkey" %in% colnames(non_dup_recs)) {
        non_dup_recs_filt <- non_dup_recs %>%
          mutate(unID = datasetkey)
      } else if ("datasetName" %in% colnames(non_dup_recs)) {
        non_dup_recs_filt <- non_dup_recs %>%
          mutate(unID = datasetName)
      } else {
        non_dup_recs_filt <- non_dup_recs %>%
          mutate(unID = datasetID)
      }
      
      non_dup_recs <- non_dup_recs_filt %>%
        mutate(unID = ifelse(is.na(unID), "NODATA", unID)) %>%
        filter(!is.na(unID))
      
      
      if (out_from_cell) {
        
        # First convert to 1 per cell (each dataset)
        per_cell <- function(dsID, dataset){
          ds <- dataset[dataset$unID == dsID,]
          
          cell_pts <- unique(terra::cellFromXY(sdm_base, as.matrix(ds[,lonlatdim])))
          
          rast_pts <- terra::xyFromCell(sdm_base, cell_pts)
          rast_pts <- data.frame(rast_pts)
          colnames(rast_pts) <- lonlatdim
          rast_pts <- cbind(rast_pts, unID = dsID)
          
          return(rast_pts)
        }
        
        non_dup_recs_cell <- lapply(unique(non_dup_recs$unID), per_cell, non_dup_recs)
        
        non_dup_recs <- do.call("rbind", non_dup_recs_cell)
        
        if (geo_out) {
          
          if (!is.null(sdm_base)) {
            base <- terra::mask(sdm_base, sdm_base, inverse = T, updatevalue = 1)
          } else {
            base <- outqc_get_base(0.05)
          }
          
          if (is.na(list.files("data/distances", pattern = ".tif", full.names = T)[1])) {
            outqc_get_distances(base, target = non_dup_recs[1:4, lonlatdim],
                                agg_res = 0.8, outfolder = "data/distances", skip_exists = T,
                                mc_cores = geo_out_mccore, verbose = F)
          }
          
          non_dup_recs$cell <- terra::cellFromXY(
            terra::rast(list.files("data/distances", pattern = ".tif", full.names = T)[1]),
            non_dup_recs[,lonlatdim]
          )
          
          non_dup_recs_sing <- non_dup_recs[!duplicated(non_dup_recs$cell),]
          
          if (is.null(geo_out_mccore)) {
            geo_out_mccore <- 1
          }
          
          outqc_get_distances(base, target = non_dup_recs_sing[, lonlatdim],
                              agg_res = 0.8, outfolder = "data/distances", skip_exists = T,
                              mc_cores = geo_out_mccore, verbose = F)
          
          non_dup_geotag <- outqc_geo(non_dup_recs_sing[, lonlatdim],
                                      dist_folder = "data/distances",
                                      mc_cores = geo_out_mccore,
                                      limit_rem = 0.005)
          
          colnames(non_dup_geotag) <- paste0("GeoTag_", colnames(non_dup_geotag))
          
          # If you want to sum other methods and get an 'ensemble' like metric
          non_dup_geotag <- cbind(non_dup_geotag,
                                  GeoTag_all = apply(non_dup_geotag, 1, function(x) sum(x)))
          
          non_dup_recs_sing <- cbind(non_dup_recs_sing, non_dup_geotag)
        }
        
        if (env_out) {
          non_dup_envtag <- outqc_env(non_dup_recs_sing[, lonlatdim],
                                      limit_rem = 0.005)
          
          colnames(non_dup_envtag) <- paste0("EnvTag_", colnames(non_dup_envtag))
          
          non_dup_envtag <- cbind(non_dup_envtag,
                                  EnvTag_all = apply(non_dup_envtag, 1, function(x) sum(x)))
          
          non_dup_recs_sing <- cbind(non_dup_recs_sing, non_dup_envtag)
        }
        
        if (geo_out | env_out) {
          non_dup_recs_sing <- non_dup_recs_sing[,4:ncol(non_dup_recs_sing)]
          non_dup_recs <- dplyr::left_join(non_dup_recs, non_dup_recs_sing, by = "cell")
          non_dup_recs <- non_dup_recs[,colnames(non_dup_recs) != "cell"]
        }
        
      }
      
      # Remove tagged outliers
      # Here we will remove based on MAD
      if (any(grepl("EnvTag", colnames(non_dup_recs)))) {
        non_dup_recs <- non_dup_recs %>%
          # Change here if you want to use other method (e.g. iqr)
          {if (!narm_out) filter(., EnvTag_bathymetry_mad != 1 | is.na(EnvTag_bathymetry_mad)) else filter(., EnvTag_bathymetry_mad != 1)}
        
      }
      if (any(grepl("GeoTag", colnames(non_dup_recs)))) {
        non_dup_recs <- non_dup_recs %>%
          # Change here if you want to use other method (e.g. iqr)
          {if (!narm_out) filter(., GeoTag_mdist_mad != 1 | is.na(GeoTag_mdist_mad)) else filter(., GeoTag_mdist_mad != 1)}
      }
      
      # Create a function that convert points to 1 per cell
      per_cell <- function(lonlat, data_type = NA, data_sel = NA){
        cell_pts <- unique(terra::cellFromXY(sdm_base, as.matrix(lonlat)))
        
        rast_pts <- terra::xyFromCell(sdm_base, cell_pts)
        rast_pts <- data.frame(rast_pts)
        rast_pts <- cbind(rast_pts, data_type = data_type, dataset_sel = data_sel)
        return(rast_pts)
      }
      
      if (nrow(non_dup_recs) > 0) {
        
        total_cell_recs <- nrow(per_cell(non_dup_recs[,lonlatdim]))
        
        get_spread <- function(long, lat){
          x <- cbind(long, lat)
          x_cords <- x
          x_cords <- apply(x_cords, 2, function(x) sum((x - mean(x))^2))
          sqrt(sum(x_cords)/nrow(x))
        }
        
        uniq_id <- unique(non_dup_recs$unID)
        
        non_dup_recs_sums <- lapply(uniq_id, function(x){
          pc <- per_cell(non_dup_recs[non_dup_recs$unID == x, lonlatdim])
          data.frame(
            unID = x,
            spreadness = get_spread(pc[,1], pc[,2]),
            n_recs = nrow(pc)
          )
        })
        
        non_dup_recs_sums <- dplyr::bind_rows(non_dup_recs_sums)
        
        non_dup_recs_sums <- non_dup_recs_sums %>%
          filter(n_recs > 10) %>%
          filter(unID != "NODATA") %>%
          # ungroup() %>%
          # filter(n_recs < max(n_recs)) %>%
          filter(!is.na(spreadness)) %>%
          mutate(impact = (n_recs / total_cell_recs) * 100)
        
        if (nrow(non_dup_recs_sums) > 1) {
          non_dup_recs_sums <- non_dup_recs_sums[order(non_dup_recs_sums$spreadness,
                                                       decreasing = T),]
          non_dup_recs_sums <- non_dup_recs_sums[1:(ceiling(nrow(non_dup_recs_sums)*.5)),]
          
          sel_dataid <- non_dup_recs_sums$unID[which.min(non_dup_recs_sums$impact)]
          
          hold_out_dataset <- non_dup_recs %>%
            filter(unID == sel_dataid)
          
          non_dup_recs <- non_dup_recs %>%
            filter(unID != sel_dataid)
        }
        
        # Convert in 1 per cell
        all_pts <- per_cell(as.matrix(non_dup_recs[,lonlatdim]),
                            data_type = "fit_points")
        
        if (exists("hold_out_dataset")) {
          all_pts <- rbind(all_pts,
                           per_cell(as.matrix(hold_out_dataset[,lonlatdim]),
                                    data_type = "eval_points"))
        }
        
        all_pts <- as.data.frame(all_pts)
        colnames(all_pts)[1:2] <- lonlatdim
        all_pts$taxonID <- species
        all_pts$species <- as.character(species_name)
        
        fs::dir_create(paste0(species_folder, "key=", species,
                              "/date=", format(Sys.Date(), "%Y%m%d"),
                              "/ftype=stdpts/"))
        
        arrow::write_parquet(all_pts,
                             paste0(species_folder, "key=", species,
                                    "/date=", format(Sys.Date(), "%Y%m%d"),
                                    "/ftype=stdpts/spdata0.parquet"))
      }
    }
  }
  
  return(invisible(NULL))
}
