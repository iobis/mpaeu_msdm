#' Standardize occurrence data for SDMs
#'
#' @param species the AphiaID for the species
#' @param species_folder the folder where the species data is located (i.e. the
#' OBIS and GBIF full exports)
#' @param species_outfolder folder to output the files
#' @param geo_out if \code{TRUE} assess geographical outliers. Note that at least 5 records are necessary
#' for any outlier removal procedure (either geographical or environmental)
#' @param env_out if \code{TRUE} assess environmental outliers
#' @param sdm_base a base \code{SpatRaster} in the target resolution. You can
#'   supply one of the environmental layers that you will use (same extent and
#'   resolution planned for the SDM)
#' @param remove_land if \code{TRUE} remove points from land. If \code{sdm_base}
#'   is supplied, will use it as reference. Otherwise, will use the
#'   [obistools::check_onland()] function
#' @param narm_out if \code{TRUE}, then the function
#' will not only remove the outliers, but also any NA value that there is (cells
#' for which no distance or environmental information was found)
#' @param add_fao if `TRUE`, add information from FishBase and SeaLifeBase of where
#' the species is native, endemic or introduced, using FAO areas
#' @param fao_areas path to parquet file with FAO areas information (see [rfishbase::fb_tbl()])
#'
#' @return saved files (see details)
#' @export
#'
#' @details This function aims to standardize (i.e. prepare) the occurrence data
#' for use on the SDM. It does not have much flexibility in terms of arguments -
#' and this is in purpose, as the aim is to produce one particular configuration
#' for the SDMs. However, all the underlying functions are available for the user
#' in the package, so you can assembly your own 'standardize' function.
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
                           species_list,
                           sdm_base,
                           species_folder = "data/raw",
                           species_outfolder = "data/species",
                           geo_out = TRUE,
                           env_out = TRUE,
                           remove_land = TRUE,
                           narm_out = TRUE,
                           add_fao = TRUE,
                           fao_areas = "data/fao_areas.parquet",
                           verbose = TRUE
) {
  
  if (remove_land & is.null(sdm_base)) {
    warning("Points from land being removed using obistools.")
  }
  
  # if (prepare_sdm & is.null(sdm_base)) {
  #   stop("To transform in 1 per cell, sdm_base should be supplied.")
  # }
  
  lonlatdim <- c("decimalLongitude", "decimalLatitude")
  
  species_code <- species
  
  species_list <- read.csv(species_list)
  species_list <- species_list %>%
    filter(taxonID == species_code)
  
  if (nrow(species_list) < 1) {
    return(invisible(NULL))
  }
  
  require("polars")
  
  # Load data
  ds_files <- list.files(species_folder, full.names = T)
  obis_ds <- ds_files[grepl("obis_", ds_files)]
  gbif_ds <- ds_files[grepl("gbif_", ds_files)]
  
  obis_ds_pl <- pl$scan_parquet(obis_ds)
  gbif_ds_pl <- pl$scan_parquet(
    file.path(gbif_ds, "**/*")
  )
  
  gbif_columns <- c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                    "decimallongitude", "depth", "depthaccuracy", "day", "month",
                    "year", "taxonkey", "specieskey", "basisofrecord", "catalognumber")
  
  obis_columns <- c("occurrenceID", "catalogNumber", "recordNumber", "fieldNumber", "AphiaID",
                    "materialSampleID", "institutionID", "collectionID", "datasetID",
                    "collectionCode", "institutionCode", "datasetName",
                    "eventID", "parentEventID", "decimalLatitude", "decimalLongitude", 
                    "species", "eventDate", "date_year", "day", "month", "year",
                    "occurrenceStatus", "flags", "depth", "maximumDepthInMeters", 
                    "minimumDepthInMeters")
  
  obis_data <- obis_ds_pl$
    select(obis_columns)$
    filter(pl$col("AphiaID") == species_code)$
    collect()
  obis_data <- obis_data$to_data_frame()
  obis_data <- obis_data %>%
    rename(taxonID = AphiaID) %>%
    mutate(day = as.numeric(day),
           month = as.numeric(month),
           year = as.numeric(year))
  
  gbif_data <- gbif_ds_pl$
    select(gbif_columns)$
    filter(pl$col("specieskey") == species_list$gbif_speciesKey)$
    collect()
  gbif_data <- gbif_data$to_data_frame()
  gbif_data <- gbif_data %>%
    mutate(taxonID = species_list$taxonID) %>%
    rename(decimalLatitude = decimallatitude,
           decimalLongitude = decimallongitude)
  
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
    
    if (!is.null(sdm_base)) {
      base <- terra::mask(sdm_base, sdm_base, inverse = T, updatevalue = 1)
    } else {
      base <- outqc_get_base(0.05)
    }
    
    # if (is.null(geo_out_mccore)) {
    #   geo_out_mccore <- 1
    # }
    
    # if (is.na(list.files("data/distances", pattern = ".tif", full.names = T)[1])) {
    #   outqc_get_distances(base, target = non_dup_recs[1:4, lonlatdim],
    #                       agg_res = 0.8, outfolder = "data/distances", skip_exists = T,
    #                       mc_cores = geo_out_mccore, verbose = F)
    # }
    
    non_dup_recs$cell <- terra::cellFromXY(
      terra::rast(list.files("data/distances/distances.zarr/", pattern = ".tif", full.names = T)[1]),
      non_dup_recs[,lonlatdim]
    )
    
    non_dup_recs_sing <- non_dup_recs[!duplicated(non_dup_recs$cell),]
    
    if (geo_out) {
      
      if (nrow(non_dup_recs_sing) >= 5) {
        
        # outqc_get_distances(base, target = non_dup_recs_sing[, lonlatdim],
        #                     agg_res = 0.8, outfolder = "data/distances", skip_exists = T,
        #                     verbose = F)
        
        non_dup_geotag <- outqc_geo(non_dup_recs_sing[, lonlatdim],
                                    dist_folder = "data/distances",
                                    limit_rem = 0.01#0.005
                                    )
        
        colnames(non_dup_geotag) <- paste0("GeoTag_", colnames(non_dup_geotag))
        
        # If you want to sum other methods and get an 'ensemble' like metric
        non_dup_geotag <- cbind(non_dup_geotag,
                                GeoTag_all = apply(non_dup_geotag, 1, function(x) sum(x)))
        
        non_dup_recs_sing <- cbind(non_dup_recs_sing, non_dup_geotag)
      } else {
        warning(paste0("Less than 5 records (per cell) for the species ", species, ", impossible to perform outlier removal. \n"))
      }
      
    }
    
    if (env_out) {
      if (nrow(non_dup_recs_sing) >= 5) {
        non_dup_envtag <- outqc_env(non_dup_recs_sing[, lonlatdim],
                                    limit_rem = 0.01)
        
        colnames(non_dup_envtag) <- paste0("EnvTag_", colnames(non_dup_envtag))
        
        non_dup_envtag <- cbind(non_dup_envtag,
                                EnvTag_all = apply(non_dup_envtag, 1, function(x) sum(x)))
        
        non_dup_recs_sing <- cbind(non_dup_recs_sing, non_dup_envtag)
      }
    }
    
    if (geo_out | env_out) {
      if (nrow(non_dup_recs_sing) >= 5) {
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
      
      # Try to add tag reference to area of occurrence
      if (add_fao) {
        fao_areas <- pl$scan_parquet(fao_areas)
        
        sp_sbase <- try(rfishbase::species(list(all_pts$species[1])))
        
        if (!inherits(sp_sbase, "try-error")) {
          if (nrow(sp_sbase) < 1) {
            sp_sbase <- try(rfishbase::species(list(all_pts$species[1]),
                                               server = "sealifebase"))
          }
        }
        
        if (!inherits(sp_sbase, "try-error") & nrow(sp_sbase) > 0) {
          fao_areas_sp <- fao_areas$filter(pl$col("SpecCode") == sp_sbase$SpecCode[1])$
            collect()
          fao_areas_sp <- fao_areas_sp$to_data_frame()
          
          fao_areas_sp <- fao_areas_sp %>%
            filter(Status %in% c("endemic", "Endemic", "native", "Native", "introduced", "Introduced"))
          
          if (nrow(fao_areas_sp) > 0) {
            fao_shp <- vect("data/shapefiles/World_Fao_Zones.shp")
            
            fao_shp_sel <- fao_shp[fao_shp$zone %in% fao_areas_sp$AreaCode]
            
            if (length(fao_shp_sel) > 0) {
              sf::sf_use_s2(FALSE)
              fao_shp_sel <- suppressMessages(
                suppressWarnings(
                  sf::st_buffer(sf::st_make_valid(sf::st_as_sf(fao_shp_sel)), 0.2)
                )
              )
              
              is_onarea <- is.related(vect(all_pts[,1:2], geom = c("decimalLongitude", "decimalLatitude")),
                                      vect(fao_shp_sel),
                                      "intersects")
              is_onarea <- ifelse(is_onarea, 1, 0)
              
              all_pts$fao_confirmed <- is_onarea
            } else {
              all_pts$fao_confirmed  <- NA
            }
            
          } else {
            warning("Impossible to retrieve area information from FishBase/SeaLifeBase or no native areas")
          }
        } else {
          warning("Impossible to retrieve species information from FishBase/SeaLifeBase")
        }
      }
      
      if (!"fao_confirmed" %in% colnames(all_pts)) {
        all_pts$fao_confirmed <- NA
      }
      
      arrow::write_parquet(all_pts,
                           paste0(species_outfolder, "/key=", species, ".parquet"))
    }
  }
  
  return(invisible(NULL))
}
