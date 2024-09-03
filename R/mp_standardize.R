#' Standardize occurrence data for SDMs
#'
#' @param species the AphiaID for the species
#' @param species_list the list of species, with species taxonID/AphiaID and
#'   GBIF keys.
#' @param sdm_base a base \code{SpatRaster} in the target resolution. You can
#'   supply one of the environmental layers that you will use (same extent and
#'   resolution planned for the SDM). Supplying this is highly recommended
#' @param species_folder the folder where the species data is located (i.e. the
#'   OBIS and GBIF full exports)
#' @param reader which package to use for reading the parquet databases. One of
#'   "polars", "duckdb" or "arrow". If you followed the project steps, then "polars" is
#'   recommended because is faster. But if you are facing any problems, just
#'   change for "arrow" or "duckdb" (the last have better memory usage)
#' @param species_outfolder folder to output the files
#' @param geo_out if \code{TRUE} assess geographical outliers. See
#'   [outqc_geo()]. Note that at least 5 records are necessary for any outlier
#'   removal procedure (either geographical or environmental)
#' @param geo_method which method to use for geographical outliers. One of
#'   "iqr", "mad", "isoforest" or "ensemble"
#' @param geo_dists_folder the folder with the distances for geographical
#'   outliers. See note on [outqc_geo()]
#' @param env_out if \code{TRUE} assess environmental outliers. See
#'   [outqc_env()]
#' @param env_variable which variable to use for environmental outliers. One of
#'   "bathymetry", "temperature", "salinity" or "shoredistance"
#' @param env_method which method to use for environmental outliers. One of
#'   "iqr", "mad", "isoforest" or "ensemble". Note that for both isoforest and
#'   ensemble, \code{env_variable} is ignored and all are used.
#' @param out_threshold the limit of points that can be removed (i.e.
#'   percentage). Should be a value between 0 and 1. The default is 0.01, which
#'   limits the maximum number of points that can be removed to 1% (as we assume
#'   outliers are extremes).
#' @param remove_land if \code{TRUE} remove points from land. If \code{sdm_base}
#'   is supplied, will use it as reference. Otherwise, will use the
#'   [obistools::check_onland()] function
#' @param approximate_land if \code{TRUE} and \code{sdm_base} is supplied, the
#'   function will attempt to approximate points on land to the nearest valid
#'   cell, based on approximate limit. This function uses adjacent cells, so it
#'   is not exactly distance based (as distances may vary according to latitude)
#' @param approximate_limit either 1 or 2, the number of adjacent cells to any
#'   direction to try when moving the point on land. 1 uses 'queen' movement,
#'   and 2 uses a 5x5 matrix (2 cells for each side/diagonal). See
#'   [terra::adjacent()] for more information
#' @param narm_out if \code{TRUE}, then the function will not only remove the
#'   outliers, but also any NA value that there is (cells for which no distance
#'   or environmental information was found)
#' @param eval_max_impact the function attempts to separate one dataset for
#'   evaluation. It will measure the spread of points (see details) and the
#'   impact of the dataset removal on the full dataset (in terms of percentage
#'   of points from the total). You can limit the maximum impact allowed when
#'   chosing a dataset
#' @param min_year minimum year to filter the data
#' @param exclude_record_type which types of observation ("basisofRecord) are 
#'   not permited (i.e. should be excluded).
#'   Should be a character vector with valid values for both OBIS and GBIF (see
#'   [rgbif::occ_download()] and [robis::occurrence()]). If NULL, it will use
#'   the default values (see details). For ignoring this parameter, set as `NA`
#' @param add_fao if `TRUE`, add information from FishBase and SeaLifeBase of
#'   where the species is native, endemic or introduced, using FAO areas
#' @param fao_areas path to parquet file with FAO areas information (see
#'   [rfishbase::fb_tbl()])
#' @param verbose if \code{TRUE}, print messages
#'
#' @return saved files (see details)
#' @export
#'
#' @details This function aims to standardize (i.e. prepare) the occurrence data
#'   for use on the SDM. Note that all the underlying functions are available
#'   for the user in the package, so you can assembly your own 'standardize'
#'   function if you need something not supported by the parameters available.
#'   The default values are the ones used on the MPA Europe project.
#'
#'   The function will perform the following steps: 1. Open the most recent
#'   files from both OBIS and GBIF (those that are available) 2. It will then
#'   verify if there are duplicated records, using GeoHash (resolution of 6) and
#'   the year of the record. 3. If \code{remove_land} is \code{TRUE}, will
#'   remove points from land. If \code{approximate_land} is \code{TRUE}, it will
#'   try to bring points on land to the closest valid cell. This is limited by
#'   \code{approximate_limit}, so that points that are above this limit are
#'   discarded. 4. If outlier removal methods were set as \code{TRUE}, it will
#'   remove the outliers according to the methods/thresholds 5. The function
#'   will try to separate one dataset for evaluation. It does that by first
#'   converting points to 1 per cell, per dataset. Then, it measure the "spread"
#'   of the points on the geographical space by applying the follow:
#' \code{
#' x_cords <- apply(x_cords, 2, function(x) sum((x - mean(x))^2))
#' sqrt(sum(x_cords)/nrow(x))
#' }
#'   Being x_cords the coordinates. Thus, the spread of points is measured by
#'   calculating the root mean square deviation of their coordinates from their
#'   respective means. It will then order the datasets from those with higher to
#'   lower spread. It will select the first 50% with highest spread. From those,
#'   it will select the one with lowest impact on the number of points ( with
#'   impact measured as the percentage from the total number which the dataset
#'   corresponds). To avoid removing a dataset with high impact on the total
#'   number of points, you can limit it through the argument
#'   \code{eval_max_impact}. If a dataset is available to be used as evaluation,
#'   those are marked as "eval_points", while those for fitting are marked as
#'   "fit_points".
#'
#'   The returned parquet file will contain:
#' - decimalLongitude
#' - decimalLatitude
#' - data_type (either fit or evaluation)
#' - dataset_sel (datasetID of the dataset used for evaluation, if available)
#' - taxonID
#' - species (scientific name)
#'   and if \code{add_fao = TRUE}, a column named "fao_confirmed" that contains
#'   1 or 0, with 1 meaning that the point is within the FAO area considered as
#'   native on FishBase or SeaLifeBase
#'   
#'   Basis of record default
#'   If nothing is supplied to `exclude_record_type`, then the function will 
#'   filter to remove occurrence records from the following types:
#'   - "FOSSIL_SPECIMEN" or "FossilSpecimen"
#'   - "LIVING_SPECIMEN" or "LivingSpecimen"
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' base_rast <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")
#' mp_standardize(21510, "data/all_splist_20240319.csv", base_rast)
#' }
mp_standardize <- function(species,
                           species_list,
                           sdm_base,
                           species_folder = "data/raw",
                           reader = "polars",
                           species_outfolder = "data/species",
                           geo_out = TRUE,
                           geo_method = "isoforest",
                           geo_dists_folder = "data/distances",
                           env_out = TRUE,
                           env_variable = "bathymetry",
                           env_method = "mad",
                           out_threshold = 0.01,
                           remove_land = TRUE,
                           approximate_land = TRUE,
                           approximate_limit = 1,
                           narm_out = TRUE,
                           eval_max_impact = 10,
                           min_year = 1950,
                           exclude_record_type = NULL,
                           add_fao = TRUE,
                           fao_areas = "data/fao_areas.parquet",
                           verbose = TRUE
) {
  
  if (!reader %in% c("polars", "arrow", "duckdb")) {
    stop("Reader should be one of `polars`, `arrow` or `duckdb`")
  }
  
  if (remove_land & is.null(sdm_base)) {
    warning("Points from land being removed using obistools.")
  }
  
  if (approximate_land & is.null(sdm_base)) {
    warning("Approximation of points in land only possible with `sdm_base`. Ignoring.")
    approximate_land <- FALSE
  }
  
  if (approximate_limit > 2) {
    warning("Maximum approximation is 2, using 2 instead.")
    approximate_limit <- 2
  } else if (is.integer(approximate_limit)) {
    approximate_limit <- as.integer(approximate_limit)
  }
  
  if (!geo_method %in% c("iqr", "mad", "isoforest", "ensemble")) {
    warning("Invalid geo_method, using IQR.")
    geo_method <- "iqr"
  }
  if (!env_method %in% c("iqr", "mad", "isoforest", "ensemble")) {
    warning("Invalid env_method, using IQR.")
    env_method <- "iqr"
  }
  if (!env_variable %in% c("bathymetry", "temperature", "salinity", "shoredistance")) {
    warning("Invalid env_variable, using bathymetry")
    env_variable <- "bathymetry"
  }
  
  if (is.null(exclude_record_type)) {
    basisrec <- c("FossilSpecimen", "FOSSIL_SPECIMEN", "LivingSpecimen", "LIVING_SPECIMEN")
  } else if (is.na(exclude_record_type)) {
    basisrec <- "excludenothing"
  }
  
  fs::dir_create(species_outfolder)
  
  lonlatdim <- c("decimalLongitude", "decimalLatitude")
  
  species_code <- species
  
  species_list <- read.csv(species_list)
  species_list <- species_list %>%
    filter(taxonID == species_code)
  
  if (nrow(species_list) < 1) {
    return(invisible(NULL))
  }
  
  if (reader == "polars") {
    require("polars")
  } else if (reader == "arrow") {
    require("arrow")
  } else if (reader == "duckdb") {
    require("DBI")
    require("duckdb")
  }
  
  # Load data
  ds_files <- list.files(species_folder, full.names = T)
  obis_ds <- ds_files[grepl("obis_", ds_files)]
  gbif_ds <- ds_files[grepl("gbif_", ds_files)]
  
  if (length(obis_ds) < 1 & length(gbif_ds) < 1) {
    stop("At least one of OBIS or GBIF data should be available.")
  }
  
  gbif_columns <- c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                    "decimallongitude", "depth", "depthaccuracy", "day", "month",
                    "year", "taxonkey", "specieskey", "basisofrecord", "catalognumber", "occurrencestatus")
  
  obis_columns <- c("occurrenceID", "catalogNumber", "recordNumber", "fieldNumber", "AphiaID",
                    "materialSampleID", "institutionID", "collectionID", "datasetID",
                    "collectionCode", "institutionCode", "datasetName",
                    "eventID", "parentEventID", "decimalLatitude", "decimalLongitude", 
                    "species", "eventDate", "date_year", "day", "month", "year",
                    "occurrenceStatus", "flags", "depth", "maximumDepthInMeters", 
                    "minimumDepthInMeters", "absence", "basisOfRecord")
  
  if (length(obis_ds) > 0) {
    
    if (reader == "polars") {
      obis_ds_pl <- pl$scan_parquet(obis_ds)
      
      obis_data <- obis_ds_pl$
        select(obis_columns)$
        filter(pl$col("AphiaID") == species_code)$
        collect()
      
      obis_data <- obis_data$to_data_frame()
      rm(obis_ds_pl)
    } else if (reader == "arrow") {
      obis_ds_pl <- open_dataset(obis_ds)
      
      obis_data <- obis_ds_pl %>%
        select(all_of(obis_columns)) %>%
        filter(AphiaID == species_code) %>%
        collect()
      rm(obis_ds_pl)
    } else if (reader == "duckdb") {
      con <- dbConnect(duckdb())
      obis_data <- dbGetQuery(con, glue::glue(
        "
        select {paste(obis_columns, collapse = ', ')}
        from read_parquet('{obis_ds}')
        where AphiaID = {species_code}
       "
      ))
      DBI::dbDisconnect(con)
    }
    
    obis_data <- obis_data %>%
      rename(taxonID = AphiaID) %>%
      mutate(day = as.numeric(day),
             month = as.numeric(month),
             year = as.numeric(year)) %>%
      filter(year >= min_year)
    
    if (basisrec[1] != "excludenothing") {
      obis_data <- obis_data %>%
        filter(!basisOfRecord %in% basisrec)
    }
    
    if (nrow(obis_data) > 0) {
      if (any(obis_data$absence)) {
        obis_data_absence <- obis_data %>%
          filter(absence)
      } else {
        obis_data_absence <- NULL
      }
    } else {
      obis_data_absence <- NULL
    }
    
  } else {
    obis_data <- NULL
    obis_data_absence <- NULL
  }
  
  if (length(gbif_ds) > 0 && !is.na(species_list$gbif_speciesKey[1])) {
    if (reader == "polars") {
      gbif_ds_pl <- pl$scan_parquet(
        file.path(gbif_ds, "**/*.parquet"), hive_partitioning = FALSE
      )
      
      gbif_data <- gbif_ds_pl$
        select(gbif_columns)$
        filter(pl$col("specieskey") == species_list$gbif_speciesKey)$
        collect()
      
      gbif_data <- gbif_data$to_data_frame()
      
      rm(gbif_ds_pl)
    } else if (reader == "arrow") {
      gbif_ds_pl <- open_dataset(gbif_ds)
      
      if (grepl("order=", gbif_ds_pl$files[1]) & 
          "gbif_order" %in% colnames(species_list)) {
        selorder <- ifelse(is.na(species_list$gbif_order[1]),
                           "noorder", species_list$gbif_order[1])
        
        gbif_data <- gbif_ds_pl %>%
          filter(order == selorder) %>%
          select(all_of(gbif_columns)) %>%
          filter(specieskey == species_list$gbif_speciesKey) %>%
          collect()
      } else {
        gbif_data <- gbif_ds_pl %>%
          select(all_of(gbif_columns)) %>%
          filter(specieskey == species_list$gbif_speciesKey) %>%
          collect()
      }
      
      rm(gbif_ds_pl)
    } else if (reader == "duckdb") {
      con <- dbConnect(duckdb())
      
      which_use <- ifelse(dir.exists(gbif_ds),
                          "folder_mode", "file_mode")
      
      if (which_use == "folder_mode") {
        
        if (grepl("order=", list.files(gbif_ds)[1])) {
          
          if ("gbif_order" %in% colnames(species_list)) {
            selorder <- ifelse(is.na(species_list$gbif_order[1]),
                               "noorder", species_list$gbif_order[1])
            selorder <- paste0("'", selorder, "'")
            selorder <- paste('"order" =', selorder, "and")
          } else {
            selorder <- ""
          }
          
          gbif_data <- dbGetQuery(con, glue::glue(
            '
             select {paste(gbif_columns, collapse = ", ")}
             from read_parquet("{gbif_ds}/*/*.parquet", hive_partitioning = true)
             where {selorder} specieskey = {species_list$gbif_speciesKey}
            '
          ))
        } else {
          gbif_data <- dbGetQuery(con, glue::glue(
            "
        select {paste(gbif_columns, collapse = ', ')}
        from read_parquet('{gbif_ds}/*')
        where specieskey = {species_list$gbif_speciesKey}
       "
          ))
        }
      } else {
        gbif_data <- dbGetQuery(con, glue::glue(
          "
        select {paste(gbif_columns, collapse = ', ')}
        from read_parquet('{gbif_ds}')
        where AphiaID = {species_list$gbif_speciesKey}
       "
        ))
      }
      
      DBI::dbDisconnect(con)
    }
    
    gbif_data <- gbif_data %>%
      mutate(taxonID = species_list$taxonID) %>%
      rename(decimalLatitude = decimallatitude,
             decimalLongitude = decimallongitude) %>%
      filter(year >= min_year) 
    
    if (basisrec[1] != "excludenothing") {
      gbif_data <- gbif_data %>%
        filter(!basisofrecord %in% basisrec)
    }
    
    if (nrow(gbif_data) > 0) {
      if (any(gbif_data$occurrencestatus == "ABSENT")) {
        gbif_data_absence <- gbif_data %>%
          filter(occurrencestatus == "ABSENT")
      } else {
        gbif_data_absence <- NULL
      }
    } else {
      gbif_data_absence <- NULL
    }
    
  } else {
    gbif_data <- NULL
    gbif_data_absence <- NULL
  }
  
  # Remove duplicates based on year and geohash
  if (!is.null(gbif_data) & !is.null(obis_data)) {
    non_dup_recs <- outqc_dup_check(list(obis = obis_data, gbif = gbif_data))
  } else if (!is.null(gbif_data)) {
    non_dup_recs <- outqc_dup_check(gbif_data)
  } else {
    non_dup_recs <- outqc_dup_check(obis_data)
  }
  
  rm(obis_data, gbif_data)
  
  # Remove from land
  if (remove_land) {
    if (!is.null(sdm_base)) {
      if (approximate_land) {
        pt_valid <- terra::extract(sdm_base, non_dup_recs[,lonlatdim])
        if (any(is.na(pt_valid[,2]))) {
          to_see <- which(is.na(pt_valid[,2]))
          to_see_crds <- non_dup_recs[to_see,lonlatdim]
          crds_to_cell <- terra::cellFromXY(sdm_base, as.data.frame(to_see_crds))
          
          near_cell <- lapply(crds_to_cell, function(cl){
            if (approximate_limit == 1) {
              direct <- "queen"
            } else {
              direct <- matrix(c(rep(1, 12), 0, rep(1,12)), 5, 5)
            }
            adj_cells <- terra::adjacent(sdm_base, cl, directions = direct)
            adj_values <- terra::extract(sdm_base, as.vector(adj_cells))
            if (all(is.na(adj_values[,1]))) {
              ret_df <- data.frame(cell = NA, x = NA, y = NA)
            } else {
              if (length(which(!is.na(adj_values[,1]))) == 1) {
                sel_cel <- which(!is.na(adj_values[,1]))
              } else {
                sel_cel <- sample(which(!is.na(adj_values[,1])), 1)
              }
              sel_cel <- as.vector(adj_cells)[sel_cel]
              sel_crds <- terra::xyFromCell(sdm_base, sel_cel)
              ret_df <- data.frame(cell = sel_cel, x = sel_crds[1,1], y = sel_crds[1,2])
            }
            return(ret_df)
          })
          
          near_cell <- do.call("rbind", near_cell)
          
          non_dup_recs[to_see,lonlatdim] <- near_cell[,2:3]
          
          pt_valid[to_see,2] <- near_cell[,1]
          
        }
      } else {
        pt_valid <- terra::extract(sdm_base, non_dup_recs[,lonlatdim])
      }
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
    non_dup_years <- non_dup_recs[,c(lonlatdim, "year")]
    
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
    if (!is.null(sdm_base)) {
      base <- terra::mask(sdm_base, sdm_base, inverse = T, updatevalue = 1)
    } else {
      sdm_base <- base <- outqc_get_base(0.05)
    }
    
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
    
    non_dup_recs$cell <- terra::cellFromXY(
      terra::rast(list.files("data/distances/", recursive = T,
                             pattern = ".tif", full.names = T)[1]),
      non_dup_recs[,lonlatdim]
    )
    
    non_dup_recs_sing <- non_dup_recs[!duplicated(non_dup_recs$cell),]
    
    if (geo_out) {
      
      if (nrow(non_dup_recs_sing) >= 5) {
        
        non_dup_geotag <- outqc_geo(non_dup_recs_sing[, lonlatdim],
                                    dist_folder = geo_dists_folder,
                                    limit_rem = out_threshold#0.005
                                    )
        
        colnames(non_dup_geotag) <- paste0("GeoTag_", colnames(non_dup_geotag))
        
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
                                    limit_rem = out_threshold)
        
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
      if (env_method == "ensemble") {
        sel_var <- colnames(non_dup_recs)[grepl("EnvTag_all", colnames(non_dup_recs))]
      } else {
        sel_var <- colnames(non_dup_recs)[grepl(paste0("EnvTag"), colnames(non_dup_recs))]  
        sel_var <- sel_var[grepl(env_method, sel_var)]
        if (env_method != "isoforest") {
          sel_var <- sel_var[grepl(env_variable, sel_var)]
        }
      }
      non_dup_recs <- non_dup_recs %>%
        {if (!narm_out) filter(., .data[[sel_var]] != 1 | is.na(.data[[sel_var]])) else filter(., .data[[sel_var]] != 1)}
      
    }
    if (any(grepl("GeoTag", colnames(non_dup_recs)))) {
      if (geo_method == "ensemble") {
        sel_var <- colnames(non_dup_recs)[grepl("GeoTag_all", colnames(non_dup_recs))]
      } else {
        sel_var <- colnames(non_dup_recs)[grepl(paste0("GeoTag"), colnames(non_dup_recs))]  
        sel_var <- sel_var[grepl(geo_method, sel_var)]
      }
      non_dup_recs <- non_dup_recs %>%
        {if (!narm_out) filter(., .data[[sel_var]] != 1 | is.na(.data[[sel_var]])) else filter(., .data[[sel_var]] != 1)}
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
        
        if (min(non_dup_recs_sums$impact) <= eval_max_impact) {
          hold_out_dataset <- non_dup_recs %>%
            filter(unID == sel_dataid)
          
          non_dup_recs <- non_dup_recs %>%
            filter(unID != sel_dataid)
        }
      }
      
      # Convert in 1 per cell
      all_pts <- per_cell(as.matrix(non_dup_recs[,lonlatdim]),
                          data_type = "fit_points")
      
      if (exists("hold_out_dataset")) {
        all_pts <- rbind(all_pts,
                         per_cell(as.matrix(hold_out_dataset[,lonlatdim]),
                                  data_type = "eval_points",
                                  data_sel = sel_dataid))
      }
      
      all_pts <- as.data.frame(all_pts)
      colnames(all_pts)[1:2] <- lonlatdim
      all_pts$taxonID <- species
      all_pts$species <- as.character(species_name)
      
      # Add min and max year information
      min_max_year <- function(target, sdm_base, rec_data) {
        cell_base <- sdm_base
        cell_base[] <- NA
        cell_base_id <- terra::cellFromXY(cell_base, as.data.frame(target[,lonlatdim]))
        target$cell <- cell_base_id
        target <- target %>%
          group_by(cell) %>%
          summarise(min_year = min(year),
                    max_year = max(year)) %>%
          filter(!is.na(cell))
        cell_base_min <- cell_base_max <- cell_base
        cell_base_min[target$cell] <- target$min_year
        cell_base_max[target$cell] <- target$max_year
        cell_base <- c(cell_base_min, cell_base_max)
        names(cell_base) <- c("min", "max")
        rec_data_ext <- terra::extract(cell_base, rec_data[,lonlatdim], ID = F)
        
        rec_data$min_year <- as.integer(rec_data_ext$min)
        rec_data$max_year <- as.integer(rec_data_ext$max)
        return(rec_data)
      }
      
      all_pts <- min_max_year(non_dup_years, sdm_base, all_pts)
      
      # Try to add tag reference to area of occurrence
      if (add_fao) {
        if (reader == "polars") {
          fao_areas <- pl$scan_parquet(fao_areas)
        } else {
          fao_areas <- arrow::open_dataset(fao_areas)
        }
        
        sp_sbase <- try(rfishbase::species(list(all_pts$species[1])))
        
        if (!inherits(sp_sbase, "try-error")) {
          if (nrow(sp_sbase) < 1) {
            sp_sbase <- try(rfishbase::species(list(all_pts$species[1]),
                                               server = "sealifebase"))
          }
        }
        
        if (!inherits(sp_sbase, "try-error") & nrow(sp_sbase) > 0) {
          if (reader == "polars") {
            fao_areas_sp <- fao_areas$filter(pl$col("SpecCode") == sp_sbase$SpecCode[1])$
              collect()
            fao_areas_sp <- fao_areas_sp$to_data_frame()
          } else {
            fao_areas_sp <- fao_areas %>%
              filter(SpecCode == sp_sbase$SpecCode[1]) %>%
              collect()
          }
          
          fao_areas_sp <- fao_areas_sp %>%
            filter(Status %in% c("endemic", "Endemic", "native", "Native", "introduced", "Introduced"))
          
          if (nrow(fao_areas_sp) > 0) {
            fao_shp <- terra::vect("data/shapefiles/World_Fao_Zones.shp")
            
            fao_shp_sel <- fao_shp[fao_shp$zone %in% fao_areas_sp$AreaCode]
            
            if (length(fao_shp_sel) > 0) {
              sf::sf_use_s2(FALSE)
              fao_shp_sel <- suppressMessages(
                suppressWarnings(
                  sf::st_buffer(sf::st_make_valid(sf::st_as_sf(fao_shp_sel)), 0.2)
                )
              )
              
              is_onarea <- terra::is.related(vect(all_pts[,1:2], geom = c("decimalLongitude", "decimalLatitude")),
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
      
      # Save absence data, if existent
      if (!is.null(gbif_data_absence) | !is.null(obis_data_absence)) {
        absence_data <- dplyr::bind_rows(obis_data_absence, gbif_data_absence)
        if (nrow(absence_data) > 0) {
          
          # Verify if it doesn't overlap with presence
          cell_base <- sdm_base
          cell_base[] <- NA
          
          cell_base[cellFromXY(cell_base, as.data.frame(absence_data[,lonlatdim]))] <- 0
          cell_base[cellFromXY(cell_base, as.data.frame(all_pts[,lonlatdim]))] <- 1
          
          absence_xy <- as.data.frame(cell_base, xy = T)
          absence_xy <- absence_xy[absence_xy[,3] == 0,]
          
          if (nrow(absence_xy) > 0) {
            colnames(absence_xy) <- c(lonlatdim, "data_type")
            absence_xy$data_type <- "absence_points"
            absence_xy$taxonID <- all_pts$taxonID[1]
            absence_xy$species <- all_pts$species[1]
            
            absence_xy <- min_max_year(absence_data, sdm_base, absence_xy)
            
            all_pts <- dplyr::bind_rows(all_pts, absence_xy)
          }
          
        }
      }
      
      arrow::write_parquet(all_pts,
                           paste0(species_outfolder, "/key=", species, ".parquet"))
    }
  }
  
  return(invisible(NULL))
}
