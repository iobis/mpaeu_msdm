### Outlier and Quality control

#' Get distances between all cells in a base raster, considering the continents
#' as barriers.
#'
#' @param base a base raster to use for calculating distances. Can be produced
#'   with [outqc_get_base()]
#' @param target optional, a dataframe containing the target x/y coordinates. If
#'   NULL, the target will be all cells.
#' @param agg_res an optional resolution to first aggregate the points. This can
#'   considerably speed up computation. Should be expressed in degrees. 
#'   This same value will be used to save a coarser resolution file. If \code{target}
#'   is not provided and \code{agg_res} is, then the base file is aggregated and
#'   all the cells for the coarser resolution are used (see details).
#' @param try_closest if \code{TRUE}, if the point lies in an invalid (empty) cell,
#' the function will try to find the closest valid point. It will look into the 8
#' adjacent cells. If there is still not one valid, this point is ignored.
#' Only relevant if target is supplied.
#' @param outfolder the folder where the files will be saved. If not existent,
#'   it will be created.
#' @param skip_exists if TRUE, it skip cells that were already calculated. Just
#'   set this to FALSE if you need to recalculate it.
#' @param mc_cores number of cores to use for parallel computation. Ignored on
#'   non-UNIX systems.
#' @param int_comp if \code{TRUE} (default), then values are first converted to
#'   integers and then the raster is saved with datatype "INT2U", what considerably
#'   reduce file size. See [terra::writeRaster()] for more details.
#' @param verbose enable messages
#'
#' @description This function will return the distances between each cell n and
#' all the other N(raster) cells. The distance is computed considering
#' continents as barriers.
#' 
#' ## Choosing the resolution
#' 
#' Choosing a very coarse resolution will be problematic, as this will cause some
#' barriers to be artificially created (e.g. small straits may become a barrier,
#' even if there is a passage). A recommended resolution is of 0.1 degree for the
#' base raster. At the same time, this will cause the computation to be much slower
#' while also increasing the size of saved cells. In this case, it's interesting to
#' aggregate the points to a coarser grid (e.g. 0.5 degrees). This will considerably
#' improve computational time while also reducing the size of saved files.
#' 
#' ## Intended use
#' 
#' This function will save files with the distance for cell i to all other cells on
#' the raster (according to the desired resolution). You can use those files directly,
#' but on the most recent version of the package the intended use is that after
#' getting all distances (by not supplying the \code{target}), you run the 
#' function `outqc_dist_tozarr`. This will create a `zarr` storage which can 
#' then be queried by the other functions.
#'
#' @return Files are saved on the disk.
#' @export
#'
#' @import parallel
#'
#' @examples
#' \dontrun{
#' # Get distance between cells in a 0.25 degrees resolution
#' outqc_get_distances(outqc_get_base(resolution = 0.25))
#' }
#' 
outqc_get_distances <- function(base,
                                target = NULL,
                                agg_res = NULL,
                                try_closest = TRUE,
                                outfolder = "distances",
                                skip_exists = TRUE,
                                mc_cores = NULL,
                                int_comp = TRUE,
                                verbose = TRUE) {
  
  # Verify if folder exists
  fs::dir_create(outfolder)
  
  if (is.null(target)) {
    if (!is.null(agg_res)) {
      if (verbose) cat("Aggregating by the resolution of", agg_res, "\n")
      new_base <- terra::aggregate(base, ceiling(agg_res/terra::res(base)[1]), na.rm = T)
      #Get all cells that are not NA
      tg_cells <- terra::as.data.frame(new_base, cells = T, xy = T)
    } else {
      #Get all cells that are not NA
      tg_cells <- terra::as.data.frame(base, cells = T, xy = T)
    }
  } else {
    
    check_valid <- function(id){
      if (is.na(is_valid[id])) {
        to_valid <- terra::adjacent(new_base, tg_cells[id,1], directions = 8)
        to_valid_vals <- terra::extract(base, terra::xyFromCell(new_base, to_valid[1,]))
        if (!all(is.na(to_valid_vals[,1]))) {
          new_valid <- to_valid[1,][which(!is.na(to_valid_vals[,1]))]
          new_valid_cell <- new_valid[1]
          new_valid_xy <- terra::xyFromCell(new_base, new_valid_cell)
          new_df <- data.frame(cell = new_valid_cell,
                               new_valid_xy)
        } else {
          new_df <- data.frame(cell = NA, x = NA, y = NA)
        }
      } else {
        new_df <- tg_cells[id,]
      }
      return(new_df)
    }
    
    if (!is.null(agg_res)) {
      if (verbose) cat("Aggregating target by the resolution of", agg_res, "\n")
      new_base <- terra::aggregate(base, ceiling(agg_res/terra::res(base)[1]))
      fgrid <- function(points, nrast){
        #nrast <- rast(resolution = agg_res)
        nrast[] <- NA
        nrast[terra::cellFromXY(nrast, data.frame(target))] <- 1
        terra::as.data.frame(nrast, cell = T, xy = T)[,1:3]
      }
      tg_cells <- fgrid(target, new_base)
      
      if (try_closest) {
        if (verbose) cat("Trying closest in case of NAs.\n")
        # Check if there are NAs
        is_valid <- unlist(terra::extract(base, tg_cells[,2:3], ID = F))
        
        new_valids <- lapply(1:nrow(tg_cells), check_valid)
        new_valids <- do.call("rbind", new_valids)
        if (verbose) cat(sum(is.na(new_valids[,1])), "NAs after check. If any, it will be disconsidered in computing distances. \n")
        tg_cells <- new_valids[!is.na(new_valids[,1]),]
      } 
      
      # Get the cells from the real base, that will be used for calculations
      tg_cells <- data.frame(cell = terra::cellFromXY(base, data.frame(tg_cells[,2:3])),
                             data.frame(tg_cells[,2:3]))
    } else {
      tg_cells <- data.frame(cell = terra::cellFromXY(base, data.frame(target)))
      
      tg_cells <- cbind(tg_cells, target)
      
      if (try_closest) {
        
        if (verbose) cat("Trying closest in case of NAs.\n")
        check_valid <- function(id){
          if (is.na(is_valid[id])) {
            to_valid <- terra::adjacent(base, tg_cells[id,1], directions = 8)
            to_valid_vals <- base[to_valid[1,]]
            if (!all(is.na(to_valid_vals[,1]))) {
              new_valid <- to_valid[1,][which(!is.na(to_valid_vals[,1]))]
              new_valid_cell <- new_valid[1]
              new_valid_xy <- terra::xyFromCell(base, new_valid_cell)
              new_df <- data.frame(cell = new_valid_cell,
                                   new_valid_xy)
            } else {
              new_df <- data.frame(cell = NA, x = NA, y = NA)
            }
          } else {
            new_df <- tg_cells[id,]
          }
          return(new_df)
        }
        
        is_valid <- unlist(terra::extract(base, tg_cells[,2:3], ID = F))
        
        new_valids <- lapply(1:nrow(tg_cells), check_valid)
        new_valids <- do.call("rbind", new_valids)
        if (verbose) cat(sum(is.na(new_valids[,1])), "NAs after check. If any, it will be disconsidered in computing distances. \n")
        tg_cells <- new_valids[!is.na(new_valids[,1]),]
        
      }
      
    }
  }
  
  .outqc_calc_distances <- function(index, # For which of the target cells get dist
                                    ag_fact = NULL,
                                    target_cells, # A data.frame with collumn 1 being the cells
                                    target_rast, # A base raster (should be the same used)
                                    outf, # Folder to save distance rasters
                                    comp_res) { # Compress results?
                                       
    
    # The distance will be computed to this cell
    tcell <- terra::cellFromXY(target_rast, target_cells[index,2:3])
    val <- target_rast[tcell][1,]
    if (is.na(val)) {
      p <- terra::vect(target_cells[index,2:3], geom = colnames(target_cells)[2:3], crs = "EPSG:4326")
      p_b <- terra::buffer(p, 20000)
      p_e <- terra::extract(target_rast, p_b, xy = T, ID = F)
      p_e <- p_e[!is.na(p_e[,1]),]
      if (nrow(p_e) > 0) {
        p_e <- terra::vect(p_e, geom = colnames(p_e)[2:3], crs = "EPSG:4326")
        near <- terra::nearest(p, p_e)
        near <- terra::as.data.frame(near, geom = "XY")
        tcell <- terra::cellFromXY(target_rast, near[1,c("x", "y")])
      } else {
        return(invisible(NULL))
      }
    }
    target_rast[tcell] <- 0
    #target_rast[target_cells[index,1]] <- 0
    
    # Compute distance and save in the outfolder
    dist_grid <- terra::gridDist(target_rast, scale = 1000)
    
    if (!is.null(ag_fact)) {
      dist_grid <- terra::aggregate(dist_grid, fact = ag_fact, fun = median, na.rm = T)
      
      new_cell <- terra::cellFromXY(dist_grid, target_cells[index,2:3])
      new_cell <- format(new_cell, trim = T, scientific = FALSE)
    } else {
      new_cell <- target_cells[index,1]
      new_cell <- format(new_cell, trim = T, scientific = FALSE)
    }
    
    if (comp_res) {
      dist_grid <- terra::as.int(dist_grid)
      
      terra::writeRaster(dist_grid,
                  filename = paste0(outf, "/cell_", new_cell, ".tif"),
                  datatype="INT2U",
                  overwrite = T)
    } else {
      terra::writeRaster(dist_grid,
                  filename = paste0(outf, "/cell_", new_cell, ".tif"),
                  overwrite = T)
    }
    
    return(invisible(NULL))
    
  }
  
  if (skip_exists) {
    exist_files <- list.files(outfolder)
    if (length(exist_files) > 0) {
      exist_files <- as.numeric(gsub(".tif", "", gsub("cell_", "", exist_files)))
    }
    if (!is.null(agg_res)) {
      new_cell <- terra::cellFromXY(new_base, tg_cells[,2:3])
    } else {
      new_cell <- tg_cells$cell
    }
    tg_cells <- data.frame(cell = tg_cells[!new_cell %in% exist_files, ])
  }
  
  if (is.null(mc_cores)) {
    if (.Platform$OS.type != "unix") {
      mc_cores <- 1
    } else {
      mc_cores <- getOption("mc.cores", 2L)
    }
  }
  
  if (nrow(tg_cells) > 0) {
    if (verbose) cat("Calculating distances...\n")
    mclapply(1:nrow(tg_cells), .outqc_calc_distances,
             target_cells = tg_cells,
             target_rast = base,
             outf = outfolder,
             ag_fact = ceiling(agg_res/terra::res(base)[1]),
             comp_res = int_comp,
             mc.cores = mc_cores)
  } else {
    if (verbose) cat("All distances already calculated and files in folder.\n")
  }
  
  return(invisible(NULL))
  
}



#' Get base raster for spatial outlier calculation
#'
#' @param resolution a numeric value indicating the resolution of the base raster (in degrees).
#' @param maskpol a shapefile (sf or terra::vect) that will be used as the mask for the barriers.
#'  Usually a world continents shapefile. If `NULL`, a world shapefile will be downloaded from
#'  NaturalEarth.
#'
#' @return a SpatRaster, with values 1 in all sea areas and NA on land.
#' @export
#'
#' @examples
#' \dontrun{
#' base <- outqc_get_base()
#' }
outqc_get_base <- function(resolution = 1,
                           maskpol = NULL){
  library(terra)
  
  if (is.null(maskpol)) {
    # Retrieve world polygon
    world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
    maskpol <- terra::vect(world)
  } else {
    if (class(maskpol)[1] != "SpatVector") {
      maskpol <- terra::vect(maskpol)
    }
  }
  
  # Gets a base raster
  base <- terra::rast(resolution = resolution)
  base[] <- 1
  names(base) <- "dist"
  terra::crs(base) <- terra::crs(maskpol)
  
  # Mask
  base <- terra::mask(base, maskpol, inverse = T)
  
  return(base)
}




#' Convert QC distances to Zarr format
#'
#' @param distfolder the folder where the distances, calculated using [outqc_get_distances()], were saved
#' @param remove_distances if \code{TRUE}, after converting to Zarr the calculated distances are deleted,
#' and only the Zarr folder is kept
#'
#' @return saved Zarr files
#' @export
#'
#' @examples
#' \dontrun{
#' outqc_dist_tozarr("distances")
#' }
outqc_dist_tozarr <- function(distfolder,
                              remove_distances = TRUE) {
  
  cli::cli_alert_info("Connecting to Python through {.code reticulate}")
  
  # Check if folder exist and have files
  if (!file.exists(distfolder) & length(list.files(distfolder)) < 1) {
    stop("Folder does not exist or does not have any file.")
  }
  
  # Save one that will be used as model
  files <- list.files(distfolder, full.names = T, pattern = "\\.tif")
  # Just to ensure that a file with valid points is used
  nr <- 0
  st <- 1
  while (nr < 1000) {
    r <- terra::rast(files[st])
    nr <- nrow(terra::as.data.frame(r))
    st <- st + 1
  }
  
  outfolder <- paste0(distfolder, "/", "distances.zarr")
  
  # Run in Python
  reticulate::py_run_string(glue::glue(
    "
import rasterio
import zarr
import xarray as xr
import os
import glob

# Function to stack GeoTIFF files and save as Zarr using Xarray
def geotiff_to_zarr(input_folder, output_path):
    # Get a list of GeoTIFF files in the input folder
    tiff_files = glob.glob(os.path.join(input_folder, '*.tif'))

    # Check if there are any files
    if not tiff_files:
        print('No GeoTIFF files found in the specified folder.')
        return

    # Read the first GeoTIFF file to get metadata
    with rasterio.open(tiff_files[0]) as src:
        meta = src.meta
        shape = src.shape

    # Create a list to store Xarray DataArray objects
    data_arrays = []
    file_basenames = []  # Store basenames for later use

    # Iterate through GeoTIFF files and create Xarray DataArray objects
    for tiff_file in tiff_files:
        file_basenames.append(int(os.path.splitext(os.path.basename(tiff_file))[0].replace('cell_', '')))
        with rasterio.open(tiff_file) as src:
            data = src.read(1)  # Assuming you want to read the first band
            data_arrays.append(xr.DataArray(data, dims=('y', 'x')))

    # Concatenate DataArray objects along a new dimension and add the basenames as a coordinate
    stacked_data = xr.concat(data_arrays, dim='time')
    stacked_data['time'] = file_basenames
    
    # Save the Xarray dataset as Zarr
    stacked_data.to_zarr(output_path, mode='w')


# Example usage:
input_folder = '<distfolder>'
output_path = '<outfolder>'

geotiff_to_zarr(input_folder, output_path)
    ", .open = "<", .close = ">"
  ))
  
  terra::writeRaster(r, paste0(outfolder, "/basedistances.tif"), overwrite = T)
  
  if (remove_distances & file.exists(outfolder)) {
    cli::cli_alert_info("Removing {.var .tif} files")
    rm(r)
    file.remove(files)
  }
  
  cli::cli_alert_success("Zarr file saved at {.path {outfolder}}")
  
  return(invisible(NULL))
  
}




#' Convert QC distances to RDS format
#'
#' @param distfolder the folder where the distances, calculated using [outqc_get_distances()], were saved
#' @param do_parallel if `TRUE` run in parallel using [furrr::future_map()]
#' @param parallel_cores number of cores to use for parallel processing. If `NULL` uses half of available cores
#' @param remove_distances if \code{TRUE}, after converting to RDS the calculated distances are deleted,
#' and only the RDS folder is kept
#'
#' @return saved RDS files
#' @export
#'
#' @examples
#' \dontrun{
#' outqc_dist_tords("distances")
#' }
outqc_dist_tords <- function(distfolder,
                             do_parallel = TRUE,
                             parallel_cores = NULL,
                             remove_distances = TRUE) {
  
  cli::cli_alert_info("Converting {.path {distfolder}} to RDS.")
  
  # Check if folder exist and have files
  if (!file.exists(distfolder) & length(list.files(distfolder)) < 1) {
    stop("Folder does not exist or does not have any file.")
  }
  
  # Save one that will be used as model
  files <- list.files(distfolder, full.names = T, pattern = "\\.tif")
  # Just to ensure that a file with valid points is used
  nr <- 0
  st <- 1
  while (nr < 1000) {
    r <- terra::rast(files[st])
    nr <- nrow(terra::as.data.frame(r))
    st <- st + 1
  }
  
  outfolder <- paste0(distfolder, "/", "distances.rds")
  fs::dir_create(outfolder)
  
  terra::writeRaster(r, paste0(outfolder, "/basedistances.tif"), overwrite = T)
  
  to_rds <- function(tf) {
    tf_r <- terra::rast(tf)
    tf_d <- terra::as.data.frame(tf_r, cell = T, na.rm = F)
    tf_d <- tf_d[order(tf_d$cell),]
    tf_v <- tf_d[,2]
    outn <- basename(tf)
    saveRDS(tf_v, file.path(outfolder, gsub(".tif", ".rds", outn)))
    return(invisible(NULL))
  }
  
  if (do_parallel) {
    require(future)
    require(furrr)
    
    if (is.null(parallel_cores)) {
      parallel_cores <- ceiling(availableCores()/2)
    }
    
    plan(multisession, workers = parallel_cores)
    si <- furrr::future_map(files, to_rds, .progress = T)
    plan(sequential)
  } else {
    for (i in cli::cli_progress_along(files)) {
      to_rds(files[i])
    }
  }
  
  available_cells <- basename(files)
  available_cells <- as.numeric(gsub("cell_", "", gsub("\\.tif", "", available_cells)))
  write.table(data.frame(cells = available_cells), 
              paste0(outfolder, "/availabledistances.txt"),
              row.names = F)
  
  if (remove_distances & file.exists(outfolder)) {
    cli::cli_alert_info("Removing {.var .tif} files")
    rm(r)
    file.remove(files)
  }
  
  cli::cli_alert_success("RDS files saved at {.path {outfolder}}")
  
  return(invisible(NULL))
  
}




#' Query the distance from points and k neighbours
#'
#' @param pts the points to be queried
#' @param kdist how much neighbours should be considered in the K nearest search
#' @param try_closest if you set \code{try_closest = TRUE} when getting the distances,
#' then you should also set it as \code{TRUE} here. Otherwise, the function will not be
#' able to find the valid cells that were computed.
#' @param distfolder the folder holding the distance layers produced with [outqc_get_distances()]
#' @param mode one of \code{rds}, \code{zarr}, \code{xarray} or \code{tif}.
#' Default is \code{zarr} which assumes that after running [outqc_get_distances()],
#' the function [outqc_dist_tords()] was used. It is recommended to use mode `rds` which
#' is less memory and time consuming
#' @param returnid if \code{TRUE} returns the ID
#' 
#' @description
#' Extract from a set of points the distance to the K nearest neighbours based on
#' distance layers generated with [outqc_get_distances()].
#' 
#' If \code{try_closest} was set to \code{TRUE} in when getting the distances,
#' then it's possible that some non-valid (empty) cells were omited. In that case
#' the function will return NA for those records.
#' 
#' Note: for \code{mode = "zarr"}, you need to have the `Rarr` package installed (through Bioconductor).
#' 
#' For mode \code{"xarray"}, you need Python and the [reticulate] package installed 
#' and the `xarray` Python package. This mode is memory intensive and is not (always)
#' returning the right results. It was kept only for development purposes.
#'
#' @return `matrix` with each row corresponding to a point, and each collumn the distance to the nearest K(k) neighbour.
#' @export
#'
#' @import parallel reticulate
#'
#' @examples
#' \dontrun{
#' outqc_query_distances(sp_pts)
#' }
outqc_query_distances <- function(pts,
                                  kdist = NULL,
                                  try_closest = TRUE,
                                  distfolder = "distances",
                                  mode = "rds",
                                  returnid = TRUE) {
  
  if (is.null(kdist)) {
    kdist <- nrow(pts)-1
  } else {
    if (kdist >= nrow(pts)) {
      warning("kdist should be lower then the number of points. Reduced to nrow(pts)-1.")
      kdist <- nrow(pts)-1
    } else {
      if (kdist == 0) {
        warning("kdist should be at least 1. Using 1 kdist")
        kdist <- 1
      }
    }
  }
  
  if (mode == "tif") {
    base <- terra::rast(list.files(distfolder, pattern = ".tif", full.names = T)[1])
  } else if (mode == "rds") {
    base <- terra::rast(paste0(distfolder, "/distances.rds/basedistances.tif"))
  } else {
    base <- terra::rast(paste0(distfolder, "/distances.zarr/basedistances.tif"))
  }
  
  pts <- as.data.frame(pts)
  
  cells_index <- terra::cellFromXY(base, pts)
  
  pts$cell <- cells_index
  
  cells_index <- unique(cells_index)
  
  
  extract_dist <- function(index, cells_index) {
    
    target <- cells_index[index]
    
    other_pts <- cells_index[-index]
    
    if (file.exists(paste0(distfolder, "/cell_", format(target, scientific = FALSE), ".tif"))) {
      sel_dist_rast <- suppressWarnings(
        try(terra::rast(paste0(distfolder, "/cell_", format(target, scientific = FALSE), ".tif")),
            silent = T)
      )
      
      if (class(sel_dist_rast)[1] == "try-error") {
        warning("File ", paste0(distfolder, "/cell_", format(target, scientific = FALSE), ".tif"), " is corrupted. Removing it and ignoring.\n")
        file.remove(paste0(distfolder, "/cell_", format(target, scientific = FALSE), ".tif"))
        return(list(dist = NA, id = NA))
      }
      
      other_dist <- sel_dist_rast[other_pts][,1]
      
      k_near <- t(as.matrix(order(other_dist)[1:kdist]))
      
      k_dist <- t(as.matrix(other_dist[k_near[1,]]))
      
      return(list(dist = k_dist,
                  id = k_near))
    } else {
      if (try_closest) {
        to_valid <- terra::adjacent(base, target, directions = 8)
        to_valid_files <- paste0(distfolder, "/cell_", format(to_valid[1,], trim = T, scientific = FALSE), ".tif")
        to_valid_files_ex <- file.exists(to_valid_files)
        if (any(to_valid_files_ex)) {
          new_file <- to_valid_files[which(to_valid_files_ex == TRUE)[1]]
          
          sel_dist_rast <- suppressWarnings(try(terra::rast(new_file), silent = T))
          
          if (class(sel_dist_rast)[1] == "try-error") {
            warning("File ", new_file, " is corrupted. Removing it and ignoring.\n")
            file.remove(new_file)
            return(list(dist = NA, id = NA))
          }
          
          other_dist <- sel_dist_rast[other_pts][,1]
          
          k_near <- t(as.matrix(order(other_dist)[1:kdist]))
          
          k_dist <- t(as.matrix(other_dist[k_near[1,]]))
          
          return(list(dist = k_dist,
                      id = k_near))
        } else {
          return(list(dist = NA, id = NA))
        }
      } else {
        return(list(dist = NA, id = NA))
      }
    }
  }
  
  
  if (mode == "tif") {
    pts_kdist <- lapply(1:length(cells_index), extract_dist, cells_index = cells_index)
  } else if (mode == "xarray") {
    
    # Import xarray
    xr <- import("xarray")
    
    # Open distances
    dists <- xr$open_zarr(paste0(distfolder, "/distances.zarr"))
    
    # Rename variable for easier handling
    dists <- dists$rename_vars('__xarray_dataarray_variable__' = 'distances')
    
    # Retrieve existing cells
    exist_cells <- dists$time$values
    # For compatibility with new Xarray versions/reticulate problems:
    if (!inherits(exist_cells, c("array", "integer", "double", "numeric"))) {
      exist_cells <- dists$time$values$tolist()
    }
    
    
    # Check if all are valid
    cell_list <- cells_index
    
    cell_list_inv <- cell_list[!cell_list %in% exist_cells]
    
    # Check for closest or remove
    if (length(cell_list_inv) > 0 & try_closest) {
      
      to_valid <- terra::adjacent(base, cell_list_inv, directions = 8)
      
      to_valid <- unname(apply(to_valid, 1, function(x) {
        x_val <- x %in% exist_cells
        if (any(x_val)) {
          x <- x[x_val]
          x[1]
        } else {
          NA
        }
      }))
      
      cell_list[!cell_list %in% exist_cells] <- to_valid
      
    } else {
      cell_list[!cell_list %in% exist_cells] <- NA
    }
    
    cell_equiv <- data.frame(cell_orig = cells_index, cell_new = cell_list)
    
    cell_to_do <- na.omit(unique(cell_list))
    
    cell_cols <- terra::colFromCell(base, cell_to_do)
    cell_rows <- terra::rowFromCell(base, cell_to_do)
    
    cell_values <- list()
    
    vals <- dists$sel(time = cell_to_do[order(cell_to_do)])
    
    for (z in 1:length(cell_to_do)) {
      vals_sel <- vals$distances$sel(x = as.integer(cell_cols[z]-1),
                                     y = as.integer(cell_rows[z]-1))
      cell_values[[z]] <- vals_sel$values
      
      # For compatibility with new Xarray versions/reticulate problems:
      if (!inherits(cell_values[[z]], c("array", "integer", "double", "numeric"))) {
        cell_values[[z]] <- vals_sel$values$tolist()
      }
      
      cell_values[[z]] <- cell_values[[z]][order(cell_values[[z]])][2:(kdist+1)]
      
    }
    
    # Join information with the original values
    all_values <- do.call("rbind", cell_values)
    
    all_values <- cbind(cell_new = cell_to_do, all_values)
    
    all_values <- dplyr::left_join(cell_equiv, as.data.frame(all_values), by = "cell_new")
    names(all_values)[1] <- "cell"
    
    # Join with the initial points list
    final_values <- dplyr::left_join(pts, all_values, by = "cell")
    
  } else if (mode == "zarr") {
    
    require(Rarr)
    
    zarr_address <- paste0(distfolder, "/distances.zarr")
    zarr_details <- zarr_overview(zarr_address, as_data_frame = T)
    
    z_time <- zarr_details[grepl("time", zarr_details$path),]
    z_data <- zarr_details[!grepl("time", zarr_details$path),]

    # Retrieve existing cells
    exist_cells <- read_zarr_array(z_time$path)
    
    # Check if all are valid
    cell_list <- cells_index
    
    cell_list_inv <- cell_list[!cell_list %in% exist_cells]
    
    # Check for closest or remove
    if (length(cell_list_inv) > 0 & try_closest) {
      
      to_valid <- terra::adjacent(base, cell_list_inv, directions = 8)
      
      to_valid <- unname(apply(to_valid, 1, function(x) {
        x_val <- x %in% exist_cells
        if (any(x_val)) {
          x <- x[x_val]
          x[1]
        } else {
          NA
        }
      }))
      
      cell_list[!cell_list %in% exist_cells] <- to_valid
      
    } else {
      cell_list[!cell_list %in% exist_cells] <- NA
    }
    
    cell_equiv <- data.frame(cell_orig = cells_index, cell_new = cell_list)
    
    cell_to_do <- na.omit(unique(cell_list))
    
    cell_cols <- terra::colFromCell(base, cell_to_do)
    cell_rows <- terra::rowFromCell(base, cell_to_do)
    
    # get_time <- function(cell) {
    #   which(exist_cells == cell)
    # }
    # 
    # time_as_cell <- unlist(lapply(cell_to_do, get_time))
    
    time_as_cell <- match(cell_to_do, exist_cells)
    
    extract_fz_info <- function(rows, cols) {
      min_r <- min(rows)
      max_r <- max(rows)
      r_seq <- min_r:max_r
      
      min_c <- min(cols)
      max_c <- max(cols)
      c_seq <- min_c:max_c
      
      rows_eq <- match(rows, r_seq)
      cols_eq <- match(cols, c_seq)
      
      list(r_seq = r_seq, c_seq = c_seq,
           req = rows_eq, ceq = cols_eq)
    }
    
    fzinfo <- extract_fz_info(cell_rows,
                              cell_cols)
    
    extract_from_zarr <- function(fz_info, t_index, kdist, data_zarr_path) {
      
      t_index_sp <- split(seq_along(t_index), ceiling(seq_along(t_index)/1000))
      
      extract_values <- function(slice, rows, cols, kdist) {
        sv <- slice[cbind(rows, cols)]
        return(sv[order(sv)][2:(kdist+1)])
      }
      
      zar_vf <- lapply(t_index_sp, function(ti){
        index <- list(t_index[ti], fz_info$r_seq, fz_info$c_seq)
        zar <- read_zarr_array(data_zarr_path, index = index)
        
        zar_v <- apply(zar, 1, extract_values,
                       rows = fz_info$req, cols = fz_info$ceq, kdist = kdist)
        return(zar_v)
      })
      
      zar_vf <- do.call("cbind", zar_vf)
      
      return(zar_vf)
    }
    
    # Join information with the original values
    all_values <- extract_from_zarr(fzinfo,
                                    time_as_cell,
                                    kdist,
                                    z_data$path)
    all_values <- t(all_values)
    
    all_values <- cbind(cell_new = cell_to_do, all_values)
    
    all_values <- dplyr::left_join(cell_equiv, as.data.frame(all_values), by = "cell_new")
    names(all_values)[1] <- "cell"
    
    # Join with the initial points list
    final_values <- dplyr::left_join(pts, all_values, by = "cell")
    
  } else if (mode == "rds") {
    
    # Retrieve existing cells
    exist_cells <- read.table(file.path(distfolder,
                                        "distances.rds",
                                        "availabledistances.txt"), header = T)
    exist_cells <- exist_cells[,1]
    
    # Check if all are valid
    cell_list <- cells_index
    
    cell_list_inv <- cell_list[!cell_list %in% exist_cells]
    
    # Check for closest or remove
    if (length(cell_list_inv) > 0 & try_closest) {
      
      to_valid <- terra::adjacent(base, cell_list_inv, directions = 8)
      
      to_valid <- unname(apply(to_valid, 1, function(x) {
        x_val <- x %in% exist_cells
        if (any(x_val)) {
          x <- x[x_val]
          x[1]
        } else {
          NA
        }
      }))
      
      cell_list[!cell_list %in% exist_cells] <- to_valid
      
    } else {
      cell_list[!cell_list %in% exist_cells] <- NA
    }
    
    cell_equiv <- data.frame(cell_orig = cells_index, cell_new = cell_list)
    
    cell_to_do <- na.omit(unique(cell_list))
    
    read_and_extract <- function(target, cells, kdist) {
      values <- readRDS(file.path(distfolder, "distances.rds", paste0("cell_", target, ".rds")))
      sv <- values[cells]
      sv[order(sv)][2:(kdist+1)]
    }
    
    all_values <- matrix(nrow = length(cell_to_do), ncol = kdist)
    
    for (kl in 1:length(cell_to_do)) {
      all_values[kl,] <- read_and_extract(cell_to_do[kl], cell_to_do, kdist)
    }
    
    all_values <- as.data.frame(all_values)
    all_values <- cbind(cell_new = cell_to_do, all_values)
    
    all_values <- dplyr::left_join(cell_equiv, as.data.frame(all_values), by = "cell_new")
    names(all_values)[1] <- "cell"
    
    # Join with the initial points list
    final_values <- dplyr::left_join(pts, all_values, by = "cell")
    
  }
  
  if (returnid & mode == "tif") {
    pts_id <- lapply(pts_kdist, function(x) x$id)
    pts_id <- data.frame(do.call("rbind", pts_id))
    
    pts_id$cell <- cells_index
    
    pts_join <- dplyr::left_join(pts, pts_id, by = "cell")
    
    pts_dist <- pts_join[,!colnames(pts_join) %in% c("cell", "decimalLongitude", "decimalLatitude", "x", "y")]
    
    colnames(pts_id) <- paste0("K", 1:kdist)
  }
  
  if (mode == "tif") {
    pts_dist <- lapply(pts_kdist, function(x) x$dist)
    pts_dist <- data.frame(do.call("rbind", pts_dist))
    
    pts_dist$cell <- cells_index
    
    pts_join <- dplyr::left_join(pts, pts_dist, by = "cell")
    
    pts_dist <- pts_join[,!colnames(pts_join) %in% c("cell", "decimalLongitude", "decimalLatitude", "x", "y")]
    
    colnames(pts_dist) <- paste0("K", 1:kdist)
  } else {
    pts_dist <- final_values[,!colnames(final_values) %in% c("cell", "cell_new", "decimalLongitude", "decimalLatitude", "x", "y")]
    
    colnames(pts_dist) <- paste0("K", 1:kdist)
  }
  
  if (returnid & mode == "tif") {
    return(list(
      dist = pts_dist,
      id = pts_id
    ))
  } else {
    return(pts_dist)
  }
}
