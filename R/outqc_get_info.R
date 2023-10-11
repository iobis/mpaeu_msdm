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
#' improve computational time while also reducing the size of save files.
#' 
#' ## Intended use
#' 
#' This function will save files with the distance for cell i to all other cells on
#' the raster (according to the desired resolution). And why to save? Because later you
#' can query the distances. The flux is for a set of points z, get all distances from
#' z{i} to z{n-i}. Then query the distances to perform outlier calculations. When you run
#' this function on the next dataset, it will run faster because potentially some of the cells 
#' were already accounted for in the previous calculation. Thus, there is a time improvement each
#' time the function is run, up to the point that all cells are calculated.
#' 
#' Alternatively you can pre-calculate the distance to all cells by not supplying the
#' \code{target}. 
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
      new_base <- terra::aggregate(base, agg_res/terra::res(base)[1])
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
      new_base <- terra::aggregate(base, agg_res/terra::res(base)[1])
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
    target_rast[target_cells[index,1]] <- 0
    
    # Compute distance and save in the outfolder
    dist_grid <- terra::gridDist(target_rast, scale = 1000)
    
    if (!is.null(ag_fact)) {
      dist_grid <- terra::aggregate(dist_grid, fact = ag_fact, fun = min, na.rm = T)

      new_cell <- terra::cellFromXY(dist_grid, target_cells[index,2:3])
    } else {
      new_cell <- target_cells[index,1]
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
             ag_fact = (agg_res/terra::res(base)[1]),
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
  crs(base) <- crs(maskpol)
  
  # Mask
  base <- terra::mask(base, maskpol, inverse = T)
  
  return(base)
}




#' Query the distance from points and k neighbours
#'
#' @param pts the points to be queried
#' @param base a base rast (of the same resolution as the one used when 
#' calculating the distances), so the points can be indexed as cells. If NULL, 
#' one of the distance layers from \code{distfolder} will be used instead.
#' @param kdist how much neighbours should be considered in the K nearest search
#' @param try_closest if you set \code{try_closest = TRUE} when getting the distances,
#' then you should also set it as \code{TRUE} here. Otherwise, the function will not be
#' able to find the valid cells that were computed.
#' @param distfolder the folder holding the distance layers produced with [outqc_get_distances()]
#' @param parallel run query in parallel (recommended for larger datasets). 
#' For now, will not work on Windows based systems.
#' @param mc_cores number of cores to use for parallel computation. 
#' 
#' @description
#' Extract from a set of points the distance to the K nearest neighbours based on
#' distance layers generated with [outqc_get_distances()].
#' 
#' If \code{try_closest} was set to \code{TRUE} in when getting the distances,
#' then it's possible that some non-valid (empty) cells were omited. In that case
#' the function will return NA for those records.
#'
#' @return `matrix` with each row corresponding to a point, and each collumn the distance to the nearest K(k) neighbour.
#' @export
#'
#' @import parallel
#'
#' @examples
#' \dontrun{
#' outqc_query_distances(sp_pts)
#' }
outqc_query_distances <- function(pts,
                                  base = NULL,
                                  kdist = NULL,
                                  try_closest = TRUE,
                                  distfolder = "distances",
                                  returnid = TRUE,
                                  parallel = TRUE,
                                  mc_cores = NULL) {
  
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
  
  if (is.null(base)) {
    base <- terra::rast(list.files(distfolder, pattern = ".tif", full.names = T)[1])
  }
  
  #pts.dist <- matrix(nrow = nrow(pts), ncol = kdist)
  
  pts <- as.data.frame(pts)
  
  cells_index <- terra::cellFromXY(base, pts)
  
  pts$cell <- cells_index
  
  cells_index <- unique(cells_index)
  
  
  extract_dist <- function(index, cells_index) {
    
    target <- cells_index[index]
    
    other_pts <- cells_index[-index]
    
    if (file.exists(paste0(distfolder, "/cell_", target, ".tif"))) {
      sel_dist_rast <- terra::rast(paste0(distfolder, "/cell_", target, ".tif"))
      
      other_dist <- sel_dist_rast[other_pts][,1]
      
      k_near <- t(as.matrix(order(other_dist)[1:kdist]))
      
      k_dist <- t(as.matrix(other_dist[k_near[1,]]))
      
      return(list(dist = k_dist,
                  id = k_near))
    } else {
      if (try_closest) {
        to_valid <- terra::adjacent(base, target, directions = 8)
        to_valid_files <- paste0(distfolder, "/cell_", to_valid[1,], ".tif")
        to_valid_files_ex <- file.exists(to_valid_files)
        if (any(to_valid_files_ex)) {
          new_file <- to_valid_files[which(to_valid_files_ex == TRUE)[1]]
          
          sel_dist_rast <- terra::rast(new_file)
          
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
  
  if (is.null(mc_cores)) {
    if (.Platform$OS.type != "unix") {
      mc_cores <- 1
    } else {
      mc_cores <- getOption("mc.cores", 2L)
    }
  }
  
  if (parallel) {
    pts_kdist <- mclapply(1:length(cells_index), extract_dist, cells_index = cells_index)
  } else {
    pts_kdist <- lapply(1:length(cells_index), extract_dist, cells_index = cells_index)
  }
  
  if (returnid) {
    pts_id <- lapply(pts_kdist, function(x) x$id)
    pts_id <- data.frame(do.call("rbind", pts_id))
    
    pts_id$cell <- cells_index
    
    pts_join <- dplyr::left_join(pts, pts_id, by = "cell")
    
    pts_dist <- pts_join[,!colnames(pts_join) %in% c("cell", "decimalLongitude", "decimalLatitude", "x", "y")]
    
    colnames(pts_id) <- paste0("K", 1:kdist)
  }
  
  pts_dist <- lapply(pts_kdist, function(x) x$dist)
  pts_dist <- data.frame(do.call("rbind", pts_dist))
  
  pts_dist$cell <- cells_index
  
  pts_join <- dplyr::left_join(pts, pts_dist, by = "cell")
  
  pts_dist <- pts_join[,!colnames(pts_join) %in% c("cell", "decimalLongitude", "decimalLatitude", "x", "y")]
  
  colnames(pts_dist) <- paste0("K", 1:kdist)
  
  if (returnid) {
    return(list(
      dist = pts_dist,
      id = pts_id
    ))
  } else {
    return(pts_dist)
  }
}
