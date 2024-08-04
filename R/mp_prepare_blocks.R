############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################

#' Get cross-validation blocks
#'
#' @param sdm_data an object generated with mp_prepare_data containing the data
#' @param method either 'blockcv' to let the [blockCV] package generate the grid
#'   or 'manual' to use a manual grid
#' @param block_types a vector of one or more of "spatial_grid", "spatial_lat"
#'   and "random". See details for more information
#' @param nfolds number of folds
#' @param auto_range if \code{TRUE}, the function
#'   [blockCV::cv_spatial_autocor()] is used to search an ideal size for the
#'   blocks. Ignored if method is 'manual'
#' @param range a range for the blocks. Ignored if method is 'manual' or
#'   'auto_range' is set to \code{TRUE}
#' @param lat_blocks the number of latitudinal blocks, if block_types
#'   'spatial_lat' is chosen
#' @param manual_shp the manual grid in sf format. Should be a named list
#'   containing the grids for one or both of spatial_grid and spatial_lat
#' @param n_iterate number of iterations for random assigning the blocks into folds
#' @param retry_if_zero if this is set to \code{TRUE}, in case any of the folds
#'   have a 0 parameter (either for training or testing presence/background),
#'   the function will retry to assign folds, by randomly decreasing or
#'   increasing the resolution of the original grid. The number of new trials is
#'   5 by default, but you can change this value by supplying an integer instead of
#'   \code{TRUE} for this parameter. You can also retry in case of a different value
#'   (e.g. only 1 test point) by changing the \code{min_class} parameter
#' @param min_block_size if \code{retry_if_zero = TRUE}, then define the minimum size
#' the block can have
#' @param max_block_size if \code{retry_if_zero = TRUE}, then define the maximum size
#' the block can have
#' @param min_class if \code{retry_if_zero = TRUE}, the function will test if any of the folds
#' have a zero in any of the classes, but instead of 0 it's possible to set a different value,
#' for example 1 to have folds with at least 2 points in each class. Note that if you change this value for
#' something too large and your points are not well distributed, then the function will probably fail
#' in assigning the groups
#' @param verbose if \code{TRUE}, print messages
#'
#' @details
#' Three types of cross-validation blocks are available:
#' - spatial_grid: a spatial grid of blocks of size \code{range} is generated
#' - spatial_lat: a spatial grid of latitudinal bins (= the number of \code{lat_blocks})
#' - random: points are randomly assigned in the folds (i.e., not a spatial cross-validation)
#' 
#'
#' @return an object of class sdm_dat containing the blocks
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
mp_prepare_blocks <- function(sdm_data, method = "blockcv", 
                              block_types = c("spatial_grid", "spatial_lat", "random"),
                              nfolds = 5, auto_range = FALSE, range = NULL, 
                              lat_blocks = 50, manual_shp = NULL,
                              n_iterate = 50, retry_if_zero = FALSE,
                              min_block_size = 0.2, max_block_size = 20,
                              min_class = 0,
                              verbose = TRUE) {
  
  if (class(sdm_data)[1] != "sdm_dat") {
    cli::cli_abort("Supply an object of class sdm_dat produced with mp_prepare_data.")
  }
  
  if (!method %in% c("manual", "blockcv")) {
    cli::cli_abort("Method should be blockcv or manual.")
  }
  
  if (!any(block_types %in% c("spatial_grid", "spatial_lat", "random"))) {
    cli::cli_abort("block_types should be one (or more) of spatial_grid, spatial_lat, random.")
  }
  
  if (!isFALSE(retry_if_zero)) {
    if (is.numeric(retry_if_zero)) {
      nretry <- retry_if_zero
      retry_if_zero <- TRUE
    } else {
      nretry <- 5
    }
  }
  
  gr <- NULL
  folds <- list()
  
  sf_pts <- sf::st_as_sf(cbind(sdm_data$training,
                               sdm_data$coord_training), 
                         coords = c("decimalLongitude", "decimalLatitude"),
                         crs = 4326)
  
  if (method == "blockcv") {
    
    if (verbose) cli::cli_inform(c("v" = "Using method blockCV"))
    
    if (auto_range) {
      sf::sf_use_s2(FALSE)
      auto_cor <- blockCV::cv_spatial_autocor(x = sf::st_make_valid(sf_pts), column = "presence")
      range = auto_cor$range
    } else {
      sf::sf_use_s2(FALSE)
      if (is.null(range)) {
        if (verbose) cli::cli_alert_danger("range was not supplied. Using 20km as range.")
        range <- 20000
      }
    }
    
    # Improve to ensure all species will be equally distributed in blocks!
    if ("spatial_grid" %in% block_types) {
      if (verbose) cli::cli_alert_info("Generating grid blocks.")
      blocks <- blockCV::cv_spatial(sf_pts, column = "presence", iteration = n_iterate, size = range,
                      hexagon = FALSE, k = nfolds, progress = F, report = T, plot = F, selection = "systematic")
      folds <- c(folds, spatial_grid = list(blocks$folds_ids))
    }
    
    if ("spatial_lat" %in% block_types) {
      if (verbose) cli::cli_alert_info("Generating latitudinal blocks.")
      blocks <- blockCV::cv_spatial(sf_pts, column = "presence", iteration = n_iterate, rows_cols = c(lat_blocks, 1),
                           hexagon = FALSE, k = nfolds, progress = F, report = T, plot = F, extend = 0.5)
      folds <- c(folds, spatial_lat = list(blocks$folds_ids))
    }
    
    if ("random" %in% block_types) {
      if (verbose) cli::cli_alert_info("Generating random samples.")
      blocks <- sample(1:5, nrow(sf_pts), replace = T)
      folds <- c(folds, random = list(blocks))
    }
    
  } else {
    
    if (verbose) cli::cli_inform(c("v" = "Using method manual"))
    
    if ("spatial_grid" %in% block_types) {
      if (verbose) cli::cli_alert_info("Generating grid blocks.")
      manual_grid <- manual_shp$spatial_grid

      blocks <- manual_block(sf_pts,
                             user_grid = manual_grid,
                             k = nfolds,
                             iterations = n_iterate)
      
      if (retry_if_zero) {
        if (any(apply(blocks$records, 1, function(x) ifelse(x <= min_class, TRUE, FALSE)))) {
          if (verbose) cli::cli_alert_warning("Empty (or different than {.var min_class}) classes in one or more folds. Retrying with different resolution...")
          init <- 1
          bottom_reached <- top_reached <- FALSE

          while(
            any(apply(blocks$records, 1, function(x) ifelse(x <= min_class, TRUE, FALSE))) &
            init <= nretry
          ) {

            if (top_reached) {
              prev_res <- max_block_size
            } else if (bottom_reached) {
              prev_res <- min_block_size
            } else {
              prev_res <- terra::res(manual_grid)[1]
            }

            if (prev_res >= max_block_size | ncell(manual_grid) < nfolds) {
              cfact <- -2
              top_reached <- TRUE
            } else if (prev_res <= min_block_size) {
              cfact <- 2
              bottom_reached <- TRUE
            } else {
              cfact <- sample(c(-2, 2), 1)
            }

            if (cfact > 0) {
              manual_grid <- terra::aggregate(manual_grid, cfact)
            } else {
              manual_grid <- terra::disagg(manual_grid, (cfact * -1))
            }

            blocks <- manual_block(sf_pts,
                                   user_grid = manual_grid,
                                   k = nfolds,
                                   iterations = n_iterate)
            init <- init + 1

            if (top_reached & bottom_reached) {
              init <- nretry + 1
            }
          }
        }
        if (any(apply(blocks$records, 1, function(x) ifelse(x <= min_class, TRUE, FALSE)))) {
          if (verbose) cli::cli_alert_danger("It was not possible to avoid empty classes in one or more folds.")
        }
      }

      folds <- c(folds, spatial_grid = list(blocks$folds))
      gr <- c(gr, spatial_grid = terra::res(manual_grid)[1])
    }

    if ("spatial_lat" %in% block_types) {
      if (verbose) cli::cli_alert_info("Generating latitudinal blocks.")
      manual_grid <- manual_shp$spatial_lat

      blocks <- manual_block(sf_pts,
                             user_grid = manual_grid,
                             k = nfolds,
                             iterations = n_iterate)
      if (retry_if_zero) {
        if (any(apply(blocks$records, 1, function(x) ifelse(x <= min_class, TRUE, FALSE)))) {
          if (verbose) cli::cli_alert_warning("Empty (or different than {.var min_class}) classes in one or more folds. Retrying with different resolution...")
          init <- 1
          bottom_reached <- top_reached <- FALSE
          
          while(
            any(apply(blocks$records, 1, function(x) ifelse(x <= min_class, TRUE, FALSE))) &
            init <= nretry
          ) {
            
            if (top_reached) {
              prev_res <- max_block_size
            } else if (bottom_reached) {
              prev_res <- min_block_size
            } else {
              prev_res <- terra::res(manual_grid)[1]
            }
            
            if (prev_res >= max_block_size | ncell(manual_grid) < nfolds) {
              cfact <- -2
              top_reached <- TRUE
            } else if (prev_res <= min_block_size) {
              cfact <- 2
              bottom_reached <- TRUE
            } else {
              cfact <- sample(c(-2, 2), 1)
            }
            
            if (cfact > 0) {
              manual_grid <- terra::aggregate(manual_grid, cfact)
            } else {
              manual_grid <- terra::disagg(manual_grid, (cfact * -1))
            }
            
            blocks <- manual_block(sf_pts,
                                   user_grid = manual_grid,
                                   k = nfolds,
                                   iterations = n_iterate)
            init <- init + 1
            
            if (top_reached & bottom_reached) {
              init <- nretry + 1
            }
          }
        }
        if (any(apply(blocks$records, 1, function(x) ifelse(x <= min_class, TRUE, FALSE)))) {
          if (verbose) cli::cli_alert_danger("It was not possible to avoid empty classes in one or more folds.")
        }
      }

      folds <- c(folds, spatial_lat = list(blocks$folds))
      gr <- c(gr, spatial_lat = terra::nrow(manual_grid))
    }

    if ("random" %in% block_types) {
      if (verbose) cli::cli_alert_info("Generating random samples.")
      blocks <- sample(1:5, nrow(sf_pts), replace = T)
      folds <- c(folds, random = list(blocks))
    }
    
  }
  
  sdm_data$blocks <- list(n = nfolds,
                          folds = folds,
                          gen_method = method,
                          grid_resolution = gr)
  
  return(sdm_data)
}



#' Fast spatial cross-validation blocks generation
#' 
#' This function will generate cross-validation blocks randomly based on a grid,
#' which may also be of latitudinal blocks. This function is simpler than the
#' [blockCV::cv_spatial()], but is aimed to be much faster.
#'
#' @param data_pts points in SF format. Should have a collumn named "presence"
#' containing the type of point (0 for background/absence and 1 for presence)
#' @param user_grid an user supplied SpatRaster ([terra]). If not supplied,
#' the function will generate a grid of resolution 1 degree based on the extent
#' of the points.
#' @param k number of folds
#' @param iterations number of iterations for sampling the points.
#'
#' @return a list containing the folds, number of folds and division of records per fold.
#' @export
#'
#' @examples
#' \dontrun{
#' user_grid <- rast(ext(-180, 180, -90, 90), resolution = 1)
#' blocks <- manual_block(species_points, user_grid, k = 5)
#' } 
#' 
manual_block <- function(data_pts, user_grid = NULL, k, iterations = 50) {
  
  # Check if grid was supplied
  if (is.null(user_grid)) {
    cli::cli_alert_warning("Grid not supplied, generating grid with resolution 1 based on points extent.")
    user_grid <- terra::rast(terra::ext(data_pts), resolution = 1)
  }
  
  # create objects to store results
  result <- list()
  folds <- list()
  
  # Get cell number
  which_cell <- terra::cellFromXY(user_grid, sf::st_coordinates(data_pts))
  
  # Sample for i iterations
  for (i in 1:iterations) {
    user_grid[] <- sample(1:k, ncell(user_grid), replace = T)
    #cell_vals <- user_grid[which_cell][,1]
    cell_vals <- terra::extract(user_grid, which_cell)[,1]
    df_res <- data.frame(fold = 1:k, pres = NA, back = NA)
    for (z in 1:k) {
      df_res$pres[z] <- length(data_pts$presence[data_pts$presence == 1 & cell_vals == z])
      df_res$back[z] <- length(data_pts$presence[data_pts$presence == 0 & cell_vals == z])
    }
    folds[[i]] <- cell_vals
    result[[i]] <- df_res
  }
  
  # Get the best one, based on the more even distribution of presence points
  best <- which.min(unlist(
    lapply(result, function(x){
      sd(x$pres)
    })
  ))
  
  # Generate a records table
  n_points <- data.frame(train_pres = rep(NA, k), train_back = NA,
                         test_pres = NA, test_back = NA)
  
  bblock <- result[[best]]
  
  for (z in 1:k) {
    n_points$train_pres[z] <- sum(bblock$pres[bblock$fold != z])
    n_points$train_back[z] <- sum(bblock$back[bblock$fold != z])
    n_points$test_pres[z] <- sum(bblock$pres[bblock$fold == z])
    n_points$test_back[z] <- sum(bblock$back[bblock$fold == z])
  }
  
  # Return object
  return(list(
    folds = folds[[best]],
    k = k,
    records = n_points
  ))
  
}


#' Plot CV folds
#'
#' Plot the cross-validation folds produced with [mp_prepare_blocks()].
#'
#' @param sdm_data an object of type \code{sdm_dat} produced with
#'   [mp_prepare_data()] and with cross-validation blocks assigned using
#'   [mp_prepare_blocks()]
#' @param block_type which block type should be ploted (one single name). If
#'   NULL the first one will be ploted
#' @param facet_run if \code{TRUE}, then facets equal the number of folds will
#'   be ploted, each showing which records will be used for training and which
#'   for testing
#' @param presence_only if \code{TRUE}, then only the presence points are ploted
#' @param include_map if \code{TRUE}, a base map of the continents is added for
#'   reference. The package [rnaturalearth] is needed for this
#'
#' @return nothing
#' @export
#'
#' @examples
#' \dontrun{
#' plot_folds(sdm_data)
#' }
plot_folds <- function(sdm_data, block_type = NULL,
                       facet_run = FALSE, presence_only = FALSE,
                       include_map = FALSE) {
  if (is.null(block_type)) {
    block_type <- names(sdm_data$blocks$folds)[1]
  }
  pdat <- cbind(sdm_data$coord_training, folds = sdm_data$blocks$folds[[block_type]])
  
  if (presence_only) {
    pdat <- pdat[sdm_data$training$presence == 1,]
  }
  
  if (include_map) {
    if (require("rnaturalearth", quietly = T)) {
      base_map <- rnaturalearth::ne_countries(returnclass = "sf")
    } else {
      cat("Package rnaturalearth is not available. Ignoring 'include_map'. \n")
      include_map <- FALSE
    }
  }
  
  if (facet_run) {
    for (z in 1:max(pdat$folds)) {
      ndat <- pdat
      ndat$folds <- ifelse(ndat$folds == z, "testing", "training")
      ndat$run <- z
      if (z == 1) {
        fdat <- ndat
      } else {
        fdat <- rbind(fdat, ndat)
      }
    }
    #pdat$folds <- as.factor(pdat$folds)
    p <- ggplot2::ggplot(fdat) +
      ggplot2::geom_point(ggplot2::aes(x = decimalLongitude, y = decimalLatitude, color = folds), size = 0.5) +
      ggplot2::scale_color_manual(values = c("#a8ddb5", "#2c7fb8")) +
      ggplot2::coord_equal() +
      ggplot2::theme_light() +
      ggplot2::facet_wrap(~run)
  } else {
    pdat$folds <- as.factor(pdat$folds)
    p <- ggplot2::ggplot(pdat) +
      ggplot2::geom_point(ggplot2::aes(x = decimalLongitude, y = decimalLatitude, color = folds), size = 0.5) +
      ggplot2::scale_color_brewer(palette = "YlGnBu") +
      ggplot2::coord_equal() +
      ggplot2::theme_light()
  }
  if (include_map) {
    p <- p + ggplot2::geom_sf(data = base_map, fill = "gray80", alpha = .7) +
      ggplot2::coord_sf() +
      ggplot2::theme_void() +
      ggplot2::xlab(NULL) + ggplot2::ylab(NULL)
  }
  print(p)
  return(invisible(NULL))
}


#' Get information on the number of points per fold
#'
#' @param sdm_data an sdm_dat object prepared using [mp_prepare_data()] and with 
#' assigned blocks using [mp_prepare_blocks()]
#' @param block_type the block type
#'
#' @return a data.frame containing number of points per class and per fold
#' @export
#'
#' @examples 
#' \dontrun{
#' folds_info(sdm_data)
#' }
folds_info <- function(sdm_data, block_type = "spatial_grid") {
  
  k <- sdm_data$blocks$n
  
  n_points <- data.frame(train_pres = rep(NA, k), train_back = NA,
                         test_pres = NA, test_back = NA)
  
  bblock <- list(fold = sdm_data$blocks$folds[[block_type]],
                 pres = sdm_data$training$presence)
  
  for (z in 1:k) {
    n_points$train_pres[z] <- length(bblock$pres[bblock$fold != z & bblock$pres == 1])
    n_points$train_back[z] <- length(bblock$pres[bblock$fold != z & bblock$pres == 0])
    n_points$test_pres[z] <- length(bblock$pres[bblock$fold == z & bblock$pres == 1])
    n_points$test_back[z] <- length(bblock$pres[bblock$fold == z & bblock$pres == 0])
  }
  
  return(n_points)
}


#' Get grid for spatial block CV
#'
#' This function uses the function [blockCV::cv_spatial_autocor()] from the package
#' [blockCV] to calculate the appropriate size of the grid
#'
#' @param sdm_data object of type `sdm_dat` with species occurrence records
#' @param env_layers environmental layers in SpatRaster format
#' @param sel_vars which variables should be used for calculating the block size
#' @param min_block_size minimum size of the block. If the calculated range is
#'   smaller than this threshold, then the minimum size is used instead
#' @param max_block_size maximum size of the block. If the calculated range is
#'   larger than this threshold, then the maximum size is used instead
#' @param verbose print messages
#'
#' @return grid (SpatRaster)
#' @export
#'
#' @examples
#' \dontrun{
#' get_block_grid(sp_data, env_data)
#' }
get_block_grid <- function(sdm_data, env_layers, sel_vars = NULL,
                           min_block_size = 0.2, max_block_size = 20,
                           verbose = FALSE) {
  
  f_data <- cbind(sdm_data$training, sdm_data$coord_training)
  f_data <- sf::st_as_sf(f_data, coords = colnames(sdm_data$coord_training),
                         crs = "EPSG:4326")
  
  if (nrow(f_data) > 10000) {
    f_data <- f_data[sample(1:nrow(f_data), 5000),]
  }
  
  vars <- colnames(sdm_data$training)
  vars <- vars[!grepl("presence", vars)]
  
  if (!is.null(sel_vars)) {
    vars <- vars[vars %in% sel_vars]
  }
  
  if (verbose) cli::cli_alert_info("Calculating grid size...")
  sf_status <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  size <- suppressMessages(suppressWarnings(
    blockCV::cv_spatial_autocor(x = f_data,
                                column = vars,
                                plot = FALSE,
                                progress = FALSE)
  ))
  sf::sf_use_s2(sf_status)
  
  lyrs_ext <- ext(env_layers)
  sel_size <- (size$range/1000)/111
  
  if (sel_size > max_block_size) {
    if (verbose) cli::cli_alert_warning("Grid size ({sel_size}) is larger than {.var max_block_size}. Using {max_block_size} instead.")
    sel_size <- max_block_size
  }
  if (sel_size < min_block_size) {
    if (verbose) cli::cli_alert_warning("Grid size ({sel_size}) is smaller than {.var min_block_size}. Using {min_block_size} instead.")
    sel_size <- min_block_size
  }
  
  xmin_ext <- round(lyrs_ext[1]-0.1, 1)
  ymax_ext <- round(lyrs_ext[4]+0.1, 1)
  
  ymin_t <- round(lyrs_ext[3]-0.1, 1)
  test_ymin <- seq(ymax_ext, ymin_t, by = -sel_size)
  ymin_ext <- ifelse(min(test_ymin) > ymin_t, round((min(test_ymin) - sel_size), 1), min(test_ymin))
  
  xmax_t <- round(lyrs_ext[2]+0.1, 1)
  test_xmax <- seq(xmin_ext, xmax_t, by = sel_size)
  xmax_ext <- ifelse(max(test_xmax) < xmax_t, round((max(test_xmax) + sel_size), 1), max(test_xmax))
  
  block_list <- list(spatial_grid = rast(ext(xmin_ext, xmax_ext, ymin_ext, ymax_ext), resolution = sel_size))
  
  return(block_list)
}



# Old version
# mp_prepare_blocks <- function(sdm_data, method = "blockcv", 
#                               block_types = c("spatial_grid", "spatial_lat", "random"),
#                               nfolds = 5, auto_range = FALSE, range = NULL, 
#                               lat_blocks = 50, manual_shp = NULL, verbose = TRUE) {
#   
#   if (class(sdm_data)[1] != "sdm_dat") {
#     cli::cli_abort("Supply an object of class sdm_dat produced with mp_prepare_data.")
#   }
#   
#   if (!method %in% c("manual", "blockcv")) {
#     cli::cli_abort("Method should be blockcv or manual.")
#   }
#   
#   if (!any(block_types %in% c("spatial_grid", "spatial_lat", "random"))) {
#     cli::cli_abort("block_types should be one (or more) of spatial_grid, spatial_lat, random.")
#   }
#   
#   folds <- list()
#   
#   sf_pts <- sf::st_as_sf(cbind(sdm_data$training,
#                                sdm_data$coord_training), coords = c("decimalLongitude", "decimalLatitude"),
#                          crs = 4326)
#   
#   if (method == "blockcv") {
#     
#     if (auto_range) {
#       sf::sf_use_s2(FALSE)
#       auto_cor <- blockCV::cv_spatial_autocor(x = sf::st_make_valid(sf_pts), column = "presence")
#       range = auto_cor$range
#     } else {
#       sf::sf_use_s2(FALSE)
#       if (is.null(range)) {
#         if (verbose) cli::cli_alert_danger("range was not supplied. Using 20km as range.")
#         range <- 20000
#       }
#     }
#     
#     # Improve to ensure all species will be equally distributed in blocks!
#     if ("spatial_grid" %in% block_types) {
#       if (verbose) cli::cli_alert_info("Generating grid blocks.")
#       blocks <- blockCV::cv_spatial(sf_pts, column = "presence", iteration = 100, size = range,
#                                     hexagon = FALSE, k = nfolds, progress = F, report = T, plot = F, selection = "systematic")
#       folds <- c(folds, spatial_grid = list(blocks$folds_ids))
#     }
#     
#     if ("spatial_lat" %in% block_types) {
#       if (verbose) cli::cli_alert_info("Generating latitudinal blocks.")
#       blocks <- blockCV::cv_spatial(sf_pts, column = "presence", iteration = 100, rows_cols = c(lat_blocks, 1),
#                                     hexagon = FALSE, k = nfolds, progress = F, report = T, plot = F, extend = 0.5)
#       folds <- c(folds, spatial_lat = list(blocks$folds_ids))
#     }
#     
#     if ("random" %in% block_types) {
#       if (verbose) cli::cli_alert_info("Generating random samples.")
#       blocks <- sample(1:5, nrow(sf_pts), replace = T)
#       folds <- c(folds, random = list(blocks))
#     }
#     
#   } else {
#     
#     if ("spatial_grid" %in% block_types) {
#       if (verbose) cli::cli_alert_info("Generating grid blocks.")
#       blocks <- blockCV::cv_spatial(sf_pts, column = "presence", iteration = 200,
#                                     user_blocks = manual_shp$spatial_grid,
#                                     k = nfolds, progress = F, report = T, plot = F)
#       folds <- c(folds, spatial_grid = list(blocks$folds_ids))
#     }
#     
#     if ("spatial_lat" %in% block_types) {
#       if (verbose) cli::cli_alert_info("Generating latitudinal blocks.")
#       blocks <- blockCV::cv_spatial(sf_pts, column = "presence", iteration = 200, 
#                                     user_blocks = manual_shp$spatial_lat,
#                                     k = nfolds, progress = F, report = T, plot = F)
#       folds <- c(folds, spatial_lat = list(blocks$folds_ids))
#     }
#     
#     if ("random" %in% block_types) {
#       if (verbose) cli::cli_alert_info("Generating random samples.")
#       blocks <- sample(1:5, nrow(sf_pts), replace = T)
#       folds <- c(folds, random = list(blocks))
#     }
#     
#   }
#   
#   sdm_data$blocks <- list(n = nfolds,
#                           folds = folds,
#                           gen_method = method)
#   
#   return(sdm_data)
# }
