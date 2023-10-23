############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################# SDM - data preparing #############################

#' Prepare data for use in SDM
#'
#' @param training a data frame containing xy coordinates  (and optionally a
#'   column named presence with the presence-absence status) for the training
#'   points (coordinates should be named decimalLongitude and decimalLatitude)
#' @param eval_data an optional data frame containing xy coordinates (and
#'   optionally a column named presence with the presence-absence status) of the
#'   points used for evaluation
#' @param species_id a name or id for the species data
#' @param native_shp an optional shapefile with the native or expert range data
#' @param env_layers environmental layers in terra rast format
#' @param gen_quad if \code{TRUE}, quadrature points are generated
#' @param quad_number the number of quadrature points to be generated. If set to
#'   \code{NULL}, then 10000 points are generated
#' @param pred_quad optionally, you can supply a data.frame of quadrature points
#'  containing the xy coordinates (named "decimalLongitude" and "decimalLatitude"),
#'  a column named presence with values of 0, and the values for the variables.
#'  Those should be named exactly as the \code{env_layers} and in the same order.
#'  Supplying the pre-defined quadrature points can speed-up the process when several
#'  runs will take place. If this is supplied, \code{gen_quad} is ignored
#' @param verbose should messages be printed
#' 
#' @return an object of class sdm_dat containing the data for SDMs
#' @export
#'
#' @examples
#' \dontrun{
#'   mp_prepare_data(training_data, env_layers = env)
#' }
mp_prepare_data <- function(training, eval_data = NULL, species_id, native_shp = NULL, env_layers,
                            gen_quad = TRUE, quad_number = NULL, pred_quad = NULL, verbose = TRUE) {
  
  if (missing(species_id)) {
    cli::cli_abort("A species_id must be supplied")
  }
  
  if (!is.null(pred_quad)) {
    gen_quad <- FALSE
  }
  
  sdm_dat <- list()
  
  training <- .check_colnames(training)
  
  training <- training[,c("decimalLongitude", "decimalLatitude")]
  
  training_data <- cbind(training, presence = 1, terra::extract(env_layers, training, ID = F))
  
  if (gen_quad) {
    if (is.null(quad_number)) {
      quad_number <- 10000
      if (verbose) cli::cli_alert_info("Number of quadrature points not supplied. Using 10000.")
    }
    
    # quad_data <- spatSample(env_layers, size = quad_number,
    #                         na.rm = T, xy = T)
    xy_d <- as.data.frame(env_layers, xy = T, na.rm = T)[,1:2]
    colnames(xy_d) <- c("decimalLongitude", "decimalLatitude")
    xy_d <- xy_d[sample(1:nrow(xy_d), size = quad_number, replace = FALSE),]
    quad_data <- cbind(xy_d, presence = 0, terra::extract(env_layers, xy_d, ID = F))
    
    training_data <- rbind(training_data, quad_data)
    
    if (any(is.na(training_data))) {
      nas <- apply(training_data, 1, sum)
      if (verbose) cli::cli_alert_warning("{sum(is.na(nas))} NA point{?s} found and removed from training data.")
      training_data <- training_data[!is.na(nas), ]
    }
  }
  
  if (!is.null(pred_quad)) {
    if (!isTRUE(all.equal(colnames(training_data), colnames(pred_quad)))) {
      cli::cli_abort("Column names of pred_quad should be exactly {.var {colnames(training_data)}}")
    }
    
    training_data <- rbind(training_data, pred_quad)
    
    if (any(is.na(training_data))) {
      nas <- apply(training_data, 1, sum)
      if (verbose) cli::cli_alert_warning("{sum(is.na(nas))} NA point{?s} found and removed from training data.")
      training_data <- training_data[!is.na(nas), ]
    }
    
    quad_number <- nrow(pred_quad)
  }
  
  if (!is.null(eval_data)) {
    
    eval_data <- .check_colnames(eval_data)
    
    if (!"presence" %in% colnames(eval_data)) {
      eval_data$presence <- 1
    }
    
    coord_eval <- eval_data[,c("decimalLongitude", "decimalLatitude")]
    
    eval_data <- as.data.frame(cbind(presence = eval_data[,"presence"],
                                     terra::extract(env_layers, coord_eval, ID = F)))
    
    if (any(is.na(eval_data))) {
      nas <- apply(eval_data, 1, sum)
      if (verbose) cli::cli_alert_warning("{sum(is.na(nas))} NA point{?s} found and removed from evaluation data.")
      training_data <- training_data[!is.na(nas), ]
    }
    
  } else {
    eval_data <- coord_eval <- NULL
  }
  
  if (!is.null(native_shp)) {
    if (class(native_shp)[1] != "sf") {
      native_shp <- st_as_sf(native_shp)
    }
  }
  
  sdm_dat <- list(training = training_data[, !colnames(training_data) %in% c("decimalLongitude", "decimalLatitude")],
                  coord_training = training_data[, c("decimalLongitude", "decimalLatitude")],
                  eval_data = eval_data,
                  coord_eval = coord_eval,
                  sp_range = native_shp,
                  quad_number = quad_number,
                  species = species_id)
  
  class(sdm_dat) <- c("sdm_dat", class(sdm_dat))
  
  return(sdm_dat)
}

.check_colnames <- function(data){
  if (!any(grepl("decimallongitude|decimalLongitude", colnames(data), ignore.case = T))) {
    stop("Impossible to identify a longitude column in the data.")
  }
  if (!any(grepl("decimallatitude|decimalLatitude", colnames(data), ignore.case = T))) {
    stop("Impossible to identify a latitude column in the data.")
  }
  
  colnames(data)[grepl("decimallongitude|decimalLongitude", colnames(data),
                       ignore.case = T)] <- "decimalLongitude"
  
  colnames(data)[grepl("decimallatitude|decimalLatitude", colnames(data),
                       ignore.case = T)] <- "decimalLatitude"
  
  return(data)
}

#' @export
print.sdm_dat <- function(x, ...) {
  
  cli::cli_h1("SDM data for species {.emph {x$species}}")
  
  cli::cli_h3("Training data")
  
  covars <- colnames(x$training)
  covars <- covars[!grepl("decimalLatitude|decimalLongitude|presence", covars)]
  
  cli::cli_bullets(c(">" = "Presence points: {nrow(x$training[x$training$presence == 1,])}",
                     ">" = "Absence points: {nrow(x$training[x$training$presence == 0,])}",
                     ">" = "Number of variables: {length(covars)}",
                     ">" = "Variable{?s}: {covars}"))
  
  cli::cli_h3("Evaluation data")
  
  if (!is.null(x$eval_data)) {
    cli::cli_bullets(c(">" = "Presence points: {nrow(x$eval_data[x$eval_data$presence == 1,])}",
                       ">" = "Absence points: {nrow(x$eval_data[x$eval_data$presence == 0,])}"))
  } else {
    cli::cli_alert_warning("Evaluation data not supplied.")
  }
  
  cli::cli_h3("Quadrature points")
  
  if (!is.null(x$quad_number)) {
    cli::cli_alert_info("{x$quad_number} quadrature points generated.")
  } else {
    cli::cli_alert_info("No quadrature points were generated.")
  }
  
  cli::cli_h3("Native range shapefile")
  
  if (!is.null(x$sp_range)) {
    cli::cli_alert_info("Native range supplied with bounding box {round(sf::st_bbox(x$sp_range), 2)}")
  } else {
    cli::cli_alert_info("No native range shapefile was supplied.")
  }
  
  if (!is.null(x$blocks)) {
    cli::cli_h3("Cross-validation blocks")
    cli::cli_bullets(c(">" = "Number of blocks: {x$blocks$n}",
                       ">" = "Types of blocks: {names(x$blocks$folds)}"))
  }
  
  return(invisible(NULL))
}


