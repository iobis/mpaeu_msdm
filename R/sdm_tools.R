#' Convert points to 1 point per cell
#'
#' @param pts dataframe with two columns named decimalLongitude and decimalLatitude
#' @param base_rast an optional raster (SpatRaster format) to be used as a model.
#' If supplied, the results are masked according to this raster, and any point that falls
#' outside a valid area will be removed.
#' @param res if \code{base_rast} is not supplied, then a raster with the given resolution
#' (in degrees) is used as a base
#'
#' @return data.frame with points
#' @export
#'
#' @examples
#' \dontrun{
#' as_cell(all_pts)
#' }
as_cell <- function(pts, base_rast = NULL, res = 0.05) {
  
  coords <- c("decimalLongitude", "decimalLatitude")
  
  if (!all(coords %in% colnames(pts))) {
    cli::cli_abort("{.var pts} should have columns named {.var {coords}}")
  }
  
  if (!is.null(base_rast) & class(base_rast)[1] != "SpatRaster") {
    cli::cli_abort("{.var base_rast} should be of type {.var SpatRaster}")
  }
  
  if (is.null(base_rast)) {
    base <- rast(res = res)
  } else {
    base <- base_rast[[1]]
  }
  
  base[] <- NA
  base[cellFromXY(base, as.data.frame(pts[,coords]))] <- 1
  
  if (!is.null(base_rast)) {
    base <- mask(base, base_rast[[1]])
  }
  
  pts_cell <- as.data.frame(base, xy = T)[,1:2]
  
  colnames(pts_cell) <- coords
  
  return(pts_cell)
  
}



#' Plot a Leaflet interactive map
#'
#' @param layer a SpatRaster with the prediction or environmental layer.
#' This can have one or more layers, preferentially named (otherwise, the
#' function will assign names as 'Layer 1', 'Layer 2', and so on)
#' @param pts an optional data.frame containing the geographic coordinates for points (two columns, x and y order)
#' @param map_palette the color palette to be used for plotting. See [leaflet::colorNumeric()] for more details
#' @param points_color the color to be used fot the points, if supplied
#'
#' @return nothing, ploted map
#' @export
#'
#' @examples
#' #' \dontrun{
#' plot_leaflet(sdm_prediction)
#' }
plot_leaflet <- function(layer, pts = NULL,
                         map_palette = "Blues", points_color = "yellow") {
  
  base_map <- leaflet::leaflet() %>% 
    leaflet::addProviderTiles("OpenStreetMap.Mapnik")
  
  if (length(unique(names(layer))) != length(names(layer))) {
    names(layer) <- paste("Layer", 1:nlyr(layer))
  }
  
  if (nlyr(layer) > 1) {
    eval(parse(text = paste0(
      "base_map <- base_map %>%",
      paste0("leaflet::addRasterImage(layer[[",1:nlyr(layer),"]], colors = '", map_palette,"', group = '", names(layer),"')", collapse = "%>%")
    )))
  }
  
  if (!is.null(pts)) {
    base_map <- base_map %>%
      leaflet::addCircles(lng = pts[, 1],
                          lat = pts[, 2],
                          color = points_color,
                          group = "Points") %>%
      leaflet::addLayersControl(
        baseGroups = names(layer),
        overlayGroups = c("Points")
      )
  } else {
    base_map <- base_map %>%
      leaflet::addLayersControl(
        baseGroups = names(layer)
      )
  }
  
  print(base_map)
  
  return(invisible(NULL))
}




#' Convert GeoTIFF files to Cloud Optimized GeoTIFF format (COG)
#'
#' @param id path to the GeoTIFF file to be optimized
#' @param remove_files if `TRUE` (default), the original file is removed if the
#'   optimization succeeds or the COG is removed if the optimization failed
#'   somehow
#' @param verbose if `TRUE` output messages
#'
#' @return optimization status and saved files
#' @export
#'
#' @details
#' There are three possible status returned by the function:
#' - optimization-succeeded: optimization was successful and file is saved. 
#' The non-optimized file is removed if `remove_files = TRUE`.
#' - optimization-failed: optimization failed.
#' - optimization-gen-failed: file was probably generated, but by some reason 
#' it was not correctly optimized. This was included more for control, as such a
#'  case should be extremely rare.
#'  
#' This function depends on the [rio-cogeo plugin](https://cogeotiff.github.io/rio-cogeo/), which
#' is a Python program with CLI interface. To install the plugin, you should use:
#' 
#' ``` 
#' $ pip install rio-cogeo
#' ```
#' 
#' @examples
#' \dontrun{
#' cogeo_optim("test_file.tif")
#' }
cogeo_optim <- function(id, remove_files = TRUE, verbose = FALSE) {
  
  if (verbose) cli::cli_alert_info("Optimizing {.path {id}} to COG using {.var rio}")
  outid <- gsub("\\.tif", "_cog.tif", id)
  
  out <- sys::exec_internal("rio",
                            c("cogeo", "create", id, outid),
                            error = F)
  if (out$status == 1) {
    if (verbose) cli::cli_alert_danger("Optimization failed")
    to_return <- data.frame(file = id, status = "optimization-failed")
  } else {
    # Validate
    out <- sys::exec_internal("rio",
                              c("cogeo", "validate", outid),
                              error = F)
    if (out$status == 1) {
      to_return <- data.frame(file = id, status = "optimization-failed")
      if (file.exists(outid) & remove_files) file.remove(outid)
      if (verbose) cli::cli_alert_danger("Optimization failed")
    } else {
      out <- sys::as_text(out$stdout)
      if (grepl("is a valid cloud optimized GeoTIFF", out)) {
        to_return <- data.frame(file = id, status = "optimization-succeeded")
        if (verbose) cli::cli_alert_success("Optimization succeeded. Output file: {.path {outid}}")
        if (remove_files) file.remove(id)
      } else {
        to_return <- data.frame(file = id, status = "optimization-gen-failed")
        if (file.exists(outid) & remove_files) file.remove(outid)
        if (verbose) cli::cli_alert_danger("Optimization is invalid")
      }
    }
  }
  return(to_return)
}


#' Generate log of modeling
#'
#' @param algos algorithms that will be used to fit
#' 
#' @details
#' The following parameters are available in the list object:
#' - taxonID
#' - scientificName
#' - group
#' - model_date
#' - model_acro
#' - n_init_points = input number of initial points before autocorrelation
#' - model_fit_points 
#' - model_eval_points
#' - algorithms
#' - algorithms_parameters = list with parameters used for model fitting (auto filled)
#' - model_result = result of model fitting
#' - model_bestparams = the best fitting parameters
#' - model_details = list with the following parameters:
#'    - ecoregions
#'    - ecoregions_included 
#'    - limited_by_depth = TRUE or FALSE, if it was limited by depth
#'    - depth_buffer
#'    - block_size
#'    - background_size
#'    - control_bias
#'    - hypothesis_tested
#'    - best_hypothesis
#'    - variables
#' - model_good = which models performed well according to the metric and threshold
#' - model_good_metric = metric used for defining "good" models
#' - model_good_threshold = threshold used for defining "good" models
#' - model_posteval = list with post-evaluation results
#' - timings
#' - obissdm_version = auto filled, version of the `obissdm` package
#'
#' Other parameters can be added as needed.
#'
#' @return a list
#' @export
#' 
#' @seealso [save_log()]
#'
#' @examples
#' sdm_log <- gen_log(c("maxent", "brt"))
#' \dontrun{
#' save_log(sdm_log, "species_log.json")
#' }
gen_log <- function(algos){
  
  log_obj <- list(
    taxonID = NULL,
    scientificName = NULL,
    group = NULL,
    model_date = NULL,
    model_acro = NULL,
    n_init_points = NULL,
    model_fit_points = NULL,
    model_eval_points = NULL,
    algorithms = algos,
    algorithms_parameters = eval(parse(text = paste0(
      "list(", paste0(algos, "= unclass(obissdm::sdm_options('", algos, "'))", collapse = ","), ")"
    ))),
    model_result = eval(parse(text = paste0(
      "list(", paste0(c(algos, "ensemble"), "=", "NULL", collapse = ","), ")"
    ))),
    model_bestparams = eval(parse(text = paste0(
      "list(", paste0(algos, "=", "NULL", collapse = ","), ")"
    ))),
    model_details = list(
      ecoregions = NULL,
      ecoregions_included = NULL,
      limited_by_depth = NULL,
      depth_buffer = NULL,
      block_size = NULL,
      background_size = NULL,
      control_bias = NULL,
      hypothesis_tested = NULL,
      best_hypothesis = NULL,
      variables = NULL
    ),
    model_good = NULL,
    model_good_metric = NULL,
    model_good_threshold = NULL,
    model_posteval = eval(parse(text = paste0(
      "list(", paste0(c(algos, "ensemble", "niche", "hyperniche"), "=", "NULL", collapse = ","), ")"
    ))),
    timings = NULL,
    obissdm_version = as.character(packageVersion("obissdm"))
  )
  
  return(log_obj)
}

#' Save log object
#'
#' @param log_object log object generated with [gen_log()]
#' @param file_path the path to save the file in .json format
#'
#' @return saved file
#' @export
#'
#' @examples
#' \dontrun{
#' save_log(sdm_log, "species_log.json")
#' }
save_log <- function(log_object, file_path) {
  
  jsonlite::write_json(log_object, file_path, pretty = T)
  
  return(invisible(NULL))
  
}