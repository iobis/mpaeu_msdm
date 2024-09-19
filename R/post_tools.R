#' Utilities for post-processing results
#'
#' @param taxon taxon ID (AphiaID) for the target species
#' @param what what you need (for example, "cvmetrics" will return files of
#'   cross-validation metrics). If `NULL` or `NA` it will return all available
#'   files
#' @param type optionally, provide the file type (example, "parquet")
#' @param results_path path to the results folder
#' @param model an optional character name for the model acronym. If not
#'   provided it will use the most recent
#' @param file path to a parquet file
#' @param to target format. Can be one of "csv", "txt" or "rds"
#' @param new_path an optional path for the new file. If NULL, it will be saved
#'   in the same path. Note that you should provide only the folder path, as the
#'   file name will be exactly the same as the original file, but with the new
#'   format
#' @param raster_obj a SpatRaster object or a path for a raster file (any format
#'   accepted by [terra::rast()])
#' @param sel_band a `character` with the band name or a `numeric` with the
#'   number of the band
#' @param path an optional file path to save the individual band. This should be
#'   a full path, including the file format (any accepted by
#'   [terra::writeRaster()]). If `NULL`, returns the SpatRaster
#' @param ... further arguments to [terra::writeRaster()]
#' 
#' @description
#' `find_file()` helps to find a specific file in the results folder for a
#' specific taxon
#'
#' `parquet_to()` is a wrapper to quickly convert `parquet` files to other
#' common formats
#'
#' `deband()` provides an easy way to extract a single band from multiband
#' rasters
#'
#' @return `find_file()` returns a character vector with file paths (full
#' address relative to the working directory) 
#' `parquet_to()` returns nothing, save the file 
#' `deband()` returns the SpatRaster if `path=NULL`, or nothing and saves the file
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' files <- find_file(127326, "cvmetrics", type = "parquet")
#' 
#' parquet_to(files[1])
#' 
#' mask_file <- find_file(127326, "mask")
#' 
#' single_band <- deband(mask_file, sel_band = "fit_region")
#' }
find_file <- function(taxon, what = NULL,
                      type = NULL,
                      results_path = "results",
                      model = NULL) {
  
  av_files <- list.files(glue::glue("../mpaeu_sdm/{results_path}/taxonid={taxon}"),
                         recursive = TRUE, full.names = TRUE)
  
  if (length(av_files) < 1) {
    cli::cli_abort("No files found for taxon {.val {taxon}} on folder {.path {results_path}}")
  } else {
    if (is.null(model)) {
      av_models <- list.files(glue::glue("../mpaeu_sdm/{results_path}/taxonid={taxon}"),
                              full.names = T)
      av_models_info <- fs::file_info(av_models)
      av_models_info <- av_models_info[order(av_models_info$modification_time,
                                             decreasing = T),]
      model <- av_models_info$path[1]
      model <- basename(model)
      model <- gsub("model=", "", model)
    }
    
    av_files <- av_files[grepl(paste0("model=", model), av_files)]
    
    if (!is.null(type)) {
      av_files <- av_files[grepl(type, av_files)]
    }
    
    if (is.null(what)) {
      what <- NA
    }
    if (is.na(what)) {
      fp <- av_files
    } else {
      fp <- av_files[grepl(what, av_files)]
      
      if (length(fp) < 1) {
        cli::cli_abort("No files found for `what` = '{what}', but there are {.val {length(av_files)}} files available.
                       Set `what` as NA or NULL to list all.")
      }
    }
    
  }
  
  return(fp)
  
}

#' @rdname find_file
#' @export
parquet_to <- function(file, to = "csv", new_path = NULL) {
  
  f <- arrow::read_parquet(file)
  
  if (!is.null(new_path)) {
    file <- basename(file)
    file <- file.path(new_path, file)
  }
  
  if (to == "csv") {
    nf <- gsub("parquet", "csv", file)
    write.csv(f, nf, row.names = F)
  }
  if (to == "txt") {
    nf <- gsub("parquet", "txt", file)
    write.table(f, nf, row.names = F)
  }
  if (to == "rds") {
    nf <- gsub("parquet", "rds", file)
    saveRDS(f, file = nf)
  }
  cli::cli_alert_success("File saved as {.val {to}} on {.path {nf}}")
  
  return(invisible(NULL))
}

#' @rdname find_file
#' @export
deband <- function(raster_obj, sel_band, path = NULL, ...) {
  if (is.character(raster_obj)) {
    raster_obj <- terra::rast(raster_obj)
  }
  if (is.numeric(sel_band)) {
    if (sel_band > nlyr(raster_obj)) {
      cli::cli_abort("There are {.val {nlyr(raster_obj)}}, but you selected band {.val {sel_band}}. Check your selection.")
    } else {
      band_res <- raster_obj[[sel_band]]
    }
  } else {
    band_names <- names(raster_obj)
    if (!sel_band %in% band_names) {
      cli::cli_abort("Band name non available. You selected {.val {sel_band}} but available values are: {.val {band_names}}")
    } else {
      band_res <- raster_obj[[sel_band]]
    }
  }
  
  if (save_file && !is.null(path)) {
    writeRaster(band_res, filename = path, ...)
    return(invisible(NULL))
  } else {
    return(band_res)
  }
  
}