#' Download data from Bio-ORACLE version 3
#'
#' @param datasets an optional character vector of datasets
#' @param future_scenarios an optional vector of future scenarios
#' @param time_steps an optional named list of time steps (see details)
#' @param variables an optional vector of the variants that you want (e.g. mean, min, max)
#' @param terrain_vars an optional vector of terrain variables
#' @param outdir the output directory
#' @param skip_exist if \code{TRUE} (default) the function will skip existent files
#' @param keep_raw wether or not should the folder with raw files be kept (recommended is FALSE to save space)
#' @param verbose print messages to track download
#'
#' @return saved files
#' @export
#' 
#' @details
#' Environmental layers are downloaded from the Bio-ORACLE ERDDAP server
#' To see all available options, visit:
#' https://erddap.bio-oracle.org/erddap/griddap/index.html?page=1&itemsPerPage=1000
#' 
#' We use the biooracler package which is still not on CRAN (but soon will be). If
#' you run the function and does not have it installed, the function will print an
#' error with the GitHub address to install it.
#' 
#' The files are saved on the following structure
#' data
#' \_ env: environmental layers
#'      \_ raw
#'          \_ env_layers: raw files downloaded from Bio-ORACLE - deleted at the end, except if \code{keep_raw = TRUE}
#'      \_ terrain: terrain layers (static)
#'      \_ current: current period layers
#'      \_ future: future period scenarios folders
#'           \_ ssp1: future layers for SSP1 scenario
#'           \_ ...: other scenarios
#'           
#' Time steps should be a named list, for each period a vector containing start and end periods (for a single step, use same value for both).
#' Datasets names should be the same one got from the ERDDAP list (see link above).
#' If a variable is not available for a period/depth but is for others, just supply it and the function will trow an error for those that it's unable to download but keep with the others.
#' See examples for details on how to set up the function.
#' 
#' Note: in the special cases that a dataset have "mean" in its name (e.g. PAR_mean), it will be renamed
#' and the _mean suppressed to avoid problems with other functions.
#'
#' @examples
#' \dontrun{
#' # List datasets to download ----
#' datasets <- c(
#'   "thetao_baseline_2000_2019_depthsurf",
#'   "so_baseline_2000_2019_depthsurf",
#'   "PAR_mean_baseline_2000_2020_depthsurf",
#'   "po4_baseline_2000_2018_depthsurf",
#'   "phyc_baseline_2000_2020_depthsurf",
#'   "sws_baseline_2000_2019_depthsurf",
#'   "siconc_baseline_2000_2020_depthsurf"
#' )
#' 
#' datasets <- c(datasets,
#'               gsub("depthsurf", "depthmean", datasets),
#'               gsub("depthsurf", "depthmax", datasets))
#' 
#' # List scenarios to download ----
#' future_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585")
#' 
#' # Define time steps ----
#' time_steps <- list(
#'   current = c("2010-01-01", "2010-01-01"),
#'   dec50 = c("2050-01-01", "2050-01-01"),
#'   dec100 = c("2090-01-01", "2090-01-01")
#' )
#' 
#' # Define variables to be download
#' # In general, available are: min, mean, max, range, ltmin, ltmax, and sd
#' vars <- c("min", "mean", "max")
#' 
#' get_env_data(datasets = datasets, future_scenarios = future_scenarios,
#'              time_steps = time_steps, variables = vars,
#'              terrain_vars = "bathymetry_mean")
#' 
#' }
get_env_data <- function(datasets = NULL,
                         future_scenarios = NULL,
                         time_steps = NULL,
                         variables = c("min", "mean", "max"),
                         terrain_vars = NULL,
                         outdir = "data/env/",
                         skip_exist = TRUE,
                         keep_raw = FALSE,
                         verbose = TRUE) {
  
  is_installed <- try(find.package("biooracler"), silent = TRUE)
  
  if (class(is_installed) == "try-error") {
    code_var <- cli::code_highlight("devtools::install_github('bio-oracle/biooracler')")
    cli::cli_abort("Package {.pkg biooracler} is not installed. Check if it's already available on CRAN or use the following code to install from GitHub: {.code {code_var}}")
  } else {
    library(biooracler)
  }
  
  # Define out directory
  rawdir <- paste0(outdir, "raw/env_layers")
  fs::dir_create(rawdir)
  fs::dir_create(paste0(outdir, c("current", paste0("future/", future_scenarios), "terrain")))
  
  # We create a function to check if the file is already downloaded
  # If you want to ignore and download again, just set ignore = T
  # in the loop functions
  check_file_exist <- function(file_name, expres, ignore = FALSE) {
    if (!file.exists(file_name) | ignore) {
      expres
    } else {
      if (verbose) cli::cli_alert_info("{file_name} already exists.")
    }
  }
  
  
  # Download all files ----
  
  if (!is.null(datasets)) {
    # We create a loop that will go through all the datasets
    for (id in datasets) {
      if (verbose) cli::cli_alert_info("Downloading {id}")
      
      # Generate the correct name for the variable, e.g. thethao_mean
      ds_vars <- gsub("_baseline(.*)$", "", id)
      ds_vars <- paste(ds_vars, variables, sep = "_")
      
      # Run for each time step
      for (ts in 1:length(time_steps)) {
        
        # Get name of the period (e.g. current)
        period <- names(time_steps)[ts]
        
        if (period == "current") {
          
          # Run for each variable
          for (sel_var in ds_vars) {
            
            # Get the final file name
            outfile <- paste0(outdir, "current/",
                              gsub("_mean", "", gsub("_2000_20[[:digit:]][[:digit:]]", "", tolower(id))),
                              "_", gsub("^(.*)_", "", sel_var), ".tif")
            
            # Check if exists
            check_file_exist(outfile, {
              
              # If it does not exists, try to download                        
              if (verbose) cli::cli_progress_step("Downloading {sel_var} - {period}", spinner = T,
                                msg_failed = "Variable {id} [{sel_var}] not available")
              
              var <- try(download_dataset(tolower(id), sel_var, list(time = time_steps[[ts]]), fmt = "raster",
                                          directory = rawdir,
                                          verbose = verbose), silent = T) # Set verbose=TRUE to debug
              
              # If download succeed, save
              if (!assertthat::is.error(var)) {
                terra::writeRaster(var, outfile, overwrite = T);rm(var)
              } else {
                if (verbose) cli::cli_progress_done(result = "failed")
              }
              # To ignore the file checking and download anyway, set this to TRUE
            }, ignore = !skip_exist)
            if (verbose) cli::cli_progress_done()
          }
          
        } else {
          
          # For the future we run for each scenario
          for (scen in future_scenarios) {
            # Get the modified dataset ID (each future have one ID)
            mod_id <- gsub("baseline_2000_20[[:digit:]][[:digit:]]",
                           paste0(scen, "_2020_2100"), id)
            
            # Run for each variable
            for (sel_var in ds_vars) {
              
              # Get the final file name
              outfile <- paste0(outdir, "future/", scen, "/",
                                gsub("_mean", "", gsub("_2020_2100", "", tolower(mod_id))),
                                "_", period, "_", gsub("^(.*)_", "", sel_var), ".tif")
              
              # Check if file exist
              check_file_exist(outfile, {
                
                # If it does not exists, try to download       
                if (verbose) cli::cli_progress_step("Downloading {scen} - {sel_var} - {period}", spinner = T,
                                  msg_failed = "Variable {mod_id}, scenario {scen} [{sel_var}], period {period} not available")
                
                var <- try(download_dataset(tolower(mod_id), sel_var, list(time = time_steps[[ts]]), fmt = "raster",
                                            directory = rawdir,
                                            verbose = verbose), silent = T) # Set verbose=TRUE to debug
                
                # If download succeed, save
                if (!assertthat::is.error(var)) {
                  terra::writeRaster(var, outfile, overwrite = T);rm(var)
                } else {
                  if (verbose) cli::cli_progress_done(result = "failed")
                }
                # To ignore the file checking and download anyway, set this to TRUE
              }, ignore = !skip_exist)
              if (verbose) cli::cli_progress_done()
            }
          }
        }
      }
      # If everything is done, conclude last message
      if (verbose) cli::cli_progress_done()
    }
  }
  
  if (!is.null(terrain_vars)) {
    for (tv in terrain_vars) {
      
      outfile <- paste0(outdir, "terrain/", tv, ".tif")
      
      check_file_exist(
        outfile,
        {
          # If it does not exists, try to download       
          if (verbose) cli::cli_progress_step("Downloading terrain {tv}", spinner = T,
                                              msg_failed = "Variable {tv} not available")
          
          var <- try(download_dataset("terrain_characteristics", tv,
                                      list(time = c("1970-01-01T00:00:00Z", "1970-01-01T00:00:00Z")), fmt = "raster",
                                      directory = rawdir,
                                      verbose = verbose), silent = T) # Set verbose=TRUE to debug
          
          # If download succeed, save
          if (!assertthat::is.error(var)) {
            terra::writeRaster(var, outfile, overwrite = T);rm(var)
          } else {
            if (verbose) cli::cli_progress_done(result = "failed")
          }
        }, ignore = !skip_exist)
    };if (verbose) cli::cli_progress_done()
  }
  
  
  # Delete raw files (optional, but recommended) ----
  if (!keep_raw) {
    if (verbose) cli::cli_alert_info("Deleting raw directory.")
    fs::dir_delete(rawdir)
  }
  
  return(invisible(NULL))
}
