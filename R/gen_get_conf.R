#' Generate configuration file in YAML format
#' 
#' This function will generate a configuration file in YAML format, which provides a 
#' convenient way to save and retrieve configurations for running the SDMs.
#' 
#' Configurations from this file can be retrieved using the function [get_conf()].
#' For details on how to fill the file, see comments on the generated file.
#'
#' @param filename a file name
#' @param outdir the output directory. If not existent, will be created
#'
#' @return saved file
#' @export
#'
#' @examples
#' \dontrun{
#' gen_configuration()
#' }
gen_configuration <- function(filename = "sdm_conf.yml", outdir = ".") {
  to_save <- 
"# Project title:
# Authors:
#
# This is a configuration file for the SDMs. It's used to configure (1) the groups
# for which to obtain the SDMs and (2) the variables to use when fitting the SDMs
# YAML files use the identation to establish list levels, so it's important to keep it consistent
# Don't use tabulation, use 2 spaces instead.
#
# File version: 1 (Y-M-D date)
---
# Fill here with the conditions to get the groups.
# For example, the group photosynthesizers have to met the condition
# kingdom == 'Chromista' | kingdom == 'Plantae'
# For character, use single ''
groups:
  photosynthesizers: kingdom == 'Chromista' | kingdom == 'Plantae'
# Put here the list of variables that will be used for each group
# The order is variables (don't edit) > group name (should be equal to the previous section) > 'hypothesis_name' > list of variables
# setting an hypothesis (by putting an 'hypothesis_name') is useful if you want to test multiple variables settings
# If you want to use a single hypothesis, simple put only one hypotehsis.
# Names of variable should be written on the format NAME_VARIANT(e.g. mean, min, max)
# So, for example, to get the mean temperature for the sea we would set: thetao_mean
# Example:
# variables:
#  fishes:
#    hypothesis1:
#      - variable1_mean
#      - variable2_max
#    hypothesis2:
#      - variable2_max
#      - variable3_min
variables:
  photosynthesizers:
    hypothesis1:
      - thetao_max
      - thetao_min
      - so_mean
"
  
  fs::dir_create(outdir)
  
  writeLines(to_save, con = paste0(outdir, "/", filename))
  
  cli::cli_alert_success("File {.file {paste0(outdir, '/', filename)}} saved")
  
  return(invisible(NULL))
}




#' Retrieve configuration settings from a configuration file
#'
#' @param conf_file the path for the configuration file
#' @param what what to retrieve (character vector)
#' 
#' This enable to retrieve configuration settings from a configuration YAML file
#' generated with [gen_configuration()]. Although this function can be used, we recommend
#' using directly the following two functions:
#' 
#' - [get_listbygroup()]: from the filters established in the "groups" section of
#' the configuration file, convert a list of species into a list with groups.
#' - [get_envofgroup()]: for a certain group, read and load the environmental files
#' saved in a folder, following the list of chosen variables.
#'
#' @return a list of configurations
#' @export
#'
#' @examples
#' \dontrun{
#' configurations <- get_conf()
#' }
get_conf <- function(conf_file = "sdm_conf.yml", what = c("groups", "variables")) {
  
  if (!file.exists(conf_file)) {
    cli::cli_abort("File {.file {conf_file}} does not exist.")
  }
  
  conf_res <- yaml::yaml.load_file(conf_file, readLines.warn = FALSE)
  
  if (!all(what %in% names(conf_res))) {
    cli::cli_abort("'what' should be one of {names(conf_res)}")
  }
  
  if (length(what) != length(names(conf_res))) {
    conf_res <- conf_res[what]
  }
  
  return(conf_res)
}



#' Includes a new column of SDM groups to a species list
#' 
#' This function will include a new grouping column to a species list (\code{data.frame} like objects)
#' based on the filters established in the configuration file.
#'
#' @param species_list the \code{data.frame} with the species list. Should have at least the
#' columns that were included as filters on the configuration file
#' @param conf_file the path for the configuration file
#' @param not_found_value which value use when the function is unable to establish a group
#' 
#'
#' @return the species list with a new column called "sdm_group"
#' @export
#'
#' @examples
#' \dontrun{
#' new_list <- get_listbygroup(sp_list)
#' }
get_listbygroup <- function(species_list, conf_file = "sdm_conf.yml",
                            not_found_value = "NOT_FOUND") {
  
  groups_conf <- get_conf(conf_file = conf_file, what = "groups")
  
  new_list <- species_list
  new_list$sdm_group <- NA
  new_list[is.na(new_list)] <- ""
  
  for (i in 1:length(groups_conf$groups)) {
    
    if (length(groups_conf$groups[[i]]) > 1) {
      cli::cli_abort("Something is wrong. Check if the condition for group {names(groups_conf$groups)[i]} is all in the same line.
                     >>> Example: {names(groups_conf$groups)[i]} == 'Value' & {names(groups_conf$groups)[i]} != 'Value <<<
                     Configuration retrieved from file {.file {conf_file}}")
    }
    
    conditions <- stringr::str_split(groups_conf$groups[[i]], "&\\|")
    # check if col exists
    to_check <- unlist(lapply(conditions[[1]], function(x){
      cond <- stringr::str_split(x, "==|!=|<=|>=|<|>|%in%")
      cond <- lapply(cond, function(z){
        z <- trimws(z[1])
        z <- gsub("!", "", z)
        z
      })
      unlist(cond)
    }))
    
    if (any(!to_check %in% colnames(species_list))) {
      cli::cli_abort("Condition columns not found in the species list. Available columns are {.emph {colnames(species_list)}}.
                     Check your configuration file: {.file {conf_file}}")
    }
    
    t_conditions <- groups_conf$groups[[i]]
    
    new_list <- eval(parse(text = glue::glue(
      "dplyr::mutate(new_list,
                     sdm_group = ifelse({t_conditions},
                                        names(groups_conf$groups)[i], sdm_group))"
    )))
    
  }
  
  new_list$sdm_group[is.na(new_list$sdm_group)] <- not_found_value
  
  return(new_list)
  
}



#' Load SpatRaster environmental layers for a group defined on a configuration
#' file
#'
#' @param group the group for which you want to retrieve the layers
#' @param depth the depth of the layer (usually one of "surf", "mean", "min" or
#'   "max")
#' @param scenario the scenario for the layers. Current or one of the
#'   future/past. See details
#' @param period optional character indicating the period
#' @param hypothesis if \code{NULL}, it will retrieve the first available
#'   hypothesis list. Otherwise, it will look for the list equivalent to the
#'   hypothesis. Should be a character (e.g. "hypothesis1")
#' @param conf_file the path for the YAML configuration file
#' @param env_folder the path for the environmental layers folder
#' @param fixed_name the name of the folder containing the fixed or non variable
#'   layers (usually terrain ones). Can also be NULL if you don't have fixed layers.
#' @param keep_in_future if \code{TRUE}, in case a future layer is not
#'   available, it will get its current period version
#' @param surface_ifnot if \code{TRUE}, it will use the surface layer in case it
#'   is not found for the required depth. If you named your surface layers in a way different than "depthsurf",
#'   then you can instead supply here the name of the surface depth.
#' @param accepted_formats a character vector of accepted formats. Those should be supported by [terra::rast()]
#' @param verbose print messages
#'
#' @return a SpatRaster with as many layers as present in the configuration file
#' @export
#'
#' @details Although you can indicate any path for the env folder, this folder
#' must be organized internally in the following way:
#' - current (for the primary fitting scenario it SHOULD be named current)
#' - \code{future_name} (the name of the folder for future layers)
#'  - {scenario1} (folder containing layers for the scenario1)
#'  - {scenario2} (folder containing layers for the scenario2)
#' - \code{fixed_name} (the name of the folder for fixed layers, i.e., layers that are non variable across periods like bathymetry)
#'
#' Note that you can use any folder name for your scenarios (e.g. ssp1, ssp2).
#' The function is not case sensitive, and will convert everything to lower
#' case.
#'
#' The function will automatically detect if the layer is fixed or not, by
#' searching in the fixed folder.
#' 
#' # Using user supplied layers vs package layers
#'
#' If you used the function [get_envlayers()] to download from the Bio-ORACLE,
#' then the function will work without any problem. However, if you added
#' additional layers take into account that the function expect that layers names contains three pieces of information 
#' separated by underscores or dashes: name of variable, depth,
#' and variant (e.g. max). One example: siconc_baseline_depthsurf_mean. Note that here you have 
#' additional information, what is fine if you have the essential information. For fixed ones,
#' you can have only the name (e.g. bathymetry).
#' 
#' Because you can have "depthmean" and a variant "mean", it is ESSENTIAL that the variant comes preceded by _ or -. For example, all those should be valid:
#' 
#' - siconc_baseline_depthmean_mean
#' - siconc_mean_baseline_depthmean
#' - siconc-baseline-mean-depthmean
#' 
#' But this one not work:
#' 
#' - mean_siconc_baseline_depthmean
#' 
#' Layers can have an additional information regarding the period. So, for example, future layers for the decade of 2050 can be named:
#' siconc_ssp126_depthsurf_dec50_mean
#' 
#' The function will do strict string match for everything, so ensure that you wrote the names correctly in your configuration file.
#' 
#' # Note
#' 
#' Layer names are changed to the retrieved combination of name + variation (e.g. thetao_mean). The old names are stored as an attribute
#' called "old_names" and can be accessed using \code{attr(raster_layers_object, "old_names")} for verification.
#' 
#' @seealso [get_listbygroup()], [get_conf()]
#'
#' @examples
#' \dontrun{
#' raster_layers <- get_envofgroup(group = "fishes")
#' }
get_envofgroup <- function(group,
                           depth = "depthsurf",
                           scenario = "current",
                           period = NULL,
                           hypothesis = NULL,
                           conf_file = "sdm_conf.yml",
                           env_folder = "data/env",
                           future_name = "future",
                           fixed_name = "terrain",
                           keep_in_future = TRUE,
                           surface_ifnot = TRUE,
                           accepted_formats = c("tif", "nc"),
                           verbose = TRUE) {
  
  if (!is.logical(surface_ifnot)) {
    surface_name <- surface_ifnot
    surface_ifnot <- TRUE
  } else {
    surface_name <- "depthsurf"
  }
  
  layers_list <- get_conf(conf_file = conf_file, what = "variables")
  
  layers_list <- layers_list$variables[[group]]
  
  if (is.null(hypothesis)) {
    layers_list <- layers_list[[1]]
  } else {
    layers_list <- layers_list[[hypothesis]]
  }
  
  layers_split <- stringr::str_split(layers_list, "_|-")
  
  layer_names <- unlist(lapply(layers_split, function(x){x[1]}))
  layer_variant <- unlist(lapply(layers_split, function(x){x[2]}))
  
  all_layers <- data.frame(name = layer_names,
                           variant = layer_variant)
  
  if (scenario == "current") {
    f_base <- list.files(paste0(env_folder, "/", scenario), full.names = T)
  } else {
    f_base <- list.files(paste0(env_folder, "/", future_name, "/", scenario), full.names = T)
  }
  
  f_base <- f_base[grepl(paste0(accepted_formats, collapse = "|"), f_base)]
  f_base <- f_base[!grepl("aux|xml", f_base)] # Remove aux or xml files
  
  f_terrain <- list.files(paste0(env_folder, "/", fixed_name), full.names = T)
  f_terrain <- f_terrain[grepl(paste0(accepted_formats, collapse = "|"), f_terrain)]
  f_terrain <- f_terrain[!grepl("aux|xml", f_terrain)] # Remove aux or xml files
  
  all_layers$file_path <- NA
  
  for (i in 1:nrow(all_layers)) {
    f_sel <- f_base[grepl(paste0(all_layers$name[i], "_"), f_base)]
    f_sel <- f_sel[grepl(paste0(c("_", "-"), all_layers$variant[i], collapse = "|"), f_sel)]
    f_final <- f_sel[grepl(depth, f_sel)]
    
    if (!is.null(period)) {
      f_final <- f_final[grepl(period, f_final)]
    }
    
    # If not found...
    if (length(f_final) < 1) {
      
      # Try surface
      if (surface_ifnot) {
        f_final <- f_sel[grepl(surface_name, f_sel)]
        
        if (!is.null(period)) {
          f_final <- f_final[grepl(period, f_final)]
        }
      }
      
      # If still not found...
      if (length(f_final) < 1) {
        
        # Try terrain
        f_sel <- f_terrain[grepl(all_layers$name[i], f_terrain)]
        if (!is.na(all_layers$variant[i])) {
          f_sel <- f_sel[grepl(paste0(c("_", "-"), all_layers$variant[i], collapse = "|"), f_sel)]
        }
        if (any(grepl("depth", f_sel))) {
          f_final <- f_sel[grepl(depth, f_sel)]
        } else {
          f_final <- f_sel
        }
        
        # If still not found...
        if (length(f_final) < 1) {
          
          # Try the current version (only non current scenarios)
          if (keep_in_future & scenario != "current") {
            
            f_base_kf <- list.files(paste0(env_folder, "/current"), full.names = T)
            f_base_kf <- f_base_kf[grepl(paste0(accepted_formats, collapse = "|"), f_base_kf)]
            f_base_kf <- f_base_kf[!grepl("aux|xml", f_base_kf)] # Remove aux or xml files
            
            f_sel <- f_base_kf[grepl(paste0(all_layers$name[i], "_"), f_base_kf)]
            f_sel <- f_sel[grepl(paste0(c("_", "-"), all_layers$variant[i], collapse = "|"), f_sel)]
            f_final <- f_sel[grepl(depth, f_sel)]
            
            if (length(f_final) < 1) {
              # Try surface
              if (surface_ifnot) {
                f_final <- f_sel[grepl(surface_name, f_sel)]
              }
            }
            
          }
          
        }
      }
      
    }
    
    if (length(f_final) > 1) {
      cli::cli_warn("More than 1 match found for {all_layers$name[i]}; the first will be used. Check your configuration file or layers! Matches found: {.emph {f_final}}. Match used: {.emph {f_final[1]}}.")
      f_final <- f_final[1]
    }
    
    if (length(f_final) < 1) {
      all_layers$file_path[i] <- NA;rm(f_final)
    } else {
      all_layers$file_path[i] <- f_final;rm(f_final)
    }
    
  }
  
  if (verbose & any(is.na(all_layers$file_path))) {
    cli::cli_alert_danger("The following variable{?s} {?was/were} not found: 
                          {paste(all_layers$name[is.na(all_layers$file_path)],
                          all_layers$variant[is.na(all_layers$file_path)],
                          sep = '_')}")
  }
  
  all_layers <- all_layers[!is.na(all_layers$file_path),]
  
  rast_loaded <- terra::rast(all_layers$file_path)
  
  attr(rast_loaded, "old_names") <- names(rast_loaded)
  
  names(rast_loaded) <- paste0(all_layers$name, ifelse(
    is.na(all_layers$variant),
    "",
    paste0("_", all_layers$variant)
  ))
  
  if (verbose) cli::cli_alert_info("Names were changed (from > to) {paste(attr(rast_loaded, 'old_names'), names(rast_loaded), sep = ' > ')}")
  
  return(rast_loaded)
}