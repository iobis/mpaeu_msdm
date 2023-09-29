############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Methods testing ################################

#' Generate virtual species for SDM testing
#'
#' Thin function is a wrapper around the [virtualspecies] package functions to
#' generate virtual species. The purpose of generating virtual species is to
#' test the capacity of distinct SDM algorithms to correctly return the "true"
#' distribution of species.
#'
#' @param layers environmental layers in raster* format (not terra compatible
#'   yet).
#' @param sim_funs functional relationships for the virtual species, generated
#'   with [virtualspecies::formatFunctions()]. All layers should have a
#'   corresponding function. If not supplied, then a PCA mode is used (then
#'   sim_pca should be supplied).
#' @param sim_pca the type of niche simulated by the PCA model. Can be "narrow"
#'   or "wide" (see details on [virtualspecies::generateSpFromPCA()]). Ignored
#'   if \code{sim_funs} is supplied.
#' @param samp_n_low  number of points for the low sampling scenario.
#' @param samp_n_high number of points for the highly sampled scenario.
#' @param samp_rep number of repetitions of the sampling.
#' @param samp_bias_layer an optional sample bias layer. Should be a raster with
#'   0-1 values.
#' @param samp_bias_only in case that a sampling bias layer is supplied, define
#'   if the sampling should be only with the bias or also without it. Default is
#'   FALSE, meaning that sampling will occur with and without bias.
#' @param samp_constr_shape an optional shapefile to constrain the sampling to
#'   that area.
#' @param samp_pa should absences also be returned? Interesting for comparisons
#'   with presence-absence methods. When this option is set to \code{TRUE},
#'   additional files with presence-absence are saved. Absences are sampled in
#'   the same number as presences (equal prevalence).
#' @param samp_replace if \code{TRUE}, then the sampling of species occurs with
#'   replacement and multiple points can be sampled on the same cell.
#' @param vsp_name a name for the virtual species. If \code{NULL} "Virtual
#'   species" will be used.
#' @param vsp_class a 'class' for the virtual species. Helpful to save
#'   information about the generating process. Can be for example
#'   "vsp_pca_narrow". If \code{NULL}, "standard_vsp" is used.
#' @param pred_fut an optional named list of raster* with future environmental
#'   conditions to predict the distribution of the species in future climates.
#'   The names of the list will be used to save the future suitability rasters.
#' @param save_key the key to save the species
#' @param save_path where to save the files. If \code{NULL} the standard format
#'   "data/virtual_species" will be used. Files will be saved according to the
#'   key provided (see details).
#' @param plot plot distributions and graphics as the virtual species is
#'   generated.
#' @param test_mode if set to \code{TRUE}, the function will stop before the
#'   sampling of points. This is mainly intended for testing if the
#'   environmental relationship functions are returning the expected suitability
#'   pattern (in this case, plot = \code{TRUE}).
#'
#' @details This function is a helper to generate virtual species for testing
#' SDM methods. Files will be saved in the path
#' data/virtual_species/key='save_key'/'file name' (note that the standard
#' folder is virtual_species, not species). A few points should be observed:
#' - the \code{save_key} should be distinct from the AphiaIDs of the project species. You can use any number, but it's important to be consistent to avoid any confusion.
#' - different than the others, the virtual species folder will contain more than only the species occurrence table.
#' - the species occurrence information will be saved in the same way as the other species, but will contain only 6 collumns: key (save_key), scientificName (which will receive the vsp_name),
#' decimalLongitude, decimalLatitude, presence (1-0, denoting presence or
#' absence) and class, which will hold the vsp_class.
#'
#' The files that will be saved in the folder are the following:
#' - parquet files containing the virtual species occurrence for each sampling schema (very low, low and high)
#' - if samp_pa is \code{TRUE}, the same files containing the absence info (saved as 'file*_pa.parquet')
#' - raster files containing the suitability, presence-absence, probabilities of occurrence
#' - if pred_fut is supplied, rasters containing the future distribution
#'
#' NOTE: the process is random and you may want to set a seed before running the
#' function.
#'
#' @return saved files in the project format
#' @export
#' 
#' @examples
#' \dontrun{
#' gen_vsp()
#' }
gen_vsp <- function(layers,
                    sim_funs = NULL, sim_pca = NULL,
                    samp_n_low = 30, samp_n_high = 150, samp_rep = 10,
                    samp_bias_layer = NULL, samp_bias_only = FALSE, samp_constr_shape = NULL, samp_pa = TRUE, samp_replace = FALSE,
                    vsp_name = NULL, vsp_class = NULL,
                    pred_fut = NULL, save_path = NULL, save_key = NULL, plot = TRUE,
                    test_mode = FALSE) {
  
  sm <- function(x){suppressMessages(x)}
  
  if (is.null(vsp_name)) {
    vsp_name <- "Virtual species"
  }
  if (is.null(vsp_class)) {
    vsp_class <- "standard_vsp"
  }
  if (is.null(sim_funs) & is.null(sim_pca)) {
    stop("If 'sim_funs' is not provided, then 'sim_pca' should be provided.")
  }
  
  if (!is.null(pred_fut)) {
    if (!is.list(pred_fut)) {
      stop("pred_fut should be a named list of RasterStack.")
    } else {
      if (is.null(names(pred_fut))) {
        stop("pred_fut should be named, as the names are used to save the files.
Example: pred_fut <- list(ssp1 = ssp1_stack, ssp2 = ssp2_stack)")
      }
    }
  }
  
  if (is.null(sim_funs)) {
    vsp_mode <- "pca"
  } else {
    vsp_mode <- "fun"
  }
  
  if (test_mode) {
    plot <- TRUE
  }
  
  if (is.null(save_path)) {
    save_path <- "data/virtual_species"
  }
  if (is.null(save_key)) {
    cli::cli_alert_danger("save_key not supplied, using 1 as key")
    save_key <- 1
  }
  
  save_path <- paste0(save_path, "/key=", save_key, "/")
  fs::dir_create(save_path)
  
  cli::cli_progress_step("Generating species name={.emph {vsp_name}} | key={save_key} | class={vsp_class}")
  
  # Generate species
  if (vsp_mode == "fun") {
    vsp_sp <- sm(virtualspecies::generateSpFromFun(layers,
                                                   parameters = sim_funs,
                                                   plot = plot))
  } else {
    vsp_sp <- sm(virtualspecies::generateSpFromPCA(raster.stack = layers, 
                                                   sample.points = TRUE,
                                                   niche.breadth = sim_pca,
                                                   plot = plot))
  }
  
  cli::cli_progress_step("Converting to presence-absence raster")
  # Convert to presence absence
  vsp_sp_pa <- sm(virtualspecies::convertToPA(vsp_sp,
                                              PA.method = "probability",
                                              prob.method = "logistic",
                                              beta = 0.5, alpha = -0.05,
                                              plot = plot))
  
  if (test_mode) {
    cli::cli_progress_done()
    cli::cli_alert_info("Early stopping (test mode).")
    return(invisible(NULL))
  }
  
  nsamp <- c(samp_n_low, samp_n_high)
  names(nsamp) <- c("low", "high")
  bias_type <- rep("no.bias", 2)
  samp_bias_layer_l <- NULL
  
  if (!is.null(samp_bias_layer)) {
    if (!samp_bias_only) {
      nsamp <- c(nsamp, nsamp)
      names(nsamp)[3:4] <- paste0(names(nsamp)[3:4], "_bias")
      samp_bias_layer_l <- list(NULL, NULL,
                                samp_bias_layer, samp_bias_layer)
      bias_type <- c(rep("no.bias", 2), rep("manual", 2))
    } else {
      names(nsamp) <- paste0(names(nsamp), "_bias")
      samp_bias_layer_l <- list(samp_bias_layer, samp_bias_layer)
      bias_type <- rep("manual", 2)
    }
  }
  cli::cli_progress_done()
  cli::cli_alert_info("Sampling points ({samp_rep} sampling{?s})")
  cli::cli_progress_bar(total = samp_rep)
  for (i in 1:samp_rep) {
    for (z in 1:length(nsamp)) {
      vsp_occ <- sm(virtualspecies::sampleOccurrences(vsp_sp_pa,
                                                      n = nsamp[z],
                                                      type = "presence only",
                                                      sampling.area = samp_constr_shape,
                                                      bias = bias_type[z],
                                                      weights = samp_bias_layer_l[[z]],
                                                      replacement = samp_replace,
                                                      plot = plot))
      
      vsp_dat <- vsp_occ$sample.points[,1:3]
      colnames(vsp_dat) <- c("decimalLongitude", "decimalLatitude", "presence")
      
      vsp_dat$key <- save_key
      vsp_dat$scientificName <- vsp_name
      vsp_dat$class <- vsp_class
      
      arrow::write_parquet(vsp_dat, paste0(save_path, "occurrences_", names(nsamp)[z], "_rep", i, ".parquet"))
      
      if (samp_pa) {
        vsp_occ_pa <- sm(virtualspecies::sampleOccurrences(vsp_sp_pa,
                                                           n = nsamp[z]*2,
                                                           type = "presence-absence",
                                                           sampling.area = samp_constr_shape,
                                                           bias = bias_type[z],
                                                           weights = samp_bias_layer_l[[z]],
                                                           replacement = samp_replace,
                                                           sample.prevalence = 0.5,
                                                           plot = plot))
        
        vsp_dat_pa <- vsp_occ_pa$sample.points[,1:3]
        colnames(vsp_dat_pa) <- c("decimalLongitude", "decimalLatitude", "presence")
        
        vsp_dat_pa$key <- save_key
        vsp_dat_pa$scientificName <- vsp_name
        vsp_dat_pa$class <- vsp_class
        
        arrow::write_parquet(vsp_dat_pa, paste0(save_path, "occurrences_", names(nsamp)[z], "_rep", i, "_pa.parquet"))
        
      }
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  cli::cli_progress_step("Saving rasters")
  
  terra::writeRaster(vsp_sp$suitab.raster, paste0(save_path, "suitability_key", save_key,".tif"),
              overwrite = T)
  
  terra::writeRaster(vsp_sp_pa$probability.of.occurrence, paste0(save_path, "probability_key", save_key,".tif"),
              overwrite = T)
  
  terra::writeRaster(vsp_sp_pa$pa.raster, paste0(save_path, "pa_key", save_key,".tif"),
              overwrite = T)
  
  saveRDS(vsp_sp, file = paste0(save_path, "vsp_key", save_key,".rds"))
  
  if (!is.null(pred_fut)) {
    cli::cli_progress_step("Making future predictions")
    for (z in 1:length(pred_fut)) {
      
      scenario <- names(pred_fut)[z]
      
      if (vsp_mode == "fun") {
        vsp_sp_fut <- sm(virtualspecies::generateSpFromFun(pred_fut[[z]],
                                                           parameters = sim_funs,
                                                           plot = plot))
      } else {
        vsp_sp_fut <- sm(virtualspecies::generateSpFromPCA(raster.stack = pred_fut[[z]],
                                                           pca = vsp_sp$details$pca,
                                                           means = vsp_sp$details$means,
                                                           sds = vsp_sp$details$sds,
                                                           plot = plot))
      }
      
      terra::writeRaster(vsp_sp_fut$suitab.raster, paste0(save_path, scenario, "_suitability_key", save_key, ".tif"),
                         overwrite = T)
    }
  }
  
  cli::cli_progress_done()
  cli::cli_alert_success("Species {vsp_name}:{save_key}:{vsp_class} concluded.")
  
  return(invisible(NULL))
  
}
