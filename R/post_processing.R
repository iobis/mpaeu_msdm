# post_list_s3 <- function(
#     bucket = "mpaeu-dist", folder = "results", outfile = "s3_list.parquet"
# ) {
#     cli::cli_alert_info("Retrieving S3 list, this may take a while...")

#     require(dplyr)

#     bucket_list <- aws.s3::get_bucket_df(
#         bucket = bucket,
#         prefix = folder,
#         use_https = TRUE,
#         max = Inf
#     )

#     cli::cli_alert_info("Processing...")

#     bucket_list <- bucket_list[grepl("results", bucket_list$Key),]

#     category <- sub("results/*/", "\\1", bucket_list$Key)
#     category <- sub("*/.*", "\\1", category)

#     bucket_list$category <- category

#     taxonid <- unlist(
#         lapply(1:nrow(bucket_list), function(x) {
#             m <- regmatches(bucket_list$Key[x], regexpr("(?<=taxonid=)\\d+", bucket_list$Key[x], perl = TRUE))
#             if (length(m) > 0) {
#                 m
#             } else {
#                 NA
#             }
#         })
#     )

#     bucket_list$taxonID <- as.integer(taxonid)

#     model_acro <- regmatches(bucket_list$Key, regexpr("model=[^/]+", bucket_list$Key))
#     model_acro <- sub("^(model=[^_]+(?:_[^_=]+)*)_.*", "\\1", model_acro)
#     model_acro <- gsub("model=", "", model_acro)
#     bucket_list$model_acro <- model_acro

#     method <- unlist(lapply(1:nrow(bucket_list), function(x) {
#         m <- regmatches(bucket_list$Key[x], regexpr("method=[^/]+", bucket_list$Key[x]))
#         if (length(m) > 0) {
#             m
#         } else {
#             NA
#         }
#     }))
#     method <- sub("^(method=[^_]+(?:_[^_=]+)*)_.*", "\\1", method)
#     method <- gsub("method=", "", method)

#     bucket_list$models <- method

#     bucket_list <- bucket_list %>%
#         group_by(taxonID, models) %>%
#         mutate(is_boot = ifelse(grepl("bootcv", Key), TRUE, FALSE))

#     arrow::write_parquet(bucket_list, outfile)

#     cli::cli_alert_success("File saved at {.file {outfile}}")

#     return(invisible(NULL))
# }


#' Prepare files for use in Zonation
#'
#' @param species a vector of species IDs (AphiaID)
#' @param source_folder the path to the results folder, typically "results"
#' @param study_area the path to the study area shapefile (any format readable by [terra::vect()])
#' @param model_acro the modelling realization acronym (e.g. "mpaeu")
#' @param target_mask which mask layer/type to use. Should be available on the mask file
#' @param target_threshold which threshold to use to remove areas with lower probability of occurrence
#' @param target_model which model (algorithm) to be used. Can be a single value (e.g. "maxent") or a vector. 
#' In case a vector is supplied, the function will check which of those are available (in the order supplied)
#' and then select the first one available. So, if you supply c("maxent", "rf", "ensemble") and "maxent" is not available,
#' but the others are, the function will use Random Forest ("rf")
#' @param future_scenarios which scenarios you want to also process (in the format "sspXXX"). Current scenario is always included.
#' @param future_periods which periods (in the format "decXX") to use. Current period is always included.
#' @param boot_weight weight to be applied to the bootstrap layer
#' @param outfolder folder to save the processed layers. If not existing, it will be created
#' @param verbose if `TRUE` print messages
#'
#' @details
#' If you followed the project structure, you should not have any problem with this function as all files will be available for use.
#' However, the function includes several flow control steps. It will hardly abort (except if no data is available at all),
#' but will instead skip species for which one problem exists. It will record the problem in a `data.frame` 
#' which is returned by the function. The possible problems are:
#' - no files available: folder for the species exist, but no prediction is available  
#' - none of the specified models available: no model matching the argument `target_model`  
#' - no bootstrap file available: predictions available, but no bootstrap file is availabe  
#' - no bootstrap for the preferred model - selecting other: if there are bootstrap, but not for the target model it will chose the next available  
#' - bootstrap not available for all predictions - skipping: bootstrap files exist but not for all scenarios  
#'
#' @return processed files and a data.frame with the processing status
#' @export
#'
#' @examples
#' \dontrun{
#' post_prepare(species = c(124287, 137098), 
#'              source_folder = "results", 
#'              study_area = "data/shapefiles/mpaeu_studyarea_v2.shp", 
#'              model_acro = "mpaeu")
#' }
post_prepare <- function(
    species,
    source_folder,
    study_area,
    model_acro = "mpaeu",
    target_mask = "fit_region",
    target_threshold = "p10",
    target_model = c("maxent", "rf", "xgboost", "ensemble", "esm"),
    future_scenarios = c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
    future_periods = c("dec50", "dec100"),
    boot_weight = 0.5,
    outfolder = "proc-layers",
    #s3_list = NULL,
    verbose = TRUE
) {

    require(dplyr)
    fs::dir_create(outfolder)
    study_area <- terra::vect(study_area)

    # Already check which folders are available
    if (verbose) cli::cli_alert_info("Checking which species are available.")
    source_folder_content <- list.files(
        source_folder, full.names = T
    )
    source_folder_base <- basename(source_folder_content)
    source_folder_base <- gsub("taxonid=", "", source_folder_base[grepl("taxonid=", source_folder_base)])
    source_folder_base <- as.numeric(source_folder_base)

    species <- species[species %in% source_folder_base]

    if (length(species) < 1) {
        cli::cli_abort("No species available. Check your source folder.")
    } else {
        rm(source_folder_content, source_folder_base)
    }

    if (verbose) cli::cli_alert_info("Checking which species are available for this model acronym.")
    source_folder_content <- list.files(
        file.path(source_folder, paste0("taxonid=", species))
    )
    species <- species[which(source_folder_content == paste0("model=", model_acro))]

    if (length(species) < 1) cli::cli_abort("No species available for that model acronym.")

    # if (!is.null(s3_list)) {
    #     s3_list <- arrow::open_dataset(s3_list)
    #     s3_av <- TRUE
    # } else {
    #     s3_av <- FALSE
    # }

    total <- length(species)
    control <- rep(NA, total)
    for (sp in seq_len(total)) {

        if (verbose) cli::cli_alert_info("Processing {sp} of {total}.")

        id <- species[sp]

        # if (s3_av) {
        #     species_dat <- s3_list |> filter(taxonID == id) |> collect() |> pull(Key)
        #     files_available <- file.path(source_folder, species_dat)
        # } else {
        #     species_dat <- list.files(file.path(source_folder, paste0("taxonid=", id)), full.names = T)
        #     species_dat <- species_dat[grepl(model_acro, species_dat)]
        #     files_available <- list.files(species_dat, full.names = T, recursive = T)
        # }
        species_dat <- list.files(file.path(source_folder, paste0("taxonid=", id)), full.names = T)
        species_dat <- species_dat[grepl(model_acro, species_dat)]
        files_available <- list.files(species_dat, full.names = T, recursive = T)

        if (length(files_available) < 1) {
            control[sp] <- "no files available"
            next
        }

        thresholds <- files_available[grepl("what=thresholds", files_available)]
        masks <- files_available[grepl("mask", files_available)]

        is_available <- lapply(target_model, \(y) any(grepl(y, files_available)))
        is_available <- unlist(is_available, use.names = F)

        if (!any(is_available)) {
            control[sp] <- "none of the specified models available"
            next
        }

        boot_files <- files_available[grepl("what=boot", files_available)]

        is_available_boot <- lapply(target_model, \(y) any(grepl(y, boot_files)))
        is_available_boot <- unlist(is_available_boot, use.names = F)

        if (!any(is_available_boot)) {
            control[sp] <- "no bootstrap file available"
            next
        } else if (target_model[is_available][1] != target_model[is_available_boot][1]) {
            control[sp] <- "no bootstrap for the preferred model - selecting other"
            best_model <- target_model[is_available_boot][1]
        } else {
            best_model <- target_model[is_available][1]
        }

        model_threshold <- arrow::read_parquet(thresholds) |> 
            filter(model == best_model) |> 
            pull(target_threshold)

        model_files <- files_available[grepl(best_model, files_available)]

        if (length(future_scenarios) > 0) {
            t_future_scenarios <- paste0(future_scenarios, "_", rep(future_periods, each = length(future_scenarios)))
            model_files <- model_files[grepl(paste(c("current", t_future_scenarios), collapse = "|"), model_files)]
        } else {
            model_files <- model_files[grepl("current", model_files)]
        }

        # Check if number of bootstrap and model files are the same
        if (length(model_files[!grepl("what=boot", model_files)]) != length(model_files[grepl("what=boot", model_files)])) {
            control[sp] <- "bootstrap not available for all predictions - skipping"
            next
        }

        model_files <- model_files[!grepl("what=boot", model_files)]

        mask_l <- terra::rast(masks)
        mask_l <- terra::subset(mask_l, target_mask)
        terra::NAflag(mask_l) <- 0

        mask_l <- terra::mask(mask_l, study_area)
        mask_l <- terra::crop(mask_l, study_area)

        if (verbose) cli::cli_progress_bar("Processing files", total = length(model_files))
        for (mf in model_files) {
            processed_rast <- terra::rast(mf) |>
                terra::crop(y = mask_l) |>
                terra::mask(mask = mask_l) |>
                terra::classify(rcl = cbind(-Inf, (model_threshold * 100), 0)) |>
                terra::classify(rcl = cbind(NA, 0)) |>
                terra::mask(mask = study_area)

            processed_boot <- gsub("_cog", "_what=bootcv_cog", mf) |>
                terra::rast() |>
                terra::crop(y = mask_l) |>
                terra::mask(mask = mask_l) |>
                terra::clamp(upper=100) |> # temporary
                terra::classify(rcl = cbind(NA, 0)) |>
                terra::mask(mask = study_area)

            boot_weighted <- processed_boot * boot_weight

            #create combined layer of uncertainty discount (=prediction-weighted bootstrap)
            proc_weighted <- processed_rast - boot_weighted 
            #if final negative value, replace it with zero
            proc_weighted <- terra::classify(proc_weighted, cbind(-Inf, 0, 0))
            proc_weighted <- terra::as.int(proc_weighted)

            terra::writeRaster(proc_weighted,
                        file.path(outfolder, gsub("cog", "proc", basename(mf))),
                        overwrite = TRUE)

            if (verbose) cli::cli_progress_update()
        }
        if (verbose) cli::cli_progress_done()
        control[sp] <- "processed"
    }

    return(data.frame(species = species, status = control))
}
