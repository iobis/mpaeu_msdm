pred_assets <- c(current = "projects/ee-silascprincipe/assets/env_data")
occurrence_asset <- "projects/ee-silascprincipe/assets/vsp_points_test"



sdm_geemod_maxent <- function(sdm_data, method = "maxnet", 
                              tune_blocks = "spatial_grid") {
  
  .gee_checks(sdm_data)
  
  timings <- .get_time()
  
  
  # Retrieve basic info ----
  nback_number <- 150000
  k_folds <- 5
  
  
  # Get species name(s) ----
  sp_exec <- ee$List(sdm_data$species)
  
  # Load species data ----
  if (mode == "local") {
    
    # TO DO!
    
  } else {
    occurrence_data <- ee$FeatureCollection(occurrence_asset)
  }
  
  
  # Load environmental data ----
  region_rast <- ee$Image("projects/ee-silascprincipe/assets/env_data")
  for (i in 1:length(pred_assets)) {
    eval(parse(text = paste0(
      names(pred_assets)[i], ' <- ee$Image("', pred_assets[i],'")'
    )))
  }
  
  eval(parse(text = paste0(
    "pred_col <- ee$ImageCollection(list(", paste0(
      names(pred_assets[i]), collapse = ","
    ), "))" 
  )))
  
  
  # Prepare folds info ----
  n_iterations <- 100
  
  # TO DO: improve for masking
  # values(sdm_data$spatial_grid) <- 1
  # spat_grid_m <- terra::mask(sdm_data$spatial_grid, )
  spat_grid_m <- sdm_data$spatial_grid
  
  spatial_grid <- terra::as.polygons(spat_grid_m)
  
  spatial_grid <- sf::st_as_sf(spatial_grid)
  spatial_grid$cell <- 1:nrow(spatial_grid)
  
  folds <- lapply(1:n_iterations,
                  function(x){sample(1:k_folds, nrow(spatial_grid), replace = T)})
  
  # Create an ee version of spatial grid
  spatial_grid_ee <- sf_as_ee(spatial_grid)
  
  
  # Prepare tune parameters ----
  # Set parameters
  tune_params <- expand.grid(remult = seq(0.5, 2, 0.5), features = c("lq", "lqh"),
                             stringsAsFactors = F)
  
  tune_params <- unique(tune_params)
  
  tune_params$features <- ifelse(grepl("h", tune_params$features), TRUE, FALSE)
  
  tune_params <- lapply(1:nrow(tune_params), function(x){
    regu <- tune_params[x, 1]
    hing <- tune_params[x, 2]
    return(list(regu, hing))
  })
  
  # Convert to ee format
  tune_list <- ee$List(tune_params)
  
  
  
  # Execute function ----
  sp_results <- sp_exec$map()
  
  
  
}


# GEE maxent function
.gee_maxent <- function(lval){
  
  # Load species data ----
  sp_code <- lval
  
  sp_points <- occurrence_data$filter(
    ee$Filter$eq("key", ee$Number(sp_code))
  )
  
  
  
  # Get best spatial grid blocks ----
  folds_ee <- ee$List(folds)
  
  it_ind <- ee$List(0:(n_iterations - 1))
  
  cells_id <- spatial_grid_ee$aggregate_array("cell")
  
  find_best_block <- function(index) {
    
    # Create folds IDs
    folds_id <- ee$List(folds_ee$get(index))
    
    
    # Remap the grid
    spatial_grid_fold <- spatial_grid_ee$remap(
      cells_id, folds_id, "cell"
    )
    
    # Extract from shapefile
    # Create a filter
    spat_filter <- ee$Filter$intersects(
      leftField = ".geo",
      rightField = ".geo",
      maxError = 10
    )
    
    # Create a join
    spat_join <- ee$Join$saveAll(
      matchesKey = "matches"
    )
    
    # Join and get only information
    points_info <- spat_join$apply(sp_points, spatial_grid_fold, spat_filter)
    
    extract_info <- function(point) {
      point <- ee$Feature(point)
      matches <- ee$List(point$get("matches"))
      properties <- matches$map(ee_utils_pyfunc(function(shape){
        ee$Feature(shape)$get("cell")
      }))
      point$set("grid_cell", properties)$select(propertySelectors = c("presence", "grid_cell"))
    }
    points_with_grid <- points_info$map(extract_info)
    
    # Calculate distribution
    folds_dist <- points_with_grid$reduceColumns(
      selectors = list("grid_cell", "presence"),
      reducer = ee$Reducer$sum()$group(
        groupField = 0,
        groupName = "fold"
      )
    )
    
    # Extract only relevant field
    folds_metrics <- ee$List(folds_dist$get('groups'))
    
    folds_metrics_final <- folds_metrics$map(ee_utils_pyfunc(
      function(feat){
        ee$Dictionary(feat)$get(key = "sum")
      }
    ))
    
    new_result <- ee$List(
      list(ee$Number(index),
           ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
           ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero()))
    )
    
    new_result
  }
  
  all_blocks <- it_ind$map(ee_utils_pyfunc(find_best_block))
  
  all_blocks_col <- ee$FeatureCollection(all_blocks$map(ee_utils_pyfunc(
    function(val){
      ee$Feature(NULL, ee$Dictionary$fromLists(
        ee$List(list("combination", "eveness", "all_non_zero")),
        val
      ))
    }
  )))
  
  all_blocks_col <- all_blocks_col$
    filter("all_non_zero == 1")
  
  best_block <- 
    ee$Number(all_blocks_col$
                filter(ee$Filter$lte("eveness", all_blocks_col$aggregate_min("eveness")))$
                aggregate_array("combination")$
                get(0))
  
  # Get final remap
  spatial_grid_final <- spatial_grid_ee$remap(
    cells_id, ee$List(folds_ee$get(best_block)), "cell"
  )
  
  
  # Sample background points ----
  zones <- region_rast$select(1)$gt(0)
  zones <- zones$updateMask(zones$neq(0))
  
  region_pol <- zones$addBands(region_rast$select(1))$reduceToVectors(
    crs = region_rast$projection(),
    geometryType = 'polygon',
    reducer = ee$Reducer$mean()
  )
  
  back_points <- zones$sample(
    region = zones$geometry(),
    geometries = TRUE
  )
  
  back_points <- back_points$randomColumn(seed = 2)
  
  perc_to_get <- ee$Number((nback_number * 100))$divide(ee$Number(back_points$size()))$ceil()
  
  back_points <- back_points$filter(ee$Filter$lte("random", ee$Number(perc_to_get$divide(100))))
  
  back_points <- back_points$map(rgee::ee_utils_pyfunc(function(feat){
    point <- ee$Feature(feat)
    point <- point$set("presence", 0L)$select(propertySelectors = list("presence"))
  }))
  
  # ee_print(back_points_b)
  
  # Merge background and points ----
  all_points <- sp_points$select("presence")$merge(back_points)
  
  spat_filter <- ee$Filter$intersects(
    leftField = ".geo",
    rightField = ".geo",
    maxError = 10
  )
  
  # Create a join
  spat_join <- ee$Join$saveAll(
    matchesKey = "matches"
  )
  
  # Join and get only information
  points_info <- spat_join$apply(all_points, spatial_grid_final, spat_filter)
  
  extract_info <- function(point) {
    point <- ee$Feature(point)
    matches <- ee$List(point$get("matches"))
    properties <- matches$map(ee_utils_pyfunc(function(shape){
      ee$Feature(shape)$get("cell")
    }))
    point$set("fold", properties)$select(propertySelectors = c("presence", "fold"))
  }
  
  all_points_with_grid <- points_info$map(extract_info)
  
  
  
  # Extract information at the points ----
  training_dat <- region_rast$sampleRegions(collection = all_points_with_grid)
  
  
  # Get best tune ----
  
  # Creates a function to get AUC
  get_auc_tun <- function(featurecol) {
    
    predicted <- ee$FeatureCollection(featurecol)$aggregate_array("probability")
    response <- ee$FeatureCollection(featurecol)$aggregate_array("presence")
    
    thresholds <- ee$List$sequence(0, 1, 0.01)
    
    thresh_res <- function(th) {
      pred_class <- ee$Array(predicted)$gte(th)
      actu_class <- ee$Array(response)#$gte(th)
      ee$List(list(pred_class, actu_class))
    }
    
    thresholded <- thresholds$map(rgee::ee_utils_pyfunc(thresh_res))
    
    get_classind <- function(classified) {
      
      predvals <- ee$Array(ee$List(classified)$get(0))
      truevals <- ee$Array(ee$List(classified)$get(1))
      
      newvar <- ee$Array(truevals$multiply(2)$add(predvals))
      
      tp <- newvar$eq(3)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
      fp <- newvar$eq(1)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
      tn <- newvar$eq(0)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
      fn <- newvar$eq(2)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
      
      tpr <- tp$divide(tp$add(fn)) # TP/TP+FN / Recall
      tnr <- tn$divide(tn$add(fp)) # TN/TN+FP
      fpr <- fp$divide(fp$add(tn)) # FP/FP+TN
      pre <- tp$divide(tp$add(fp)) # TP/TP+FP
      
      ee$Dictionary(list(tp = tp, fp = fp, tn = tn, fn = fn, tpr = tpr, tnr = tnr,
                         fpr = fpr, pre = pre))
      
    }
    
    classif <- thresholded$map(rgee::ee_utils_pyfunc(get_classind))
    
    FPR <- classif$map(rgee::ee_utils_pyfunc(function(dict){
      value = ee$Dictionary(dict)$get('fpr')
      ee$Array(value)$get(ee$List(list(0)))
    }))
    
    TPR <- classif$map(rgee::ee_utils_pyfunc(function(dict){
      value = ee$Dictionary(dict)$get('tpr')
      ee$Array(value)$get(ee$List(list(0)))
    }))
    
    # Concatenate FPR and TPR arrays along axis 1
    X <- ee$Array(FPR) # tentar nao usar array, ver amanha
    Y <- ee$Array(TPR)
    X1 <- X$slice(0,1)$subtract(X$slice(0,0,-1))
    Y1 <- Y$slice(0,1)$add(Y$slice(0,0,-1))
    auc <- X1$multiply(Y1)$multiply(0.5)$reduce('sum', ee$List(list(0)))$abs()$toList()$get(0)
    
    #auc_feat <- ee$Feature(NULL, list(auc = auc))
    
    #return(ee$Number(auc))
    ee$Feature(NULL, list(auc = auc))
    
  }
  
  # Creates function to perform tuning
  tune_model <- function(listval) {
    
    # Get parameters
    all_parameters <- ee$List(listval)
    
    regmult <- all_parameters$get(0)
    hingefeat <- all_parameters$get(1)
    
    cv_model <- function(fold) {
      
      training = training_dat$filter(ee$String(ee$String("fold != ")$cat(ee$Number(fold)$format())))
      testing = training_dat$filter(ee$String(ee$String("fold == ")$cat(ee$Number(fold)$format())))
      
      maxclassifier <- ee$Classifier$amnhMaxent(
        autoFeature = FALSE,
        linear = TRUE,
        quadratic = TRUE,
        product = FALSE,
        threshold = FALSE,
        hinge = hingefeat, # Only hinge variable
        addSamplesToBackground = FALSE,
        betaMultiplier = ee$Number(regmult) # regularization multiplier
      )$train(
        features = training,
        classProperty = 'presence',
        inputProperties = region_rast$bandNames()
      )
      
      testpred <- ee$FeatureCollection(testing$classify(maxclassifier))
      
      auc_val <- get_auc_tun(testpred)
      
      return(auc_val)
      
    }
    
    fold_list <- ee$List$sequence(1, k_folds, 1)
    
    auc_cv <- fold_list$map(rgee::ee_utils_pyfunc(cv_model))
    
    mean_auc <- ee$FeatureCollection(auc_cv)$
      aggregate_mean("auc")
    
    return(ee$Feature(NULL, list(auc = mean_auc)))
  }
  
  tune_results <- tune_list$map(rgee::ee_utils_pyfunc(tune_model))
  
  tune_results_fc <- ee$FeatureCollection(tune_results)
  
  tr_list <- ee$List(tune_results_fc$aggregate_array("auc"))
  
  assign_list <- ee$List$sequence(0, 7, 1) # If more tuning parameters, upd here
  
  tune_res_fin <- assign_list$map(
    rgee::ee_utils_pyfunc(function(index){
      sel_list <- ee$List(tune_list$get(index))
      regmult <- sel_list$get(0)
      hingefeat <- sel_list$get(1)
      auc_val <- tr_list$get(index)
      return(ee$Feature(NULL, list(auc = auc_val,
                                   regmult = regmult,
                                   hingefeat = hingefeat)))
    })
  )
  
  best_tune <- ee$FeatureCollection(tune_res_fin)$
    filter(ee$Filter$eq("auc", ee$FeatureCollection(tune_res_fin)$
                          aggregate_min("auc")))
  
  
  # Cross validate best model ----
  regmult <- ee$List(best_tune$aggregate_array("regmult"))$get(0)
  hingefeat <- ee$List(best_tune$aggregate_array("hingefeat"))$get(0)
  
  cv_model_final <- function(fold) {
    
    training = training_dat$filter(ee$String(ee$String("fold != ")$cat(ee$Number(fold)$format())))
    testing = training_dat$filter(ee$String(ee$String("fold == ")$cat(ee$Number(fold)$format())))
    
    maxclassifier <- ee$Classifier$amnhMaxent(
      autoFeature = FALSE,
      linear = TRUE,
      quadratic = TRUE,
      product = FALSE,
      threshold = FALSE,
      hinge = hingefeat, # Only hinge variable
      addSamplesToBackground = FALSE,
      betaMultiplier = ee$Number(regmult) # regularization multiplier
    )$train(
      features = training,
      classProperty = 'presence',
      inputProperties = region_rast$bandNames()
    )
    
    testpred <- ee$FeatureCollection(testing$classify(maxclassifier))
    
    # testpred_sel <- testpred$select(
    #   propertySelectors = list("fold", "presence", "probability"),
    #   retainGeometry = FALSE
    # )
    testpred_sel <- testpred$select(
      list("fold", "presence", "probability")
    )
    
    return(ee$FeatureCollection(testpred_sel))
    
  }
  
  fold_list <- ee$List$sequence(1, k_folds, 1)
  
  final_m_cv <- fold_list$map(rgee::ee_utils_pyfunc(cv_model_final))
  
  final_m_cv <- ee$FeatureCollection(final_m_cv)$flatten()
  
  
  # Train full model ----
  maxclassifier_final <- ee$Classifier$amnhMaxent(
    autoFeature = FALSE,
    linear = TRUE,
    quadratic = TRUE,
    product = FALSE,
    threshold = FALSE,
    hinge = hingefeat, # Only hinge variable
    addSamplesToBackground = FALSE,
    betaMultiplier = ee$Number(regmult) # regularization multiplier
  )$train(
    features = training_dat,
    classProperty = 'presence',
    inputProperties = region_rast$bandNames()
  )
  
  # Predict to full data ----
  full_data_pred <- ee$FeatureCollection(training_dat$classify(maxclassifier_final))
  
  full_data_pred <- full_data_pred$select(
    list("fold", "presence", "probability")
  )
  
  
  # Predict to scenarios ----
  predictions <- pred_col$map(
    rgee::ee_utils_pyfunc(
      function(img){
        pred_img <- img$classify(maxclassifier_final)
        pred_img
      }
    )
  )
  
  # Export predictions ----
  for (i in 1:length(pred_assets)) {
    task_img <- ee_image_to_drive(
      image = ee$Image(predictions$get((i-1))),
      fileFormat = "GEO_TIFF",
      fileNamePrefix = paste0(sp_code$getInfo(), "_geemax_", names(pred_assets)[i])
    )
    
    task_img$start()
    ee_monitoring(task_img)
  }
  
  
  # Save objects for further processing ----
  # Save best tune
  bt_vector <- ee_table_to_drive(
    collection = best_tune,
    fileFormat = "CSV",
    fileNamePrefix = paste0(sp_code, "_geemax_besttune")
  )
  bt_vector$start()
  ee_monitoring(bt_vector)
  
  # Save cross-validation metrics
  finalcv_vector <- ee_table_to_drive(
    collection = final_m_cv,
    fileFormat = "CSV",
    fileNamePrefix = paste0(sp_code, "_geemax_crossvalid")
  )
  finalcv_vector$start()
  ee_monitoring(finalcv_vector)
  
  # Save full model metrics 
  fmod_vector <- ee_table_to_drive(
    collection = full_data_pred,
    fileFormat = "CSV",
    fileNamePrefix = paste0(sp_code, "_geemax_fullmodel")
  )
  fmod_vector$start()
  ee_monitoring(fmod_vector)
  
  # Save training points (optional)
  tdat_vector <- ee_table_to_drive(
    collection = all_points_with_grid,
    fileFormat = "CSV",
    fileNamePrefix = paste0(sp_code, "_geemax_trainingdata")
  )
  tdat_vector$start()
  ee_monitoring(tdat_vector)
  
  return(paste0(sp_code, "_done"))
  
}


# Additional functions for GEE modules
# Check if object is of type sdm_dat_gee
#' @export
.gee_checks <- function(x) {
  if (class(x)[1] != "sdm_dat_gee") {
    cli::cli_abort("sdm_data should be of type sdm_dat_gee (generated with {.fun obissdm::mp_prepare_geedat})")
  } else {
    # Check if Earth Engine was correctly started
    obj_try <- try(ee$Number(1), silent = T)
    if (class(obj_try)[1] == "try-error") {
      cli::cli_abort("Unable to call rgee functions. Have you authenticated and started the session with {.fun obissdm::start_ee}?")
    }
  }
  return(invisible(NULL))
}



#' Start rgee (Google Earth Engine) session
#'
#' @param init_only if \code{TRUE}, just initialize session. Useful if you already authenticated but restarted the R session
#' @param ... optional parameters passed to [rgee::ee_Authenticate())]
#'
#' @return nothing, just starts session
#' @export
#'
#' @details
#' Soon more details.
#' 
#'
#' @examples
#' \dontrun{
#' start_ee()
#' }
start_ee <- function(...) {
  
  rgee_available <- require("rgee", quietly = T)
  
  if (!rgee_available) {
    if (interactive()) {
      cli::cli_alert_danger("Package {.pkg rgee} is not installed. Do you want to try to install it?")
      install_it <- tolower(readline("y/n \t"))
      install_it <- ifelse(install_it == "y" | install_it == "yes", TRUE, FALSE)
      if (install_it) {
        install.packages("rgee")
      } else {
        cli::cli_abort("Package {.pkg rgee} is not installed. Install it before using GEE modules.")
      }
    } else {
      cli::cli_abort("Package {.pkg rgee} is not installed. Install it before using GEE modules.")
    }
  }
  
  cli::cli_alert_info("Starting authentication process")
  # # Authenticate to Google Earth Engine
  # ee_Authenticate(...)
  
  cli::cli_alert_info("Initializing session")
  # Initialize the Google Earth Engine session
  ee_Initialize()
  
  return(invisible(NULL))
  
}