############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################


# Main function(s)----

#' Fit (tune) species distribution model
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param sdm_method algorithm for model fiting. See [sdm_options()] for available
#' @param options list with options from [sdm_options()] to override the default
#' @param tune_blocks which blocks to use for tuning
#' @param metric which metric to use for model optimization
#' @param verbose if `TRUE` print messages
#' @param ... additional parameters (none implemented)
#'
#' @return `sdm_result` object containing the model
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_fit(sp_data)
#' }
sdm_fit <- function(sdm_data,
                    sdm_method = "maxent",
                    options = NULL,
                    tune_blocks = "spatial_grid",
                    metric = "cbi",
                    verbose = TRUE,
                    ...) {
  
  if (!sdm_method %in% names(sdm_options())) {
    cli::cli_abort("sdm_method {.var {sdm_method}} is not a recognized algorithm.")
  }
  
  arguments <- list(
    sdm_data = sdm_data,
    options = options,
    verbose = verbose,
    tune_blocks = tune_blocks,
    metric = metric
  )
  
  model_fit <- switch(
    sdm_method,
    gam = rlang::exec(sdm_module_gam, !!!arguments),
    glm = rlang::exec(sdm_module_glm, !!!arguments),
    maxent = rlang::exec(sdm_module_maxent, !!!arguments),
    lasso = rlang::exec(sdm_module_lasso, !!!arguments),
    rf = rlang::exec(sdm_module_rf, !!!arguments),
    brt = rlang::exec(sdm_module_brt, !!!arguments),
    xgboost = rlang::exec(sdm_module_xgboost, !!!arguments),
    lgbm = rlang::exec(sdm_module_lgbm, !!!arguments),
  )
  
  return(model_fit)
  
}


#' Retrieve options used for model tuning
#'
#' @param sdm_method `NULL` for returning full list, or an algorithm name
#'
#' @return list of options for the chosen algorithm or all
#' @export
#'
#' @examples
#' # Print all information
#' sdm_options()
#' 
#' # Retrieve for a single algorithm
#' sdm_options("maxent")
#' 
#' # Change details of an algorithm
#' max_opt <- sdm_options("maxent")
#' max_opt$features <- c("lq", "lqh")
#' 
#' # Supply to fit function
#' \dontrun{
#' sdm_species <- sdm_fit(sp_data, options = max_opt)
#' }
#' 
sdm_options <- function(sdm_method = NULL) {
  
  opt <- list(
    # Maxent (maxnet)
    maxent = list(
      features = "lq",
      remult = seq(0.5, 2, 0.5)
    ),
    # LASSO
    lasso = list(
      method = "naive",
      alpha_param = c(0, 0.5, 1),
      total_area = NULL,
      weight_resp = FALSE
    ),
    # BRT
    brt = list(
      method = "naive",
      n_trees = 100,
      t_depth = 1,
      lr = 0.01,
      bag_fraction = 0.75,
      weight_resp = FALSE
    ),
    # RF
    rf = list(
      n_trees = 500,
      mtry = c("default", "double", "total"),
      method = "classification",
      type = "down-sampled"
    ),
    # GLM
    glm = list(
      method = "naive",
      total_area = NULL,
      weight_resp = TRUE,
      quadratic = TRUE
    ),
    # GAM
    gam = list(
      method = "naive",
      total_area = NULL,
      weight_resp = TRUE,
      k_val = c(5, 10),
      select = c(FALSE, TRUE)
    ),
    # XGBoost
    xgboost = list(
      shrinkage = c(0.1, 0.3, 0.5),
      gamma = seq(from = 0, to = 4, by = 2),
      depth = c(3, 5),
      rounds = c(10, 50, 100),
      scale_pos_weight = 1
    ),
    # LightGBM
    lgbm = list(
      objective = "regression",
      max_dep = seq(2, 4, 1), # Deeper, more complex models
      lr = c(0.1, 0.05, 0.01, 0.005, 0.001), # controls the step size, lower learns "more"
      weight_resp = FALSE,
      num_leaves = 31, #c(10, 30) # controls the maximum number of leaves. High, more complex but risky overfit
      feature_fraction = 1, #c(0.8) # controls the fraction of features to consider for each tree, can reduce overfit
      bagging_fraction = 1,#c(0.8) # fraction of data to be used for each iteration, reduces overfit
      early_stopping_round = 0 # number of rounds after which the training will stop if there's no improvement in the validation metric
    ),
    # ESM
    esm = list(
      features = "lq",
      remult = seq(2, 3, 1)
    )
  )
  
  if (!is.null(sdm_method)) {
    if (!sdm_method %in% names(opt)) {
      cli::cli_abort("sdm_method {.var {sdm_method}} is not a recognized algorithm. Run {.fun sdm_options} with no arguments to see available options.")
    }
  }
  
  if (!is.null(sdm_method)) {
    opt <- opt[[sdm_method]]
    
    class(opt) <- c("sdm_options", "sdm_sel_options", class(opt))
    attr(opt, "method") <- sdm_method
  } else {
    class(opt) <- c("sdm_options", class(opt))
  }
  
  return(opt)
}



#' Test multiple hypothesis
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param variables list of hypothesis
#' @param fast_mode if `TRUE`, only the first model is tuned, and subsequent hypothesis are
#' tested using the best tune parameters
#' @param sdm_method algorithm to fit
#' @param compare_metric which metric to compare between hypothesis
#' @param return_all if `TRUE` return all models (default is `FALSE`)
#' @param verbose turn on messages
#' @param options options for model fitting
#' @param ... additional parameters (none implemented)
#'
#' @return fitted model, with best model information
#' @export
#'
#' @examples
#' \dontrun{
#' best_hypothesis <- sdm_multhypo(sp_data, hypo_list)
#' }
sdm_multhypo <- function(sdm_data,
                         variables,
                         fast_mode = TRUE,
                         sdm_method = "maxent",
                         compare_metric = "cbi",
                         return_all = FALSE,
                         verbose = TRUE,
                         options = NULL,
                         ...) {
  
  if (!inherits(variables, "list")) {
    cli::cli_abort("{.arg variables} should be of type {.val list}.")
  } else if (is.null(names(variables))) {
    if (verbose) cli::cli_warn("{.arg variables} should be a named list. Using generic names instead.")
    
    names(variables) <- paste0("hypothesis", seq(1, length(variables)))
  }
  
  if (is.null(options)) {
    options <- sdm_options(sdm_method)
  }
  
  n_vars <- length(variables)
  
  if (verbose) cli::cli_alert_info("Testing {.val {n_vars}} hypothesis using {.var {sdm_method}} and metric {.var {compare_metric}}.")
  
  result_metric <- rep(NA, n_vars)
  
  to_keep <- list()
  
  for (hyp in 1:n_vars) {
    
    if (verbose) cli::cli_alert_info("Fitting model {hyp} out of {n_vars}.")
    
    new_sdm_data <- sdm_data
    new_sdm_data$training <- new_sdm_data$training[,c("presence", variables[[hyp]])]
    
    if (!is.null(new_sdm_data$eval_data)) {
      new_sdm_data$eval_data <- new_sdm_data$eval_data[,c("presence", variables[[hyp]])]
    }
    
    if (fast_mode) {
      
      if (verbose & hyp == 1) cli::cli_alert_info("Using {.var fast_mode}. Only first model will be tunned, and others will use best tuning parameters.")
      
      mod <- sdm_fit(new_sdm_data, sdm_method = sdm_method, verbose = verbose, options = options, ...)
      
      if (hyp == 1) {
        for (nopt in names(mod$parameters)) {
          options[[nopt]] <- mod$parameters[[nopt]]
        }
      }
      
    } else {
      mod <- sdm_fit(sdm_data, sdm_method = sdm_method, verbose = verbose, options = options, ...)
    }
    
    to_keep[[hyp]] <- mod
    
    result_metric[hyp] <- mean(mod$cv_metrics[,compare_metric], na.rm = T)
    
  }
  
  best_model <- which.max(result_metric)[1]
  
  names(result_metric) <- names(variables)
  
  return_obj <- list(
    model = NULL,
    method = sdm_method,
    best_model = names(variables)[best_model],
    best_variables = variables[[best_model]],
    metrics = result_metric
  )
  
  if (return_all) {
    return_obj$model <- to_keep
  } else {
    return_obj$model <- to_keep[[best_model]]
  }
  
  class(return_obj) <- c("sdm_mult_result", class(return_obj))
  
  return(return_obj)
  
}


# SDM modules ----

#' Fit MAXENT SDM using maxnet
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_maxent(sp_data)
#' }
sdm_module_maxent <- function(sdm_data, options = NULL, verbose = TRUE,
                              tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("maxent")
  }
  
  features <- options[["features"]]
  remult <- options[["remult"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  # Tune model
  # Create grid for tuning
  tune_grid <- expand.grid(
    remult = remult,
    features = features,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- .maxent_cv(p, dat, b_index,
                             features = tune_grid$features[k],
                             regmult = tune_grid$remult[k])
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    tune_block <- apply(tune_block, 2, mean, na.rm = T)
    
    tune_test[k] <- tune_block[metric]
  }
  
  # Get best tune
  best_tune <- tune_grid[which.max(tune_test)[1],]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  full_fit <- maxnet::maxnet(p = p,
                             data = dat,
                             f = maxnet::maxnet.formula(p = p,
                                                        data = dat,
                                                        classes = best_tune$features),
                             regmult = best_tune$remult,
                             addsamplestobackground = T) # Change to F
  
  pred_full <- predict(full_fit, dat, type = "cloglog")
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = "maxent",
    model = full_fit,
    variables = colnames(dat),
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = best_tune,
    cv_metrics = cv_results[[which.max(tune_test)[1]]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    pred_eval <- predict(full_fit, dat_eval, type = "cloglog")
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "Maxent model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.maxent_cv <- function(p, dat, blocks, features, regmult){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- dat[blocks != id,]
    
    mfit <- try(maxnet::maxnet(p = train_p,
                               data = train_dat,
                               f = maxnet::maxnet.formula(p = train_p,
                                                          data = train_dat,
                                                          classes = features),
                               regmult = regmult,
                               addsamplestobackground = T), # change to F
                silent = T)
    
    if (inherits(mfit, "try-error")) { # Still to be fixed/verified
      NULL
    } else {
      pred <- predict(mfit, test_dat, type = "cloglog")
      
      eval_metrics(test_p, pred)
    }
    
  })
  
  return(do.call("rbind", blocks_results))
  
}





#' Fit Elasticnet/LASSO SDM using glmnet
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_lasso(sp_data)
#' }
sdm_module_lasso <- function(sdm_data, options = NULL, verbose = TRUE,
                             tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("lasso")
  }
  
  method <- options[["method"]]
  alpha_param <- options[["alpha_param"]]
  total_area <- options[["total_area"]]
  weight_resp <- options[["weight_resp"]]
  
  if (method == "dwpr" & is.null(total_area)) {
    stop("When method is 'dwpr' total_area should be supplied.")
  }
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  if (method == "iwlr" | method == "dwpr") {
    weight_resp <- TRUE
    to_norm <- TRUE
    glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  } else {
    to_norm <- FALSE
    glmnet::glmnet.control(pmin = 1e-09, fdev = 1e-05)
  }
  
  meas <- "auc"
  fam <- "binomial"
  resp_vec <- p
  
  if (weight_resp) {
    if (method == "iwlr") {
      wt <- (1e3)^(1 - p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(p))
      wt[p == 0] <- total_area/sum(p == 0)
      resp_vec <- p/wt
      fam <- "poisson"
      meas <- "default"
    } else {
      pres <- length(p[p == 1])
      bkg <- length(p[p == 0])
      wt <- ifelse(p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  # Tune model
  .cat_sdm(verbose, "Tuning model")
  
  forms <- as.formula(paste("~ 1",
                            paste(colnames(dat), collapse = "+"),
                            paste(paste("I(", colnames(dat), "^2", ")", sep = ""), 
                                  collapse = "+"), sep = "+"))
  
  training_poly <- model.matrix(forms, data = dat) 
  training_poly <- training_poly[,-1]
  
  b_index <- sdm_data$blocks$folds[[tune_blocks]]
  
  if (length(alpha_param) > 1) {
    
    alpha_error <- rep(NA, length(alpha_param))
    
    for (ap in 1:length(alpha_param)) {
      lasso_cv_a <- try(glmnet::cv.glmnet(x = training_poly,
                                          y = resp_vec,
                                          family = fam,
                                          alpha = alpha_param[ap],
                                          weights = wt,
                                          nfolds = length(unique(b_index)),
                                          foldid = b_index,
                                          type.measure = meas))
      
      if (!inherits(lasso_cv_a, "try-error")) {
        alpha_error[ap] <- lasso_cv_a$cvm[lasso_cv_a$index["1se",]]
      } else {
        alpha_error[ap] <- NA
      }
    }
    
    if (all(is.na(alpha_error))) {
      stop("Failed to fit glmnet with the supplied alpha parameters. Try a different range of values, or 1 for LASSO and 0 for Ridge.\n")
    }
    
    alpha_error <- data.frame(error = alpha_error, aparam = alpha_param)
    alpha_error <- alpha_error[!is.na(alpha_error$error), ]
    final_alpha_param <- alpha_error$aparam[which.min(alpha_error$error)]
  } else {
    final_alpha_param <- alpha_param
  }
  
  lasso_cv <- glmnet::cv.glmnet(x = training_poly,
                                y = resp_vec,
                                family = fam,
                                alpha = final_alpha_param,
                                weights = wt,
                                nfolds = length(unique(b_index)),
                                foldid = b_index,
                                type.measure = meas,
                                keep = TRUE)
  
  timings <- .get_time(timings, "tuning")
  
  best_lambda <- lasso_cv$lambda.1se
  
  # Get data for CV
  pred_test <- lasso_cv$fit.preval[,lasso_cv$index["1se",]]
  
  if (fam == "binomial") {
    pred_test <- 1/(1+exp(-pred_test))
  } else {
    pred_test <- exp(pred_test)
  }
  
  if (to_norm) {
    pred_test <- .normalize_res(pred_test)
  }
  
  cv_results <- lapply(unique(b_index), function(fid){
    eval_metrics(p[b_index == fid], pred_test[b_index == fid])
  })
  cv_results <- as.data.frame(do.call("rbind", cv_results))
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  pred_full <- predict(lasso_cv, training_poly, s = c("lambda.1se"), type = "response")
  
  if (to_norm) {
    pred_final <- .normalize_res(pred_full)
  }
  
  metrics_full <- eval_metrics(p, pred_full[,1])
  
  lasso_cv$fit.preval <- NULL
  
  timings <- .get_time(timings, "evaluate final")
  
  # Get model type to save
  if (final_alpha_param > 0 & final_alpha_param < 1) {
    model_type <- "elasticnet"
  } else if (final_alpha_param == 1) {
    model_type <- "lasso"
  } else {
    model_type <- "ridge"
  }
  
  result <- list(
    name = model_type,
    model = lasso_cv,
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(method = method,
                      alpha_param = final_alpha_param,
                      total_area = total_area,
                      weight_resp = weight_resp,
                      lambda_1se = lasso_cv$lambda.1se),
    cv_metrics = cv_results,
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    eval_poly <- model.matrix(forms, data = dat_eval) 
    eval_poly <- eval_poly[,-1]
    
    pred_eval <- predict(lasso_cv, eval_poly, s = c("lambda.1se"), type = "response")
    
    if (to_norm) {
      pred_eval <- .normalize_res(pred_eval)
    }
    
    eval_metric <- eval_metrics(p_eval, pred_eval[,1])
    
    result$eval_metrics <- eval_metric
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "LASSO model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}



#' Fit Boosted Regression Trees SDM using gbm
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_brt(sp_data)
#' }
sdm_module_brt <- function(sdm_data, options = NULL, verbose = TRUE,
                           tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("brt")
  }
  
  method <- options[["method"]]
  n_trees <- options[["n_trees"]]
  t_depth <- options[["t_depth"]]
  lr <- options[["lr"]]
  bag_fraction <- options[["bag_fraction"]]
  weight_resp <- options[["weight_resp"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training
  
  # Settings
  if (method == "iwlr" | method == "dwpr") {
    weight_resp <- TRUE
    to_norm <- TRUE
  } else {
    to_norm <- FALSE
  }
  
  if (weight_resp) {
    if (method == "iwlr") {
      wt <- (10^6)^(1 - p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(p))
      wt[p == 0] <- total_area/sum(p == 0)
      dat$presence <- p/wt
    } else {
      pres <- length(p[p == 1])
      bkg <- length(p[p == 0])
      wt <- ifelse(p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  # Tune model
  # Create grid for tuning
  tune_grid <- expand.grid(
    lr = lr,
    n_trees = n_trees,
    t_depth = t_depth,
    bag_fraction = bag_fraction,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))

    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- try(.brt_cv(p, dat, b_index,
                              n_trees = tune_grid$n_trees[k],
                              t_depth = tune_grid$t_depth[k], 
                              lr = tune_grid$lr[k],
                              bag_fraction = tune_grid$bag_fraction[k],
                              wt = wt, to_norm = to_norm, method = method),
                      silent = T)
    
    if (inherits(tune_block, "try-error")) tune_block <- NULL
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    if (!is.null(tune_block)) {
      tune_block <- apply(tune_block, 2, mean, na.rm = T)
      
      tune_test[k] <- tune_block[metric]
    } else {
      tune_test[k] <- NA
    }
  }
  
  # Get best tune
  bfit <- which.max(tune_test)[1]
  best_tune <- tune_grid[bfit,]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  full_fit <- gbm::gbm(presence ~ ., 
                       distribution = "poisson",
                       data = dat,
                       weights = wt,
                       n.trees = best_tune$n_trees,
                       interaction.depth = best_tune$t_depth,
                       shrinkage = best_tune$lr, 
                       bag.fraction = best_tune$bag_fraction,
                       cv.folds = 0,
                       n.minobsinnode = 5, 
                       verbose = FALSE)
  
  pred_full <- predict(full_fit, dat, type = "response")
  
  if (to_norm) pred_full <- .normalize_res(pred_full)
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = "brt",
    model = full_fit,
    variables = colnames(dat)[colnames(dat) != "presence"],
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(
      method = method,
      n_trees = best_tune$n_trees,
      t_depth = best_tune$t_depth,
      lr = best_tune$lr,
      bag_fraction = best_tune$bag_fraction,
      weight_resp = weight_resp
    ),
    cv_metrics = cv_results[[bfit]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    pred_eval <- predict(full_fit, dat_eval, type = "response")
    
    if (to_norm) pred_eval <- .normalize_res(pred_eval)
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "BRT model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.brt_cv <- function(p, dat, blocks, n_trees, t_depth, lr, bag_fraction,
                    wt, method, to_norm){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- dat[blocks != id,]
    
    nwt <- wt[blocks != id]
    
    if (method == "dwpr") {
      mfit <- gbm::gbm(presence ~ ., 
                  distribution = "poisson",
                  data = train_dat,
                  weights = nwt,
                  n.trees = n_trees,
                  interaction.depth = t_depth,
                  shrinkage = lr, 
                  bag.fraction = bag_fraction,
                  cv.folds = 0,
                  n.minobsinnode = 5, 
                  verbose = FALSE)
    } else {
      mfit <- gbm::gbm(presence ~ ., 
                  distribution = "bernoulli",
                  data = train_dat,
                  weights = nwt,
                  n.trees = n_trees,
                  interaction.depth = t_depth,
                  shrinkage = lr, 
                  bag.fraction = bag_fraction,
                  cv.folds = 0,
                  n.minobsinnode = 5, 
                  verbose = FALSE)
    }
    
    pred <- predict(mfit, test_dat, type = "response")
    
    if (to_norm) pred <- .normalize_res(pred)

    eval_metrics(test_p, pred)
    
  })
  
  return(do.call("rbind", blocks_results))
  
}





#' Fit Random Forest SDM using randomForest
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_rf(sp_data)
#' }
sdm_module_rf <- function(sdm_data, options = NULL, verbose = TRUE,
                          tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("rf")
  }
  
  n_trees <- options[["n_trees"]]
  method <- options[["method"]]
  type <- options[["type"]]
  mtry <- options[["mtry"]]
  mtry_names <- mtry
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training
  
  # Check if regression or classification
  if (method != "regression" | type == "down-sampled") {
    dat$presence <- as.factor(dat$presence)
  }
  
  # Check mtry
  if (any(mtry != "default")) {
    preds <- (ncol(dat) - 1)
    mtry <- unlist(lapply(mtry, function(x){
      if (method == "classification") {
        switch(x,
               default = floor(sqrt(preds)),
               double = floor(sqrt(preds)) * 2,
               total = preds
        )
      } else {
        switch(x,
               default = max(floor(ncol(preds)/3), 1),
               double = max(floor(ncol(preds)/3), 1) * 2,
               total = preds
        )
      }
    }))
    if (any(mtry > preds)) {
      mtry <- mtry[-which(mtry > preds)]
    }
  } else {
    preds <- (ncol(dat) - 1)
    if (method == "classification") {
      mtry <- floor(sqrt(preds))
    } else {
      mtry <- max(floor(ncol(preds)/3), 1)
    }
  }
  
  # Tune model
  # Create grid for tuning
  tune_grid <- expand.grid(
    n_trees = n_trees,
    mtry = mtry,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- .rf_cv(p, dat, b_index,
                         ntrees = tune_grid$n_trees[k],
                         mtry = tune_grid$mtry[k],
                         type = type,
                         method = method)
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    tune_block <- apply(tune_block, 2, mean, na.rm = T)
    
    tune_test[k] <- tune_block[metric]
  }
  
  # Get best tune
  best_tune <- tune_grid[which.max(tune_test)[1],]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  if (type == "down-sampled") {
    pres <- sum(p) 
    smpsize <- c("0" = pres, "1" = pres)
  }
  
  full_fit <- randomForest::randomForest(formula = presence ~ .,
                                         data = dat,
                                         ntree = best_tune$n_trees, 
                                         mtry = best_tune$mtry,
                                         sampsize = smpsize,
                                         replace = TRUE)
  
  if (type == "down-sampled" | method == "classification") {
    pred_full <- predict(full_fit, dat, type = "prob")
    pred_full <- pred_full[,2]
  } else {
    pred_full <- predict(full_fit, dat, type = "response")
  }
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = paste0("rf_",
                  method, "_",
                  ifelse(type == "down-sampled", "ds", "normal")),
    model = full_fit,
    variables = colnames(dat)[colnames(dat) != "presence"],
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(
      n_trees = best_tune$n_trees,
      mtry = mtry_names[which(mtry == best_tune$mtry)],
      method = method,
      type = type
    ),
    cv_metrics = cv_results[[which.max(tune_test)[1]]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    if (type == "down-sampled" | method == "classification") {
      pred_eval <- predict(full_fit, dat_eval, type = "prob")
      pred_eval <- pred_eval[,2]
    } else {
      pred_eval <- predict(full_fit, dat_eval, type = "response")
    }
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "RF model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.rf_cv <- function(p, dat, blocks, ntrees, mtry, type, method){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- dat[blocks != id,]
    
    if (type == "down-sampled" | method == "classification") {
      pres <- sum(train_p) 
      smpsize <- c("0" = pres, "1" = pres)
    }
    
    mfit <- randomForest::randomForest(formula = presence ~ .,
                                       data = train_dat,
                                       ntree = ntrees, 
                                       mtry = mtry,
                                       sampsize = smpsize,
                                       replace = TRUE)
    
    if (type == "down-sampled" | method == "classification") {
      pred <- predict(mfit, test_dat, type = "prob")
      pred <- pred[,2]
    } else {
      pred <- predict(mfit, test_dat, type = "response")
    }
    
    eval_metrics(test_p, pred)
    
  })
  
  return(do.call("rbind", blocks_results))
  
}




#' Fit GLM SDM
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_glm(sp_data)
#' }
sdm_module_glm <- function(sdm_data, options = NULL, verbose = TRUE,
                           tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("glm")
  }
  
  method <- options[["method"]]
  weight_resp <- options[["weight_resp"]]
  quadratic <- options[["quadratic"]]
  total_area <- options[["total_area"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training
  
  # Settings
  if (method == "iwlr" | method == "dwpr") {
    weight_resp <- TRUE
    to_norm <- TRUE
  } else {
    to_norm <- FALSE
  }
  
  if (weight_resp) {
    if (method == "iwlr") {
      wt <- (10^6)^(1 - p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(p))
      wt[p == 0] <- total_area/sum(p == 0)
      dat$presence <- p/wt
    } else {
      pres <- length(p[p == 1])
      bkg <- length(p[p == 0])
      wt <- ifelse(p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  # Tune model
  .cat_sdm(verbose, "Tuning model")
  
  var_names <- colnames(dat)[colnames(dat) != "presence"]
  
  forms <- as.formula(
    paste("presence ~",
          paste(
            var_names,
            collapse = "+"
          ),
          ifelse(quadratic,  paste0("+ ", paste0(
            "I(", var_names, "^2)",
            collapse = " + "
          ), "")))
  )
  
  b_index <- sdm_data$blocks$folds[[tune_blocks]]
  
  cv_results <- .glm_cv(p, dat, b_index, forms, wt = wt, method, to_norm)
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  if (method == "dwpr") {
    full_fit <- glm(forms, family = poisson(), data = dat, weights = wt)
  } else {
    full_fit <- glm(forms, family = binomial(), data = dat, weights = wt)
  }
  
  pred_full <- predict(full_fit, dat, type = "response")
  
  if (to_norm) pred_full <- .normalize_res(pred_full)
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = paste0("glm_", ifelse(method == "naive", "normal", method)),
    model = full_fit,
    variables = var_names,
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(weight_resp = weight_resp,
                      total_area = total_area,
                      quadratic = quadratic,
                      method = method),
    cv_metrics = as.data.frame(cv_results),
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    pred_eval <- predict(full_fit, dat_eval, type = "response")
    
    if (to_norm) pred_eval <- .normalize_res(pred_eval)
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "GLM model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.glm_cv <- function(p, dat, blocks, forms, wt, method, to_norm){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- dat[blocks != id,]
    
    nwt <- wt[blocks != id]
    
    train_dat$nwt <- nwt
    
    rm(nwt)
    
    if (method == "dwpr") {
      mfit <- glm(forms, family = poisson(), data = train_dat, weights = nwt)
    } else {
      mfit <- glm(forms, family = binomial(), data = train_dat, weights = nwt)
    }
    
    pred <- predict(mfit, test_dat, type = "response")
    
    if (to_norm) pred <- .normalize_res(pred)
    
    eval_metrics(test_p, pred)
    
  })
  
  return(do.call("rbind", blocks_results))
  
}


#' Fit GAM SDM using mgcv
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_gam(sp_data)
#' }
sdm_module_gam <- function(sdm_data, options = NULL, verbose = TRUE,
                           tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("gam")
  }
  
  method <- options[["method"]]
  weight_resp <- options[["weight_resp"]]
  k_val <- options[["k_val"]]
  total_area <- options[["total_area"]]
  select <- options[["select"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training
  
  # Settings
  if (method == "iwlr") {
    fam <- binomial(link = "logit")
    weight_resp <- TRUE
    to_norm <- TRUE
  } else if (method == "dwpr") {
    fam <- poisson()
    weight_resp <- TRUE
    to_norm <- TRUE
  } else {
    fam <- binomial(link = "logit")
    to_norm <- FALSE
  }
  
  if (weight_resp) {
    if (method == "iwlr") {
      wt <- (10^6)^(1 - p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(p))
      wt[p == 0] <- total_area/sum(p == 0)
      dat$presence <- p/wt
    } else {
      pres <- length(p[p == 1])
      bkg <- length(p[p == 0])
      wt <- ifelse(p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  # Tune model
  .cat_sdm(verbose, "Tuning model")
  
  var_names <- colnames(dat)[colnames(dat) != "presence"]
  
  
  # Create grid for tuning
  tune_grid <- expand.grid(
    k_val = k_val,
    select = select,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))
    
    tune_forms <- as.formula(
      paste("presence ~", paste(
        "s(", var_names, ", k=", tune_grid$k_val[k], ", bs='cr')",
        collapse = "+"
      ))
    )
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- .gam_cv(p, dat, b_index,
                          forms = tune_forms,
                          wt = wt, family = fam,
                          to_norm = to_norm,
                          select = tune_grid$select[k])
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    tune_block <- apply(tune_block, 2, mean, na.rm = T)
    
    tune_test[k] <- tune_block[metric]
  }
  
  # Get best tune
  best_tune <- tune_grid[which.max(tune_test)[1],]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  forms <- as.formula(
    paste("presence ~", paste(
      "s(", var_names, ", k=", best_tune$k_val, ", bs='cr')",
      collapse = "+"
    ))
  )
  
  full_fit <- mgcv::bam(forms, family = fam,
                        data = dat,
                        weights = wt,
                        method = "fREML",
                        select = best_tune$select)
  
  pred_full <- predict(full_fit, dat, type = "response")
  
  if (to_norm) pred_full <- .normalize_res(pred_full)
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = paste0("gam_", ifelse(method == "naive", "normal", method)),
    model = full_fit,
    variables = var_names,
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(weight_resp = weight_resp,
                      total_area = total_area,
                      k_val = best_tune$k_val,
                      select = best_tune$select,
                      method = method),
    cv_metrics = cv_results[[which.max(tune_test)[1]]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    pred_eval <- predict(full_fit, dat_eval, type = "response")
    
    if (to_norm) pred_eval <- .normalize_res(pred_eval)
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "GAM model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.gam_cv <- function(p, dat, blocks, forms, wt, family, to_norm, select){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- dat[blocks != id,]
    
    nwt <- wt[blocks != id]
    
    train_dat$nwt <- nwt
    
    rm(nwt)
    
    mfit <- mgcv::bam(forms, family = family, data = train_dat,
                      weights = nwt,
                      method = "fREML",
                      select = select)
    
    pred <- predict(mfit, test_dat, type = "response")
    
    if (to_norm) pred <- .normalize_res(pred)
    
    eval_metrics(test_p, pred)
    
  })
  
  return(do.call("rbind", blocks_results))
  
}





#' Fit XGBoost SDM
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_xgboost(sp_data)
#' }
sdm_module_xgboost <- function(sdm_data, options = NULL, verbose = TRUE,
                               tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("xgboost")
  }
  
  shrinkage <- options[["shrinkage"]]
  gamma <- options[["gamma"]]
  depth <- options[["depth"]]
  rounds <- options[["rounds"]]
  scale_pos_weight <- options[["scale_pos_weight"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  
  # Tune model
  # Create grid for tuning
  tune_grid <- expand.grid(
    shrinkage = shrinkage,
    gamma = gamma,
    depth = depth,
    rounds = rounds,
    scale_pos_weight = scale_pos_weight,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- .xgb_cv(p, dat, b_index,
                          shrinkage = tune_grid$shrinkage[k],
                          gamma = tune_grid$gamma[k],
                          depth = tune_grid$depth[k],
                          rounds = tune_grid$rounds[k],
                          scale_pos_weight = tune_grid$scale_pos_weight[k])
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    tune_block <- apply(tune_block, 2, mean, na.rm = T)
    
    tune_test[k] <- tune_block[metric]
  }
  
  # Get best tune
  best_tune <- tune_grid[which.max(tune_test)[1],]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  full_fit <- xgboost::xgboost(
    data = as.matrix(dat),
    label = p,
    max_depth = best_tune$depth,
    gamma = best_tune$gamma, 
    scale_pos_weight = best_tune$scale_pos_weight,
    nrounds = best_tune$rounds,
    learning_rate = best_tune$shrinkage,
    verbose = 0,
    nthread = 1, 
    objective = "binary:logistic")
  
  pred_full <- predict(full_fit, as.matrix(dat), type = "response")
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = "xgboost",
    model = full_fit,
    variables = colnames(dat)[colnames(dat) != "presence"],
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(
      depth = best_tune$depth,
      gamma = best_tune$gamma, 
      scale_pos_weight = best_tune$scale_pos_weight,
      rounds = best_tune$rounds,
      shrinkage = best_tune$shrinkage,
      objective = "binary:logistic"
    ),
    cv_metrics = cv_results[[which.max(tune_test)[1]]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    pred_eval <- predict(full_fit, as.matrix(dat_eval), type = "response")
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "XGBoost model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.xgb_cv <- function(p, dat, blocks, shrinkage, gamma, depth, rounds,
                    scale_pos_weight){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- dat[blocks != id,]
    
    mfit <-  xgboost::xgboost(
      data = as.matrix(train_dat),
      label = train_p,
      max_depth = depth,
      gamma = gamma, 
      scale_pos_weight = scale_pos_weight,
      nrounds = rounds,
      learning_rate = shrinkage,
      verbose = 0,
      nthread = 1, 
      objective="binary:logistic")
    
    pred <- predict(mfit, as.matrix(test_dat), type = "response")
    
    eval_metrics(test_p, pred)
    
  })
  
  return(do.call("rbind", blocks_results))
  
}




#' Fit LightGBM (gradient-boosting machine) SDM
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_lgbm(sp_data)
#' }
sdm_module_lgbm <- function(sdm_data, options = NULL, verbose = TRUE,
                            tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("lgbm")
  }
  
  objective <- options[["objective"]]
  max_dep <- options[["max_dep"]]
  lr <- options[["lr"]]
  num_leaves <- options[["num_leaves"]]
  feature_fraction <- options[["feature_fraction"]]
  bagging_fraction <- options[["bagging_fraction"]]
  early_stopping_round <- options[["early_stopping_round"]]
  weight_resp <- options[["weight_resp"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  # Settings
  if (weight_resp) {
    pres <- length(train_p[train_p == 1])
    bkg <- length(train_p[train_p == 0])
    wt <- ifelse(train_p == 1, 1, pres / bkg)
  } else {
    wt <- NULL
  }
  
  # Tune model
  # Create grid for tuning
  tune_grid <- expand.grid(
    objective = objective,
    max_dep = max_dep,
    lr = lr,
    num_leaves = num_leaves,
    feature_fraction = feature_fraction,
    bagging_fraction = bagging_fraction,
    early_stopping_round = early_stopping_round,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- .lgbm_cv(p, dat, b_index,
                           objective = tune_grid$objective[k],
                           max_dep = tune_grid$max_dep[k],
                           lr = tune_grid$lr[k],
                           num_leaves = tune_grid$num_leaves[k],
                           feature_fraction = tune_grid$feature_fraction[k],
                           bagging_fraction = tune_grid$bagging_fraction[k],
                           early_stopping_round = tune_grid$early_stopping_round[k],
                           wt = wt)
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    tune_block <- apply(tune_block, 2, mean, na.rm = T)
    
    tune_test[k] <- tune_block[metric]
  }
  
  # Get best tune
  best_tune <- tune_grid[which.max(tune_test)[1],]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  to_hold <- sample(1:nrow(dat), nrow(dat)*0.2)
  
  to_include <- 1:nrow(dat)
  to_include <- to_include[-to_hold]
  
  full_fit <- lightgbm::lgb.train(
    params = list(
      objective = best_tune$objective,
      metric = "auc",
      max_depth = best_tune$max_dep,
      num_iterations = 1000,
      early_stopping_rounds = best_tune$early_stopping_round,
      learning_rate = best_tune$lr,
      bagging_fraction = best_tune$bagging_fraction,
      bagging_freq = 1,
      feature_fraction = best_tune$feature_fraction,
      num_leaves = best_tune$num_leaves
      #is_unbalance = T
      #feature_fraction = .8
    ),
    valids = list(test = lightgbm::lgb.Dataset(as.matrix(dat)[to_hold,], label = p[to_hold], weight = wt[to_hold])),
    data = lightgbm::lgb.Dataset(as.matrix(dat)[to_include,], label = p[to_include], weight = wt[to_include])
  )
  
  if (objective == "regression") {
    pred_full <- predict(full_fit, as.matrix(dat), type = "response")
    pred_full <- exp(pred_full)/(1+exp(pred_full))
  } else {
    pred_full <- predict(full_fit, as.matrix(dat), type = "response")
  }
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = "lgbm",
    model = full_fit,
    variables = colnames(dat)[colnames(dat) != "presence"],
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = list(
      objective = best_tune$objective,
      max_dep = best_tune$max_dep,
      lr = best_tune$lr,
      num_leaves = best_tune$num_leaves,
      feature_fraction = best_tune$feature_fraction,
      bagging_fraction = best_tune$bagging_fraction,
      early_stopping_round = best_tune$early_stopping_round,
      weight_resp = weight_resp
    ),
    cv_metrics = cv_results[[which.max(tune_test)[1]]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    if (objective == "regression") {
      pred_eval <- predict(full_fit, as.matrix(dat_eval), type = "response")
      pred_eval <- exp(pred_eval)/(1+exp(pred_eval))
    } else {
      pred_eval <- predict(full_fit, as.matrix(dat_eval), type = "response")
    }
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "LightGBM model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}

#' @export
.lgbm_cv <- function(p, dat, blocks, objective, max_dep,
                     lr, num_leaves, feature_fraction, bagging_fraction,
                     early_stopping_round, wt){
  
  blocks_results <- lapply(1:length(unique(blocks)), function(id){
    
    test_p <- p[blocks == id]
    test_dat <- dat[blocks == id,]
    
    train_p <- p[blocks != id]
    train_dat <- as.matrix(dat[blocks != id,])
    
    nwt <- wt[blocks != id]
    
    to_hold <- sample(1:nrow(train_dat), nrow(train_dat)*0.2)
    
    to_include <- 1:nrow(train_dat)
    to_include <- to_include[-to_hold]
    
    mfit <- lightgbm::lgb.train(
      params = list(
        objective = objective,
        metric = "auc",
        max_depth = max_dep,
        num_iterations = 1000,
        early_stopping_rounds = early_stopping_round,
        learning_rate = lr,
        bagging_fraction = bagging_fraction,
        bagging_freq = 1,
        feature_fraction = feature_fraction,
        num_leaves = num_leaves
        #is_unbalance = T
        #feature_fraction = .8
      ),
      valids = list(test = lightgbm::lgb.Dataset(train_dat[to_hold,], label = train_p[to_hold], weight = nwt[to_hold])),
      data = lightgbm::lgb.Dataset(train_dat[to_include,], label = train_p[to_include], weight = nwt[to_include])
    )
    
    if (objective == "regression") {
      pred <- predict(mfit, as.matrix(test_dat))
      pred <- exp(pred)/(1+exp(pred))
    } else {
      pred <- predict(mfit, as.matrix(test_dat))
    }
    
    eval_metrics(test_p, pred)
    
  })
  
  return(do.call("rbind", blocks_results))
  
}


#' Fit Ensemble of Small Models using maxnet
#'
#' @param sdm_data `sdm_dat` object containing occurrences and environmental data
#' @param options list with options from [sdm_options()] to override the default
#' @param verbose if `TRUE` display messages
#' @param tune_blocks which blocks to use for cross-validation
#' @param metric which metric to optimize in tuning using cross-validation
#'
#' @return fitted model (`sdm_result` object)
#' @export
#' 
#' @references Breiner, F. T., Nobis, M. P., Bergamini, A., & Guisan, A. (2018). 
#' Optimizing ensembles of small models for predicting the distribution of species 
#' with few occurrences. In N. Isaac (Ed.), Methods in Ecology and Evolution 
#' (Vol. 9, Issue 4, pp. 802–808). Wiley. https://doi.org/10.1111/2041-210x.12957
#'
#' @examples
#' \dontrun{
#' sdm_species <- sdm_module_esm(sp_data)
#' }
sdm_module_esm <- function(sdm_data, options = NULL, verbose = TRUE,
                              tune_blocks = "spatial_grid", metric = "cbi") {
  
  # Checkings
  .check_type(sdm_data)
  .cat_sdm(verbose, "Preparing data")
  timings <- .get_time()
  
  # Get options
  if (is.null(options)) {
    options <- sdm_options("esm")
  }
  
  features <- options[["features"]]
  remult <- options[["remult"]]
  
  # Separate data
  p <- sdm_data$training$presence
  dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  # Tune model
  # Create grid for tuning
  tune_grid <- expand.grid(
    remult = remult,
    features = features,
    stringsAsFactors = FALSE
  )
  
  .cat_sdm(verbose, "Tuning model")
  
  tune_test <- rep(NA, nrow(tune_grid))
  cv_results <- list()
  
  for (k in 1:nrow(tune_grid)) {
    
    .cat_sdm(verbose, glue::glue("Tuning option {k} out of {nrow(tune_grid)}"))
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    tune_block <- .maxent_cv(p, dat, b_index,
                             features = tune_grid$features[k],
                             regmult = tune_grid$remult[k])
    
    cv_results[[k]] <- as.data.frame(tune_block)
    
    tune_block <- apply(tune_block, 2, mean, na.rm = T)
    
    tune_test[k] <- tune_block[metric]
  }
  
  # Get best tune
  best_tune <- tune_grid[which.max(tune_test)[1],]
  
  timings <- .get_time(timings, "tuning")
  
  # Fit full model
  .cat_sdm(verbose, "Training and evaluating final model")
  
  full_fit <- maxnet::maxnet(p = p,
                             data = dat,
                             f = maxnet::maxnet.formula(p = p,
                                                        data = dat,
                                                        classes = best_tune$features),
                             regmult = best_tune$remult,
                             addsamplestobackground = T) # Change to F
  
  pred_full <- predict(full_fit, dat, type = "cloglog")
  
  metrics_full <- eval_metrics(p, pred_full)
  
  timings <- .get_time(timings, "evaluate final")
  
  # Prepare returning object
  result <- list(
    name = "maxent",
    model = full_fit,
    variables = colnames(dat),
    n_pts = c(presence = sum(p),
              background = (length(p) - sum(p))),
    timings = timings,
    cv_method = tune_blocks,
    parameters = best_tune,
    cv_metrics = cv_results[[which.max(tune_test)[1]]],
    full_metrics = metrics_full,
    eval_metrics = NULL
  )
  
  
  # If evaluation dataset is available, evaluate
  if (!is.null(sdm_data$eval_data)) {
    
    # Separate data
    p_eval <- sdm_data$eval_data$presence
    dat_eval <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    # Predict
    pred_eval <- predict(full_fit, dat_eval, type = "cloglog")
    
    metrics_eval <- eval_metrics(p_eval, pred_eval)
    
    result$eval_metrics <- metrics_eval
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
    
  }
  
  .cat_sdm(verbose, "Maxent model concluded", bg = T)
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}