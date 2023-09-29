############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################



#' SDM module: Generalized Linear Model
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param method which method to use for the GLM. One of normal, iwlr or dwpr
#'   (see details)
#' @param weight_resp wether to apply weights or not for the normal type (see
#'   details)
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#' @param total_area for the DWPR, a numeric value with the total area of the
#'   study region
#'
#' @details The following methods are available
#' - \code{normal}: a standard binomial GLM. If weight_resp is \code{TRUE},
#'  then the response is weighted in a way that the total weight of presences
#'  is the same as the weight of background points (naive weighting).
#' - \code{iwlr}: Infinitely Weighted Linear Regression
#' - \code{dwpr}: Down Weighted Poisson Regression
#'
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_glm(species_data)
#' }
sdm_module_glm <- function(sdm_data, method = "iwlr", weight_resp = TRUE,
                           tune_blocks = "spatial_grid", blocks_all = FALSE,
                           total_area = NULL) {
  
  if (method == "dwpr" & is.null(total_area)) {
    stop("When method is 'dwpr' total_area should be supplied.")
  }
  
  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training
  
  if (method == "iwlr" | method == "dwpr") {
    weight_resp <- TRUE
    to_norm <- TRUE
  } else {
    to_norm <- FALSE
  }
  
  if (weight_resp) {
    if (method == "iwlr") {
      wt <- (10^6)^(1 - train_p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(train_p))
      wt[train_p == 0] <- total_area/sum(train_p == 0)
      train_dat$presence <- train_p/wt
    } else {
      pres <- length(train_p[train_p == 1])
      bkg <- length(train_p[train_p == 0])
      wt <- ifelse(train_p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  ### Tune model
  
  # Create a function of the model
  glm_mod <- function(labels, training, weights = NULL) {
    
    var_names <- colnames(training)
    var_names <- var_names[!grepl("presence", var_names)]
    
    # forms <- as.formula(
    #   paste("presence ~", paste(
    #     var_names, "+ poly(", var_names, ", degree = 2)",
    #     collapse = "+"
    #   ))
    # )
    
    forms <- as.formula(
      paste("presence ~", paste(
        var_names,
        collapse = "+"
      ))
    )
    
    if (method == "dwpr") {
      glm_m <- glm(forms, family = poisson(), data = training, weights = weights)
    } else {
      glm_m <- glm(forms, family = binomial(), data = training, weights = weights)
    }
    
    eval(parse(text = paste0(
      "nscope <- list(",
      paste0("'", var_names, "' = ~ 1 +", var_names, "+ poly(", var_names, ", 2)", collapse = ","),
      ")"
    )))
    
    fmod <- gam::step.Gam(glm_m,
                  scope = nscope,
                  direction = "both",
                  data = training, 
                  trace = F)
    
    if (is.null(fmod)) {
      warning("step.Gam failed, returning the full model.")
      return(glm_m)
    } else {
      return(fmod)
    }
  }
  
  glm_sel_mod <- glm_mod(train_p, train_dat, wt)
  
  timings <- .get_time(timings, "tuning")
  
  ### Evaluate final model through CV
  
  if (blocks_all) {
    final_cv_metric <- list()
    
    for (b in 1:length(sdm_data$blocks$folds)) {
      
      b_index <- sdm_data$blocks$folds[[b]]
      
      final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                     folds = b_index, labels = train_p,
                                     training = train_dat, metrics = "all",
                                     normalize = to_norm,
                                     pred_type = "response", model_fun = glm_mod,
                                     weights = wt)
      
      final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
      
    }
    
    names(final_cv_metric) <- names(sdm_data$blocks$folds)
    cv_method <- names(sdm_data$blocks$folds)
    
  } else {
    
    final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                       folds = sdm_data$blocks$folds[[tune_blocks]], labels = train_p,
                       training = train_dat, metrics = "all",
                       normalize = to_norm,
                       pred_type = "response", model_fun = glm_mod,
                       weights = wt)
    
    final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
    
    cv_method <- tune_blocks
  }
  
  timings <- .get_time(timings, "cv")
  
  ### Train final model
  m_final <- glm_sel_mod
  
  ### Evaluate final model (full data)
  
  pred_final <- predict(m_final, train_dat, type = "response")
  
  if (to_norm) {
    pred_final <- .normalize_res(pred_final)
  }
  
  final_full_metric <- eval_metrics(train_p, as.vector(pred_final))
  
  timings <- .get_time(timings, "evaluate final")
  
  ### Prepare returning object
  result <- list(
    name = paste0("glm_", method),
    model = m_final,
    timings = timings,
    tune_cv_method = tune_blocks,
    cv_method = cv_method,
    parameters = as.character(m_final$formula)[3],
    cv_metrics = final_cv_metric,
    full_metrics = final_full_metric,
    eval_metrics = NULL
  )
  
  ### If testing dataset is available, use it to evaluate
  
  if (!is.null(sdm_data$eval_data)) {
    eval_p <- sdm_data$eval_data$presence
    eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    pred_eval <- predict(m_final, eval_dat, type = "response")
    
    if (to_norm) {
      pred_eval <- .normalize_res(pred_eval)
    }
    
    eval_metric <- eval_metrics(eval_p, as.vector(pred_eval))
    
    result$eval_metrics <- eval_metric
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
  }
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}




#' SDM module: Generalized Additive Model
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param method which method to use for the GAM. One of normal, iwlr or dwpr
#'   (see details)
#' @param weight_resp wether to apply weights or not for the normal type (see
#'   details)
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#' @param total_area for the DWPR, a numeric value with the total area of the
#'   study region
#'
#' @details The following methods are available
#' - \code{normal}: a standard binomial GAM. If weight_resp is \code{TRUE},
#'  then the response is weighted in a way that the total weight of presences
#'  is the same as the weight of background points (naive weighting).
#' - \code{iwlr}: Infinitely Weighted Linear Regression
#' - \code{dwpr}: Down Weighted Poisson Regression
#'
#' @note
#' Although the DWPR option is available for GAM, it will usually take long time
#' to run. Users are advised to avoid DWPR for this algorithm.
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_gam(species_data)
#' }
sdm_module_gam <- function(sdm_data, method = "iwlr", weight_resp = TRUE,
                           tune_blocks = "spatial_grid", blocks_all = FALSE,
                           total_area = NULL) {
  
  if (method == "dwpr" & is.null(total_area)) {
    stop("When method is 'dwpr' total_area should be supplied.")
  }
  
  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training
  
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
      wt <- (10^6)^(1 - train_p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(train_p))
      wt[train_p == 0] <- total_area/sum(train_p == 0)
      train_dat$presence <- train_p/wt
    } else {
      pres <- length(train_p[train_p == 1])
      bkg <- length(train_p[train_p == 0])
      wt <- ifelse(train_p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  ### Tune model
  s_param <- c(5, 10)
  
  auc_test <- rep(NA, length(s_param))
  
  # Create a function of the model
  gam_mod <- function(labels, training, sval, weights = NULL) {
    
    var_names <- colnames(training)
    var_names <- var_names[!grepl("presence", var_names)]
    
    forms <- as.formula(
      paste("presence ~", paste(
        "s(", var_names, ", k=", sval, ")",
        collapse = "+"
      ))
    )
    
    mgcv::gam(formula = forms,
              data = training,
              family = fam,
              weights = weights,
              method = "REML")
  }
  
  for (z in 1:length(auc_test)) {
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    auc_block <- cv_mod(n_folds = sdm_data$blocks$n,
                        folds = b_index, labels = train_p,
                        training = train_dat, metrics = "auc",
                        normalize = to_norm,
                        pred_type = "response", model_fun = gam_mod,
                        sval = s_param[z], weights = wt)
    
    auc_test[z] <- mean(unlist(auc_block), na.rm = T)
  }
  
  timings <- .get_time(timings, "tuning")
  
  ### Evaluate final model through CV
  
  if (blocks_all) {
    final_cv_metric <- list()
    
    for (b in 1:length(sdm_data$blocks$folds)) {
      
      b_index <- sdm_data$blocks$folds[[b]]
      
      final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                     folds = b_index, labels = train_p,
                                     training = train_dat, metrics = "all",
                                     normalize = to_norm,
                                     pred_type = "response", model_fun = gam_mod,
                                     sval = s_param[which.min(auc_test)], weights = wt)
      
      final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
      
    }
    
    names(final_cv_metric) <- names(sdm_data$blocks$folds)
    cv_method <- names(sdm_data$blocks$folds)
    
  } else {
    
    final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                       folds = b_index, labels = train_p,
                       training = train_dat, metrics = "all",
                       normalize = to_norm,
                       pred_type = "response", model_fun = gam_mod,
                       sval = s_param[which.min(auc_test)], weights = wt)
    
    final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
    
    cv_method <- tune_blocks
  }
  
  timings <- .get_time(timings, "cv")
  
  ### Train final model
  m_final <- gam_mod(labels = train_p, training = train_dat,
                     sval = s_param[which.min(auc_test)],
                     weights = wt)
  
  ### Evaluate final model (full data)
  
  pred_final <- predict(m_final, train_dat, type = "response")
  
  if (to_norm) {
    pred_final <- .normalize_res(pred_final)
  }
  
  final_full_metric <- eval_metrics(train_p, as.vector(pred_final))
  
  timings <- .get_time(timings, "evaluate final")
  
  ### Prepare returning object
  result <- list(
    name = paste0("gam_", method),
    model = m_final,
    timings = timings,
    tune_cv_method = tune_blocks,
    cv_method = cv_method,
    parameters = c(s = s_param[which.min(auc_test)]),
    cv_metrics = final_cv_metric,
    full_metrics = final_full_metric,
    eval_metrics = NULL
  )
  
  ### If testing dataset is available, use it to evaluate
  
  if (!is.null(sdm_data$eval_data)) {
    eval_p <- sdm_data$eval_data$presence
    eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    pred_eval <- predict(m_final, eval_dat, type = "response")
    
    if (to_norm) {
      pred_eval <- .normalize_res(pred_eval)
    }
    
    eval_metric <- eval_metrics(eval_p, as.vector(pred_eval))
    
    result$eval_metrics <- eval_metric
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
  }
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
  
}



#' SDM module: LASSO (Regularized Generalized Linear Models)
#' 
#' Implements the Lasso Regularized Generalized Linear Models through the
#' [glmnet] package.
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param weight_resp if weight_resp is \code{TRUE}, then the occurrence points
#'   are weighted according to weighting method (see below)
#' @param method the method for the weighting (if \code{weight_resp = TRUE}).
#'   Should be one of "naive", "iwlr" (for Infinitely Weighted Regression) or 
#'   "dwpr" (for Down-Weighted Poisson Regression). If "naive" then the response 
#'   is weighted in a way that the total weight of presences is the same as the 
#'   weight of background points
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#' @param total_area if weighting method is the DWPR (Down-Weighted Poisson 
#'   Regression), a numeric value with the total area of the study region
#'
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_lasso(species_data)
#' }
sdm_module_lasso <- function(sdm_data, weight_resp = TRUE, method = "naive",
                             tune_blocks = "spatial_grid", blocks_all = FALSE,
                             total_area = NULL) {
  
  if (method == "dwpr" & is.null(total_area)) {
    stop("When method is 'dwpr' total_area should be supplied.")
  }
  
  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
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
  resp_vec <- train_p
  
  if (weight_resp) {
    if (method == "iwlr") {
      wt <- (1e3)^(1 - train_p)
    } else if (method == "dwpr") {
      wt <- rep(1e-6, length(train_p))
      wt[train_p == 0] <- total_area/sum(train_p == 0)
      resp_vec <- train_p/wt
      fam <- "poisson"
      meas <- "default"
    } else {
      pres <- length(train_p[train_p == 1])
      bkg <- length(train_p[train_p == 0])
      wt <- ifelse(train_p == 1, 1, pres / bkg)
    }
  } else {
    wt <- NULL
  }
  
  ### Tune model
  forms <- as.formula(paste("~ 1",
                            paste(colnames(train_dat), collapse = "+"),
                            paste(paste("I(", colnames(train_dat), "^2", ")", sep = ""), 
                                  collapse = "+"), sep = "+"))
  
  training_poly <- model.matrix(forms, data = train_dat) 
  training_poly <- training_poly[,-1]
  
  lasso_cv <- glmnet::cv.glmnet(x = training_poly,
                                y = resp_vec,
                                family = fam,
                                alpha = 1,
                                weights = wt,
                                nfolds = 5,
                                foldid = sdm_data$blocks$folds[[tune_blocks]],
                                type.measure = meas)
  
  timings <- .get_time(timings, "tuning")
  
  # Create a function of the model
  glmnet_mod <- function(labels, training, weights) {
    if (method != "dwpr") {
      glmnet::glmnet(x = training,
                     y = labels,
                     family = fam,
                     alpha = 1,
                     weights = weights,
                     lambda = lasso_cv$lambda.1se)
    } else {
      glmnet::glmnet(x = training,
                     y = (labels / weights),
                     family = fam,
                     alpha = 1,
                     weights = weights,
                     lambda = lasso_cv$lambda.1se)
    }
  }
  
  ### Evaluate final model through CV
  
  if (blocks_all) {
    final_cv_metric <- list()
    
    for (b in 1:length(sdm_data$blocks$folds)) {
      
      b_index <- sdm_data$blocks$folds[[b]]
      
      final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                     folds = b_index, labels = train_p,
                                     training = training_poly, metrics = "all",
                                     pred_type = "response", model_fun = glmnet_mod,
                                     weights = wt, normalize = to_norm)
      
      final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
      
    }
    
    names(final_cv_metric) <- names(sdm_data$blocks$folds)
    cv_method <- names(sdm_data$blocks$folds)
    
  } else {
    
    final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                       folds = sdm_data$blocks$folds[[tune_blocks]],
                       labels = train_p,
                       training = training_poly, metrics = "all",
                       pred_type = "response", model_fun = glmnet_mod,
                       weights = wt, normalize = to_norm)
    
    final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
    
    cv_method <- tune_blocks
  }
  
  timings <- .get_time(timings, "cv")
  
  ### Train final model
  m_final <- glmnet::glmnet(x = training_poly,
                            y = resp_vec,
                            family = fam,
                            alpha = 1,
                            weights = wt,
                            lambda = lasso_cv$lambda.1se)
  
  ### Evaluate final model (full data)
  
  pred_final <- predict(m_final, training_poly, type = "response")
  
  if (to_norm) {
    pred_final <- .normalize_res(pred_final)
  }
  
  final_full_metric <- eval_metrics(train_p, pred_final[,1])
  
  timings <- .get_time(timings, "evaluate final")
  
  ### Prepare returning object
  result <- list(
    name = "lasso",
    model = m_final,
    timings = timings,
    tune_cv_method = tune_blocks,
    cv_method = cv_method,
    parameters = c("lambda_1se" = lasso_cv$lambda.1se),
    cv_metrics = final_cv_metric,
    full_metrics = final_full_metric,
    eval_metrics = NULL
  )
  
  ### If testing dataset is available, use it to evaluate
  
  if (!is.null(sdm_data$eval_data)) {
    eval_p <- sdm_data$eval_data$presence
    eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    eval_poly <- model.matrix(forms, data = eval_dat) 
    eval_poly <- eval_poly[,-1]
    
    pred_eval <- predict(m_final, eval_poly, type = "response")
    
    if (to_norm) {
      pred_eval <- .normalize_res(pred_eval)
    }
    
    eval_metric <- eval_metrics(eval_p, pred_eval[,1])
    
    result$eval_metrics <- eval_metric
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
  }
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
  
}





#' SDM module: Random Forest
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param method one of "classification" or "regression." If type is "down-sampled",
#' then method have to be "classification"
#' @param type either "normal" or "down-sampled" (see details)
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#'   
#' @details
#' Random forest is specially sensitive for class-imbalance (even for regression).
#' An alternative is to run a down-sampled classification, which will draw
#' equal sized samples of presence/background.
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_rf(species_data)
#' }
sdm_module_rf <- function(sdm_data, method = "classification", type = "down-sampled",
                          tune_blocks = "spatial_grid", blocks_all = FALSE) {
  
  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training
  
  if (method != "regression" | type == "down-sampled") {
    train_dat$presence <- as.factor(train_dat$presence)
  }
  
  if (type == "down-sampled") {
    
    method <- "classification"
    
    pres <- sum(train_p) 
    smpsize <- c("0" = pres, "1" = pres)
    
    ### Tune model
    
    ntree <- c(500, 750, 1000)
    
    tune_auc <- rep(NA, length(ntree))
    
    rf_fun <- function(labels, training, ntree) {
      pres <- sum(labels) 
      smpsize <- c("0" = pres, "1" = pres)
      randomForest::randomForest(formula = presence ~.,
                                 data = training,
                                 ntree = ntree, 
                                 sampsize = smpsize,
                                 replace = TRUE)
    }
    
    for (i in 1:length(ntree)) {
      m_train <- cv_mod(n_folds = sdm_data$blocks$n,
                        folds = sdm_data$blocks$folds[[tune_blocks]],
                        labels = train_p, training = train_dat,
                        prob_mode = T,
                        metrics = "auc", pred_type = "prob",
                        model_fun = rf_fun,
                        ntree = ntree[i])
      
      tune_auc[i] <- mean(unlist(m_train), na.rm = T)
    }
    
    timings <- .get_time(timings, "tuning")
    
    ### Evaluate final model through CV
    
    if (blocks_all) {
      final_cv_metric <- list()
      
      for (b in 1:length(sdm_data$blocks$folds)) {
        
        b_index <- sdm_data$blocks$folds[[b]]
        
        final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                       folds = b_index,
                                       labels = train_p, training = train_dat,
                                       prob_mode = T,
                                       metrics = "all", pred_type = "prob",
                                       model_fun = rf_fun,
                                       ntree = ntree[which.min(tune_auc)])
        
        final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
        
      }
      
      names(final_cv_metric) <- names(sdm_data$blocks$folds)
      cv_method <- names(sdm_data$blocks$folds)
      
    } else {
      
      final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                         folds = sdm_data$blocks$folds[[tune_blocks]],
                         labels = train_p, training = train_dat,
                         prob_mode = T,
                         metrics = "all", pred_type = "prob",
                         model_fun = rf_fun,
                         ntree = ntree[which.min(tune_auc)])
      
      final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
      
      cv_method <- tune_blocks
    }
    
    timings <- .get_time(timings, "cv")
    
    ### Get final model
    m_final <- randomForest::randomForest(formula = presence ~.,
                                          data = train_dat,
                                          ntree = ntree[which.min(tune_auc)], 
                                          sampsize = smpsize,
                                          replace = TRUE)
    
    ### Evaluate final model (full data)
    pred_final <- predict(m_final, train_dat, type = "prob")
    
    final_full_metric <- eval_metrics(train_p, pred_final[,2])
    
    timings <- .get_time(timings, "evaluate final")
    
    ### Prepare returning object
    result <- list(
      name = paste0("rf_", method, "_ds"),
      model = m_final,
      timings = timings,
      tune_cv_method = tune_blocks,
      cv_method = cv_method,
      parameters = c(
        n_trees = ntree[which.min(tune_auc)]
      ),
      cv_metrics = final_cv_metric,
      full_metrics = final_full_metric,
      eval_metrics = NULL
    )
    
    ### If testing dataset is available, use it to evaluate
    
    if (!is.null(sdm_data$eval_data)) {
      eval_p <- sdm_data$eval_data$presence
      eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
      
      pred_eval <- predict(m_final, eval_dat, type = "prob")
      
      eval_metric <- eval_metrics(eval_p, pred_eval[,2])
      
      result$eval_metrics <- eval_metric
      
      result$timings <- .get_time(result$timings, "evaluation dataset")
    }
    
    class(result) <- c("sdm_result", class(result))
    
    return(result)
    
    
  }
  
  if (type == "normal") {
    
    ### Tune model
    
    if (method == "classification") {
      prob <- TRUE
      predt <- "prob"
    } else {
      prob <- FALSE
      predt <- "response"
    }
    
    ntree <- c(500, 750, 1000)
    
    tune_auc <- rep(NA, length(ntree))
    
    rf_fun <- function(labels, training, ntree) {
      suppressWarnings(
        randomForest::randomForest(formula = presence ~.,
                                   data = training,
                                   ntree = ntree)
      )
    }
    
    for (i in 1:length(ntree)) {
      m_train <- cv_mod(n_folds = sdm_data$blocks$n,
                        folds = sdm_data$blocks$folds[[tune_blocks]],
                        labels = train_p, training = train_dat,
                        prob_mode = prob,
                        metrics = "auc", pred_type = predt,
                        model_fun = rf_fun,
                        ntree = ntree[i])
      
      tune_auc[i] <- mean(unlist(m_train), na.rm = T)
    }
    
    timings <- .get_time(timings, "tuning")
    
    ### Evaluate final model through CV
    
    if (blocks_all) {
      final_cv_metric <- list()
      
      for (b in 1:length(sdm_data$blocks$folds)) {
        
        b_index <- sdm_data$blocks$folds[[b]]
        
        final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                       folds = b_index,
                                       labels = train_p, training = train_dat,
                                       prob_mode = prob,
                                       metrics = "all", pred_type = predt,
                                       model_fun = rf_fun,
                                       ntree = ntree[which.min(tune_auc)])
        
        final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
        
      }
      
      names(final_cv_metric) <- names(sdm_data$blocks$folds)
      cv_method <- names(sdm_data$blocks$folds)
      
    } else {
      
      final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                         folds = sdm_data$blocks$folds[[tune_blocks]],
                         labels = train_p, training = train_dat,
                         prob_mode = prob,
                         metrics = "all", pred_type = predt,
                         model_fun = rf_fun,
                         ntree = ntree[which.min(tune_auc)])
      
      final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
      
      cv_method <- tune_blocks
    }
    
    timings <- .get_time(timings, "cv")
    
    ### Get final model
    m_final <- randomForest::randomForest(formula = presence ~.,
                                          data = train_dat,
                                          ntree = ntree[which.min(tune_auc)])
    
    ### Evaluate final model (full data)
    if (prob) {
      pred_final <- predict(m_final, train_dat, type = "prob")
      pred_final <- pred_final[,2]
    } else {
      pred_final <- predict(m_final, train_dat, type = "response")
    }
    
    final_full_metric <- eval_metrics(train_p, pred_final)
    
    timings <- .get_time(timings, "evaluate final")
    
    ### Prepare returning object
    result <- list(
      name = paste0("rf_", method, "_normal"),
      model = m_final,
      timings = timings,
      tune_cv_method = tune_blocks,
      cv_method = cv_method,
      parameters = c(
        n_trees = ntree[which.min(tune_auc)]
      ),
      cv_metrics = final_cv_metric,
      full_metrics = final_full_metric,
      eval_metrics = NULL
    )
    
    ### If testing dataset is available, use it to evaluate
    
    if (!is.null(sdm_data$eval_data)) {
      eval_p <- sdm_data$eval_data$presence
      eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
      
      if (prob) {
        pred_eval <- predict(m_final, eval_dat, type = "prob")
        pred_eval <- pred_eval[,2]
      } else {
        pred_eval <- predict(m_final, eval_dat, type = "response")
      }
      
      eval_metric <- eval_metrics(eval_p, pred_eval)
      
      result$eval_metrics <- eval_metric
      
      result$timings <- .get_time(result$timings, "evaluation dataset")
    }
    
    class(result) <- c("sdm_result", class(result))
    
    return(result)
    
    
  }
}



#' SDM module: Boosted Regression Trees
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param weight_resp if weight_resp is \code{TRUE},
#'  then the response is weighted in a way that the total weight of presences
#'  is the same as the weight of background points (naive weighting)
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#'
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_brt(species_data)
#' }
sdm_module_brt <- function(sdm_data, weight_resp = TRUE,
                           tune_blocks = "spatial_grid", blocks_all = FALSE) {
  
  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training
  
  b_index <- sdm_data$blocks$folds[[tune_blocks]]
  
  if (weight_resp) {
    pres <- length(train_p[train_p == 1])
    bkg <- length(train_p[train_p == 0])
    wt <- ifelse(train_p == 1, 1, pres / bkg)
  } else {
    wt <- rep(1, length(train_p))
  }
  
  ### Tune model
  
  m_train <- dismo::gbm.step(data = train_dat,
                             gbm.x = 2:ncol(train_dat),
                             gbm.y = 1,
                             fold.vector = b_index,
                             family = "bernoulli", 
                             tree.complexity = 1, #5
                             learning.rate = 0.01, #0.001
                             step.size = 100, #50
                             bag.fraction = 0.75,
                             max.trees = 10000,
                             n.trees = 100, #50
                             n.folds = 5,
                             site.weights = wt,
                             plot.main = F,
                             verbose = F,
                             silent = T)
  
  timings <- .get_time(timings, "tuning")
  
  ### Evaluate final model through CV
  
  gbm_fmod <- function(labels, training, weights, ...) {
    fdata <- cbind(labels, training)
    dismo::gbm.fixed(data = fdata, gbm.x = 2:ncol(fdata), gbm.y = 1, verbose = F,
              site.weights = weights, ...)
  }
  
  if (blocks_all) {
    final_cv_metric <- list()
    
    for (b in 1:length(sdm_data$blocks$folds)) {
      
      b_index <- sdm_data$blocks$folds[[b]]
      
      final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                     folds = b_index,
                                     labels = train_p, training = train_dat[,-1],
                                     weights = wt, metrics = "all", pred_type = "response",
                                     model_fun = gbm_fmod,
                                     tree.complexity = m_train$interaction.depth,
                                     learning.rate = m_train$shrinkage,
                                     n.trees = m_train$n.trees,
                                     bag.fraction = m_train$bag.fraction)
      
      final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
      
    }
    
    names(final_cv_metric) <- names(sdm_data$blocks$folds)
    cv_method <- names(sdm_data$blocks$folds)
    
  } else {
    
    final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                       folds = sdm_data$blocks$folds[[tune_blocks]],
                       labels = train_p, training = train_dat[,-1],
                       weights = wt, metrics = "all", pred_type = "response",
                       model_fun = gbm_fmod,
                       tree.complexity = m_train$interaction.depth,
                       learning.rate = m_train$shrinkage,
                       n.trees = m_train$n.trees,
                       bag.fraction = m_train$bag.fraction)
    
    final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
    
    cv_method <- tune_blocks
  }
  
  timings <- .get_time(timings, "cv")
  
  ### Get final model
  m_final <- m_train
  
  ### Evaluate final model (full data)
  
  pred_final <- predict(m_final, train_dat[,-1], type = "response")
  
  final_full_metric <- eval_metrics(train_p, pred_final)
  
  timings <- .get_time(timings, "evaluate final")
  
  ### Prepare returning object
  result <- list(
    name = "brt",
    model = m_final,
    timings = timings,
    tune_cv_method = tune_blocks,
    cv_method = cv_method,
    parameters = c(
      tree_complexity = m_train$interaction.depth,
      learning_rate = m_train$shrinkage,
      n_trees = m_train$n.trees,
      bag_fraction = m_train$bag.fraction
    ),
    cv_metrics = final_cv_metric,
    full_metrics = final_full_metric,
    eval_metrics = NULL
  )
  
  ### If testing dataset is available, use it to evaluate
  
  if (!is.null(sdm_data$eval_data)) {
    eval_p <- sdm_data$eval_data$presence
    eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    pred_eval <- predict(m_final, eval_dat, type = "response")
    
    eval_metric <- eval_metrics(eval_p, pred_eval)
    
    result$eval_metrics <- eval_metric
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
  }
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
    
}


#' SDM module: MAXENT (Maximum Entropy)
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param method either "maxnet" to use the [maxnet::maxnet()] implementation or
#' "maxent" to use the Java implementation (see note).
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#' @param maxent_pred_layers optional. If method is "maxent", there is no predict
#' method available. In this case, if you want the predictions, you need to supply 
#' the prediction layers (in SpatRaster format) here
#' 
#' @note
#' The dismo implementation of the Maxent (through rJava) is causing R to crash.
#' Thus we use here a workaround which call the java program through a system call.
#' It's advisable to use the maxnet implementation which produces basically identical
#' outputs and have a predict method.
#'
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_maxent(species_data)
#' }
sdm_module_maxent <- function(sdm_data, method = "maxnet",
                              tune_blocks = "spatial_grid", blocks_all = FALSE,
                              maxent_pred_layers = NULL) {
  
  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  ### Tune model
  
  tune_grid <- expand.grid(
    remult = seq(0.5, 2, 0.5),
    features = c("lq", "lqh"),
    stringsAsFactors = F
  )
  
  tune_grid <- unique(tune_grid)
  
  auc_test <- rep(NA, nrow(tune_grid))
  
  if (method == "maxnet") {
    
    # Create a function of the model
    maxnet_mod <- function(labels, training, features, rm) {
      maxnet::maxnet(p = labels,
                     data = training,
                     f = maxnet::maxnet.formula(p = labels,
                                                data = training,
                                                classes = features),
                     regmult = rm,
                     addsamplestobackground = F)
    }
    
    for (z in 1:nrow(tune_grid)) {
     
      b_index <- sdm_data$blocks$folds[[tune_blocks]]
      
      auc_block <- cv_mod(n_folds = sdm_data$blocks$n,
                          folds = b_index, labels = train_p,
                          training = train_dat, metrics = "auc",
                          pred_type = "cloglog", model_fun = maxnet_mod,
                          features = tune_grid$features[z], rm = tune_grid$remult[z])
      
      auc_test[z] <- mean(unlist(auc_block), na.rm = T)
    }
    
    timings <- .get_time(timings, "tuning")
    
    ### Evaluate final model through CV
    
    best_tune <- tune_grid[which.min(auc_test),]
    
    if (blocks_all) {
      final_cv_metric <- list()
      
      for (b in 1:length(sdm_data$blocks$folds)) {
        
        b_index <- sdm_data$blocks$folds[[b]]
        
        final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                       folds = b_index, labels = train_p,
                                       training = train_dat, metrics = "all",
                                       pred_type = "cloglog", model_fun = maxnet_mod,
                                       features = best_tune$features, rm = best_tune$remult)
        
        final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
        
      }
      
      names(final_cv_metric) <- names(sdm_data$blocks$folds)
      cv_method <- names(sdm_data$blocks$folds)
      
    } else {
      
      final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                         folds = b_index, labels = train_p,
                         training = train_dat, metrics = "all",
                         pred_type = "cloglog", model_fun = maxnet_mod,
                         features = best_tune$features, rm = best_tune$remult)
      
      final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
      
      cv_method <- tune_blocks
    }
    
    timings <- .get_time(timings, "cv")
    
    ### Train final model
    m_final <- maxnet_mod(labels = train_p, training = train_dat,
                          features = best_tune$features, rm = best_tune$remult)
    
    ### Evaluate final model (full data)
    
    pred_final <- predict(m_final, train_dat, type = "cloglog")
    
    final_full_metric <- eval_metrics(train_p, pred_final)
    
    timings <- .get_time(timings, "evaluate final")
    
    ### Prepare returning object
    result <- list(
      name = method,
      model = m_final,
      timings = timings,
      tune_cv_method = tune_blocks,
      cv_method = cv_method,
      parameters = best_tune,
      cv_metrics = final_cv_metric,
      full_metrics = final_full_metric,
      eval_metrics = NULL
    )
    
    ### If testing dataset is available, use it to evaluate
    
    if (!is.null(sdm_data$eval_data)) {
      eval_p <- sdm_data$eval_data$presence
      eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
      
      pred_eval <- predict(m_final, eval_dat, type = "cloglog")
      
      eval_metric <- eval_metrics(eval_p, pred_eval)
      
      result$eval_metrics <- eval_metric
      
      result$timings <- .get_time(result$timings, "evaluation dataset")
    }
    
    class(result) <- c("sdm_result", class(result))
    
    return(result)
    
  }
  
  if (method == "maxent") {
    
    ncords <- sdm_data$coord_training
    colnames(ncords) <- c("longitude", "latitude")
    
    train_dat <- cbind(ncords, train_dat)
    
    # rJava and dismo crashing on R version 4.3.2 with RStudio
    # maxmod <- dismo::maxent(x = train_dat,
    #                         p = train_p,
    #                         removeDuplicates = FALSE,
    #                         path = "output/maxent_files",
    #                         args = c("betamultiplier=1", "noautofeature", "nothreshold", "nohinge", "noproduct"))
    # Attempting alternative approach
    
    # Create a function of the model
    maxent_mod <- function(labels, training, test, features, rm, return_pred = NULL) {
      
      tdir <- tempdir()
      
      if (dir.exists(paste0(tdir, "/output"))) {
        unlink(paste0(tdir, "/output"), recursive = T)
      }
      
      dir.create(paste0(tdir, "/output"))
      
      train <- cbind(species = "species", training)
      
      back <- train[labels == 0,]
      train <- train[labels == 1,]
      
      test <- cbind(species = "species", test)
      
      write.csv(train, paste0(tdir, "/species_swd.csv"), row.names = F)
      write.csv(back, paste0(tdir, "/back_swd.csv"), row.names = F)
      write.csv(test, paste0(tdir, "/test_swd.csv"), row.names = F)
      
      max_path <- find.package("dismo")
      max_path <- paste0(max_path, "/java/maxent.jar")
      
      ft <- "noautofeature nothreshold noproduct"
      
      if (!grepl("q", features)) {
        ft <- paste(ft, "noquadratic")
      }
      if (!grepl("h", features)) {
        ft <- paste(ft, "nohinge")
      }
      
      if (!is.null(return_pred)) {
        if (dir.exists(paste0(tdir, "/proj"))) {
          unlink(paste0(tdir, "/proj"), recursive = T)
        }
        
        dir.create(paste0(tdir, "/proj"))
        
        writeRaster(return_pred, paste0(tdir, "/proj/", names(return_pred), ".asc"),
                    NAflag=-128, overwrite = T)
        
        pred_ob <- paste0("projectionlayers=", tdir, "/proj")
      } else {
        pred_ob <- paste0("projectionlayers=", tdir, "/test_swd.csv")
      }
      
      system(paste("java -mx512m -jar", max_path,
                   paste(
                     paste0("samplesfile=", tdir, "/species_swd.csv"),
                     paste0("environmentallayers=", tdir, "/back_swd.csv"),
                     pred_ob,
                     paste0("outputdirectory=", tdir, "/output"),
                     paste0("betamultiplier=", rm),
                     ft, "redoifexists autorun visible=False"
                   )))
      
      if (is.null(return_pred)) {
        pred <- read.csv(paste0(tdir, "/output/species_test_swd.csv"))
        
        return(pred[,3])
      } else {
        pred <- rast(paste0(tdir, "/output/species_proj.asc"))
        
        return(pred)
      }
    }
    
    cv_mod_max <- function(n_folds, folds, labels, training, metrics, model_fun, ...) {
      
      if (metrics == "all") {
        metrics <- c("auc", "tss", "prg", "cbi")
      }
      
      cv_metrics <- lapply(1:n_folds, function(index){
        m_block <- model_fun(labels[folds != index], training[folds != index,],
                             test = training[folds == index,],
                             ...)
        
        em <- eval_metrics(labels[folds == index], as.vector(m_block), metrics)
        
        return(em)
      })
      
      return(cv_metrics)
    }
    
    
    for (z in 1:nrow(tune_grid)) {
      
      b_index <- sdm_data$blocks$folds[[tune_blocks]]
      
      auc_block <- cv_mod_max(n_folds = sdm_data$blocks$n,
                          folds = b_index, labels = train_p,
                          training = train_dat, metrics = "auc",
                          model_fun = maxent_mod,
                          features = tune_grid$features[z], rm = tune_grid$remult[z])
      
      auc_test[z] <- mean(unlist(auc_block), na.rm = T)
    }
    
    timings <- .get_time(timings, "tuning")
    
    ### Evaluate final model through CV
    
    best_tune <- tune_grid[which.min(auc_test),]
    
    if (blocks_all) {
      final_cv_metric <- list()
      
      for (b in 1:length(sdm_data$blocks$folds)) {
        
        b_index <- sdm_data$blocks$folds[[b]]
        
        final_cv_metric[[b]] <- cv_mod_max(n_folds = sdm_data$blocks$n,
                                           folds = b_index, labels = train_p,
                                           training = train_dat, metrics = "all",
                                           model_fun = maxent_mod,
                                           features = best_tune$features, rm = best_tune$remult)
        
        final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
        
      }
      
      names(final_cv_metric) <- names(sdm_data$blocks$folds)
      cv_method <- names(sdm_data$blocks$folds)
      
    } else {
      
      final_cv <- cv_mod_max(n_folds = sdm_data$blocks$n,
                             folds = b_index, labels = train_p,
                             training = train_dat, metrics = "all",
                             model_fun = maxent_mod,
                             features = best_tune$features, rm = best_tune$remult)
      
      final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
      
      cv_method <- tune_blocks
    }
    
    timings <- .get_time(timings, "cv")
    
    ### Train final model
    m_final <- maxent_mod(labels = train_p, training = train_dat, test = NULL,
                          features = best_tune$features, rm = best_tune$remult,
                          return_pred = maxent_pred_layers)
    
    ### Evaluate final model (full data)
    
    pred_final <- terra::extract(m_final, sdm_data$coord_training)
    
    final_full_metric <- eval_metrics(train_p, pred_final[,2])
    
    timings <- .get_time(timings, "evaluate final")
    
    ### Prepare returning object
    result <- list(
      name = method,
      model = m_final,
      timings = timings,
      tune_cv_method = tune_blocks,
      cv_method = cv_method,
      parameters = best_tune,
      cv_metrics = final_cv_metric,
      full_metrics = final_full_metric,
      eval_metrics = NULL
    )
    
    ### If testing dataset is available, use it to evaluate
    
    # if (!is.null(sdm_data$eval_data)) {
    #   eval_p <- sdm_data$eval_data$presence
    #   eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    #   
    #   pred_eval <- predict(m_final, eval_dat, type = "cloglog")
    #   
    #   eval_metric <- eval_metrics(eval_p, pred_eval)
    #   
    #   result$eval_metrics <- eval_metric
    #   
    #   result$timings <- .get_time(result$timings, "evaluation dataset")
    # }
    
    class(result) <- c("sdm_result", class(result))
    
    return(result)
    
  }
  
}



#' SDM module: LightGBM (Light Gradient Boosting Model)
#'
#' @param sdm_data object of class sdm_dat returned by [mp_prepare_data()]
#' @param weight_resp if weight_resp is \code{TRUE},
#'  then the response is weighted in a way that the total weight of presences
#'  is the same as the weight of background points (naive weighting)
#' @param tune_blocks the name of the blocks (added with [mp_prepare_blocks()])
#'   which should be used for tuning the models
#' @param blocks_all if \code{TRUE}, then after tuning the model, it will
#'   produce cross-validated metrics for all types of blocks available (e.g.
#'   latitudinal)
#'   
#' @details
#' LightGBM is a powerfull tool for machine learning. It was never applied to SDMs,
#' but its robustness and speed makes it a good candidate for distribution modelling.
#' We implement a module here mainly for testing purposes.
#' 
#'
#'
#' @return an sdm_result object (a list) containing:
#' - name = name of method
#' - model = the model object
#' - timings = a character vector of the timings for each step
#' - parameters = the best tuned parameters
#' - cv_metrics = the cross-validated metrics
#' - full_metrics = the metrics of the final model fitted using the whole data.
#'  Just for reference, as this is not a fair scenario
#'  (same data used for training and testing).
#' - eval_metrics = metrics for the evaluation dataset (if this is supplied).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sdm_module_brt(species_data)
#' }
sdm_module_lgbm <- function(sdm_data, weight_resp = TRUE,
                            tune_blocks = "spatial_grid", blocks_all = FALSE) {

  .check_type(sdm_data)
  
  timings <- .get_time()
  
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  train_dat <- as.matrix(train_dat)
  
  if (weight_resp) {
    pres <- length(train_p[train_p == 1])
    bkg <- length(train_p[train_p == 0])
    wt <- ifelse(train_p == 1, 1, pres / bkg)
  } else {
    wt <- NULL
  }
  
  ### Tune model
  lgbm_mod <- function(labels, training, weights = NULL, max_dep, lr) {
    
    to_hold <- sample(1:nrow(training), nrow(training)*0.2)
    
    to_include <- 1:nrow(training)
    to_include <- to_include[-to_hold]
    
    lightgbm::lgb.train(
      params = list(
        objective = "regression",
        metric = "auc",
        max_depth = max_dep,
        num_iterations = 1000,
        early_stopping_rounds = 100,
        learning_rate = lr,
        bagging_fraction = 0.8,
        bagging_freq = 1,
        feature_fraction = 0.8
        #num_leaves = nl
        #is_unbalance = T
        #feature_fraction = .8
      ),
      valids = list(test = lightgbm::lgb.Dataset(training[to_hold,], label = labels[to_hold], weight = weights[to_hold])),
      data = lightgbm::lgb.Dataset(training[to_include,], label = labels[to_include], weight = weights[to_include])
    )
  }
  
  tune_grid <- expand.grid(
    max_dep = seq(2, 4, 1), # Deeper, more complex models
    lr = c(0.1, 0.05, 0.01, 0.005, 0.001) # controls the step size, lower learns "more"
    #, num_leaves = c(10, 30) # controls the maximum number of leaves. High, more complex but risky overfit
    # , feature_fraction = c(0.8) # controls the fraction of features to consider for each tree, can reduce overfit
    # , bagging_fraction = c(0.8) # fraction of data to be used for each iteration, reduces overfit
    # , early_stopping_round = 20 # number of rounds after which the training will stop if there's no improvement in the validation metric
  )
  
  tune_grid <- unique(tune_grid)
  
  auc_test <- rep(NA, nrow(tune_grid))
  
  for (z in 1:nrow(tune_grid)) {
    
    b_index <- sdm_data$blocks$folds[[tune_blocks]]
    
    auc_block <- cv_mod(n_folds = sdm_data$blocks$n,
                        folds = b_index, labels = train_p,
                        training = train_dat, metrics = "auc",
                        weights = wt, max_dep = tune_grid$max_dep[z],
                        lr = tune_grid$lr[z], 
                        pred_type = NULL, model_fun = lgbm_mod)
    
    auc_test[z] <- mean(unlist(auc_block), na.rm = T)
  }
  
  timings <- .get_time(timings, "tuning")
  
  
  ### Evaluate final model through CV
  
  if (blocks_all) {
    final_cv_metric <- list()
    
    for (b in 1:length(sdm_data$blocks$folds)) {
      
      b_index <- sdm_data$blocks$folds[[b]]
      
      final_cv_metric[[b]] <- cv_mod(n_folds = sdm_data$blocks$n,
                                     folds = b_index, labels = train_p,
                                     training = train_dat, metrics = "all",
                                     weights = wt,
                                     max_dep = tune_grid$max_dep[which.min(auc_test)],
                                     lr = tune_grid$lr[which.min(auc_test)],
                                     pred_type = NULL, model_fun = lgbm_mod)
      
      final_cv_metric[[b]] <- as.data.frame(do.call("rbind", final_cv_metric[[b]]))
      
    }
    
    names(final_cv_metric) <- names(sdm_data$blocks$folds)
    cv_method <- names(sdm_data$blocks$folds)
    
  } else {
    
    final_cv <- cv_mod(n_folds = sdm_data$blocks$n,
                       folds = b_index, labels = train_p,
                       training = train_dat, metrics = "all",
                       weights = wt,
                       max_dep = tune_grid$max_dep[which.min(auc_test)],
                       lr = tune_grid$lr[which.min(auc_test)],
                       pred_type = NULL, model_fun = lgbm_mod)
    
    final_cv_metric <- as.data.frame(do.call("rbind", final_cv))
    
    cv_method <- tune_blocks
  }
  
  timings <- .get_time(timings, "cv")
  
  ### Train final model
  # m_final <- lgb.train(
  #   params = list(
  #     objective = "binary",
  #     metric = "auc",
  #     max_depth = tune_grid$max_dep[which.min(auc_test)],
  #     num_iterations = 200,
  #     boosting_rounds = 50,
  #     #early_stopping_rounds = 40,
  #     learning_rate = tune_grid$lr[which.min(auc_test)],
  #     is_unbalance = T
  #     #feature_fraction = .8
  #   ),
  #   data = lgb.Dataset(as.matrix(train_dat), label = train_p, weight = wt)
  # )
  
  m_final <- lgbm_mod(train_p, train_dat, wt,
                      tune_grid$max_dep[which.min(auc_test)],
                      tune_grid$lr[which.min(auc_test)])
  
  ### Evaluate final model (full data)
  
  pred_final <- predict(m_final, as.matrix(train_dat))
  
  final_full_metric <- eval_metrics(train_p, pred_final)
  
  timings <- .get_time(timings, "evaluate final")
  
  ### Prepare returning object
  result <- list(
    name = "lgbm",
    model = m_final,
    timings = timings,
    tune_cv_method = tune_blocks,
    cv_method = cv_method,
    parameters = tune_grid[which.min(auc_test),],
    cv_metrics = final_cv_metric,
    full_metrics = final_full_metric,
    eval_metrics = NULL
  )
  
  ### If testing dataset is available, use it to evaluate
  
  if (!is.null(sdm_data$eval_data)) {
    eval_p <- sdm_data$eval_data$presence
    eval_dat <- sdm_data$eval_data[, !colnames(sdm_data$eval_data) %in% "presence"]
    
    pred_eval <- predict(m_final, as.matrix(eval_dat))
    
    eval_metric <- eval_metrics(eval_p, pred_eval)
    
    result$eval_metrics <- eval_metric
    
    result$timings <- .get_time(result$timings, "evaluation dataset")
  }
  
  class(result) <- c("sdm_result", class(result))
  
  return(result)
}




# Internal functions

# Check if object is of type sdm_dat
#' @export
.check_type <- function(x) {
  if (class(x)[1] != "sdm_dat") {
    cli::cli_abort("sdm_data should be of type sdm_dat (generated with mp_prepare_data)")
  } else {
    return(invisible(NULL))
  }
}

# Get evaluation times
#' @export
.get_time <- function(previous = NULL, name = NULL) {
  if (is.null(previous)) {
    the_time <- c(NA)
    attributes(the_time) <- list(start = Sys.time())
    # the_time <- c(Sys.time())
    # names(the_time) <- "start"
  } else {
    n <- difftime(Sys.time(), attr(previous, "start"), units='mins')
    names(n) <- name
    if (is.na(previous[1])) {
      the_time <- n
      attr(the_time, "start") <- attr(previous, "start")
    } else{
      the_time <- c(previous, n)
      attr(the_time, "start") <- attr(previous, "start")
    }
  }
  return(the_time)
}

# Extract maximum time of sdm_result object
#' @export
.extract_time <- function(model) {
  max(model$timings)
}

# Normalize result (0-1)
#' @export
.normalize_res <- function(pred_vals) {
  (pred_vals - min(pred_vals)) / (max(pred_vals) - min(pred_vals))
}
