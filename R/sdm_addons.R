############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
######################## SDM modules accessory functions #######################

#' Cross-validate model
#'
#' @param n_folds the number of folds
#' @param folds folds index
#' @param labels the presence-absence vector
#' @param training training data
#' @param metrics which metrics return
#' @param normalize rescale the predicted values to the scale 0-1. Relevant for
#'   down-weighted or infinitely weighted models.
#' @param pred_type type of prediction
#' @param model_fun the function of the model
#' @param return_block_pred if \code{TRUE}, then the predictions for each block
#'   are returned.
#' @param ... additional parameters for the model function
#'
#' @return a list with the metrics for each fold
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
cv_mod <- function(n_folds, folds, labels, training, weights = NULL, 
                   prob_mode = F, metrics, normalize = FALSE,
                   pred_type, model_fun, return_block_pred = FALSE, ...) {
  
  if (metrics == "all") {
    metrics <- c("auc", "cbi", "pr", "prg", "tss", "spec", "sens", "kap", "fmeas", "opr", "upr")
  }
  
  if (!is.null(weights)) {
    cv_metrics <- lapply(1:n_folds, function(index){
      m_block <- model_fun(labels[folds != index], training[folds != index,],
                           weights = weights[folds != index], ...)
      
      pred_test <- predict(m_block, training[folds == index,], type = pred_type)
      
      if (normalize) {
        pred_test <- .normalize_res(pred_test)
      }
      
      em <- eval_metrics(labels[folds == index], as.vector(pred_test), metrics)
      
      if (return_block_pred) {
        em <- list(block_pred = pred_test,
                   metric = em)
      }
      
      return(em)
    })
  } else {
    cv_metrics <- lapply(1:n_folds, function(index){
      m_block <- model_fun(labels[folds != index], training[folds != index,],
                           ...)
      
      pred_test <- predict(m_block, training[folds == index,], type = pred_type)
      
      if (normalize) {
        pred_test <- .normalize_res(pred_test)
      }
      
      if (prob_mode) {
        em <- eval_metrics(labels[folds == index], as.vector(pred_test[,2]), metrics)
      } else {
        em <- eval_metrics(labels[folds == index], as.vector(pred_test), metrics)
      }
      
      if (return_block_pred) {
        em <- list(block_pred = pred_test,
                   metric = em)
      }
      
      return(em)
    })
  }
  
  return(cv_metrics)
}



#' Evaluate models with common metrics
#'
#' @param original a vector with the original values
#' @param predicted a vector with predicted values
#' @param metrics a character vector with one or more metrics (see details)
#' @param thresholds a character vector with one or more thresholds to be used
#'   for those metrics that need thresholding (TSS, Specificity, Sensitivity and
#'   Kappa; see details). Accepted values are "p10", "maxsss" and "mtp".
#' 
#' @details
#' Eight metrics are available:
#' - \code{auc}: Area Under the Curve
#' - \code{cbi}: Continuous Boyce Index
#' - \code{pr}: Area under the Precision-Recall curve
#' - \code{prg}: Precision-Recall-Gain curve
#' - \code{tss}: True Skill Statistics
#' - \code{spec}: Specificity (True Negative Rate)
#' - \code{sens}: Sensitivity (True Positive Rate)
#' - \code{kap}: Kappa statistics
#' - \code{fmeas}: F-measure (Sorensen's similarity index)
#' - \code{opr}: Overprediction rate
#' - \code{upr}: (underprediction rate)
#' 
#' For the thresholding methods, three options are available:
#' - \code{p10}: The value at which 90% of presence points are included
#' - \code{mtp}: The minimum training presence, i.e., the value at which all 
#'   training presences are included.
#' - \code{maxsss}: The value that maximizes the sum of specificity plus sensitivity.
#' 
#'
#' @return a data.frame with the chosen metrics
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
eval_metrics <- function(original, predicted,
                         metrics = c("auc", "cbi", "pr", "prg", "tss", "spec", "sens", "kap", "fmeas", "opr", "upr"),
                         thresholds = c("maxsss", "mtp", "p10")) {
  
  th_based <- c("tss", "spec", "sens", "kap", "fmeas", "opr", "upr")
  
  if (any(th_based %in% metrics)) {
    no_thresh <- metrics[!metrics %in% th_based]
    thresh <- metrics[metrics %in% th_based]
    
    thresholds <- rep(thresholds, each = length(thresh))
    
    thresh <- paste(thresh, thresholds, sep = "_")
    
    metrics <- c(no_thresh, thresh)
    
    thresh_res <- list()
    
    thresh_names <- unique(thresholds)
    
    for (z in 1:length(thresh_names)) {
      
      if (thresh_names[z] == "maxsss") {
        to_get <- "maxSSS"
        q <- 0
      }
      if (thresh_names[z] == "mtp") {
        to_get <- "MTP"
        q <- 0
      }
      if (thresh_names[z] == "p10") {
        to_get <- "MTP"
        q <- 0.1
      }
      
      th_val <- modEvA::getThreshold(obs = original,
                                     pred = predicted,
                                     threshMethod = to_get,
                                     quant = q)
      
      suppressWarnings(th_res <- modEvA::threshMeasures(obs = original,
                                                        pred = predicted,
                                                        thresh = th_val,
                                                        plot = F, standardize = F))
      
      thresh_res[[z]] <- th_res$ThreshMeasures
      
    }
    
    names(thresh_res) <- thresh_names
  }
  
  if (any(c("auc", "pr") %in% metrics)) {
    auc_m <- precrec::auc(precrec::evalmod(scores = predicted, 
                                           labels = original))
  }
  
  res <- rep(NA, length(metrics))
  metric_name <- gsub("_[^_]*", "", metrics)
  thresh_type <- gsub("([^_]+)_", "", metrics)
  
  # To improve: measures extraction
  for (i in 1:length(metrics)) {
    
    res[i] <- switch(metric_name[i],
                     auc = auc_m[auc_m$curvetypes == "ROC", 4],
                     cbi = modEvA::Boyce(obs = original, pred = predicted, plot = F)$Boyce,
                     pr = auc_m[auc_m$curvetypes == "PRC", 4],
                     prg = prg::calc_auprg(prg::create_prg_curve(labels = original, pos_scores = predicted)),
                     tss = thresh_res[[thresh_type[i]]]["TSS",],
                     spec = thresh_res[[thresh_type[i]]]["Specificity",],
                     sens = thresh_res[[thresh_type[i]]]["Sensitivity",],
                     kap = thresh_res[[thresh_type[i]]]["kappa",],
                     fmeas = thresh_res[[thresh_type[i]]]["F1score",],
                     opr = thresh_res[[thresh_type[i]]]["OPR",],
                     upr = thresh_res[[thresh_type[i]]]["UPR",]
                     )
  }
  
  names(res) <- metrics
  
  return(res)
}



# Print the SDM result returned by the SDM modules
#' @export
print.sdm_result <- function(x, print_all = FALSE) {
  
  cli::cli_h1(x$name)
  
  cli::cli_alert_info("Time (in minutes, from start): {paste(names(x$timings), '=', round(x$timings, 2))}")
  cli::cli_alert_info("Best tune parameters: {paste(names(x$parameters), '=', x$parameters)}")
  cli::cli_alert_info("Tuned using [{x$tune_cv_method}] CV")
  
  cli::cli_h3("Metrics for the full model")
  if (length(x$full_metrics) > 3 & !print_all) {
    print_note <- TRUE
    x$full_metrics <- x$full_metrics[1:3]
    if (!is.data.frame(x$cv_metrics)) {
      for (i in 1:length(x$cv_metrics)) {
        x$cv_metrics[[i]] <- x$cv_metrics[[i]][,1:3]
      }
    } else {
      x$cv_metrics <- x$cv_metrics[,1:3]
    }
  }
  m <- paste(names(x$full_metrics), "=", round(x$full_metrics, 2))
  names(m) <- rep(">", length(names(m)))
  
  cli::cli_bullets(m)
  
  cli::cli_h3("Metrics for the CV model")
  
  if (!is.data.frame(x$cv_metrics)) {
    for (i in 1:length(x$cv_metrics)) {
      cli::cat_line("Cross-validation [", x$cv_method[i], "]")
      cvm <- apply(x$cv_metrics[[i]], 2, mean)
      cvsd <- apply(x$cv_metrics[[i]], 2, sd)
      
      m <- paste(names(cvm), "=", round(cvm, 2), "±", round(cvsd, 2))
      names(m) <- rep(">", length(names(m)))
      
      cli::cli_bullets(m)
      cli::cat_line()
    }
  } else {
    cli::cat_line("Cross-validated with [", x$cv_method, "]")
    
    cvm <- apply(x$cv_metrics, 2, mean)
    cvsd <- apply(x$cv_metrics, 2, sd)
    
    m <- paste(names(cvm), "=", round(cvm, 2), "±", round(cvsd, 2))
    names(m) <- rep(">", length(names(m)))
    
    cli::cli_bullets(m)
  }
  if (print_note) {
    cli::cli_inform("{.emph More than 3 metrics available, only the first 3 were printed. Set argument {.code print_all = TRUE} to print all.}")
  }
  
  return(invisible(NULL))
}



# Predict methods for the models returned by the SDM modules
#' @export
setOldClass("sdm_result")

#' @export
setMethod("predict", signature(object = "sdm_result"),
          function(object, layers, ...) {
            x <- object
            if (class(layers)[1] == "SpatRaster") {
              switch (x$name,
                      maxnet = predict(layers, x$model, type = "cloglog", na.rm = T, ...),
                      brt = predict(layers, x$model, type = "response", na.rm = T, ...),
                      rf_classification_ds = {
                        p <- predict(layers, x$model, type = "prob")
                        p[[2]]
                      },
                      rf_classification_normal = {
                        p <- predict(layers, x$model, type = "prob")
                        p[[2]]
                      },
                      rf_regression_normal = predict(layers, x$model, type = "response", ...),
                      gam_iwlr = predict(layers, x$model, type = "response", ...),
                      gam_dwpr = predict(layers, x$model, type = "response", ...),
                      gam_normal = predict(layers, x$model, type = "response", ...),
                      glm_iwlr = predict(layers, x$model, type = "response", ...),
                      glm_dwpr = predict(layers, x$model, type = "response", ...),
                      glm_normal = predict(layers, x$model, type = "response", ...),
                      lasso = {
                        lay_vals <- values(layers)
                        forms <- as.formula(paste("~ 1",
                                                  paste(colnames(lay_vals), collapse = "+"),
                                                  paste(paste("I(", colnames(lay_vals), "^2", ")", sep = ""), 
                                                        collapse = "+"), sep = "+"))
                        layers_poly <- model.matrix(forms,
                                                    model.frame(forms, as.data.frame(lay_vals),
                                                                na.action=na.pass)) 
                        layers_poly <- layers_poly[,-1]
                        p <- predict(x$model, as.matrix(layers_poly), s = "lambda.1se", type = "response", ...)
                        f <- layers[[1]]
                        values(f) <- p
                        f
                      },
                      elasticnet = {
                        lay_vals <- values(layers)
                        forms <- as.formula(paste("~ 1",
                                                  paste(colnames(lay_vals), collapse = "+"),
                                                  paste(paste("I(", colnames(lay_vals), "^2", ")", sep = ""), 
                                                        collapse = "+"), sep = "+"))
                        layers_poly <- model.matrix(forms,
                                                    model.frame(forms, as.data.frame(lay_vals),
                                                                na.action=na.pass)) 
                        layers_poly <- layers_poly[,-1]
                        p <- predict(x$model, as.matrix(layers_poly), s = "lambda.1se", type = "response", ...)
                        f <- layers[[1]]
                        values(f) <- p
                        f
                      },
                      ridge = {
                        lay_vals <- values(layers)
                        forms <- as.formula(paste("~ 1",
                                                  paste(colnames(lay_vals), collapse = "+"),
                                                  paste(paste("I(", colnames(lay_vals), "^2", ")", sep = ""), 
                                                        collapse = "+"), sep = "+"))
                        layers_poly <- model.matrix(forms,
                                                    model.frame(forms, as.data.frame(lay_vals),
                                                                na.action=na.pass)) 
                        layers_poly <- layers_poly[,-1]
                        p <- predict(x$model, as.matrix(layers_poly), s = "lambda.1se", type = "response", ...)
                        f <- layers[[1]]
                        values(f) <- p
                        f
                      },
                      lgbm = {
                        vals <- values(layers)
                        pred <- predict(x$model, as.matrix(vals))
                        p <- layers[[1]]
                        values(p) <- pred
                        p <- mask(p, layers[[1]])
                        p
                      }
              )
            } else {
              switch (x$name,
                      maxnet = predict(x$model, layers, type = "cloglog", na.rm = T, ...),
                      brt = predict(x$model, layers, type = "response", na.rm = T, ...),
                      rf_classification_ds = {
                        p <- predict(x$model, layers, type = "prob")
                        p[,2]
                      },
                      rf_classification_normal = {
                        p <- predict(x$model, layers, type = "prob")
                        p[,2]
                      },
                      rf_regression_normal = predict(x$model, layers, type = "response", ...),
                      gam_iwlr = predict(x$model, layers, type = "response", ...),
                      gam_dwpr = predict(x$model, layers, type = "response", ...),
                      gam_normal = predict(x$model, layers, type = "response", ...),
                      glm_iwlr = predict(x$model, layers, type = "response", ...),
                      glm_dwpr = predict(x$model, layers, type = "response", ...),
                      glm_normal = predict(x$model, layers, type = "response", ...),
                      lasso = {
                        forms <- as.formula(paste("~ 1",
                                                  paste(colnames(layers), collapse = "+"),
                                                  paste(paste("I(", colnames(layers), "^2", ")", sep = ""), 
                                                        collapse = "+"), sep = "+"))
                        
                        to_pred <- model.matrix(forms, data = as.data.frame(layers)) 
                        to_pred <- to_pred[,-1]
                        # layers_poly <- apply(layers, 2, function(x) x^2)
                        # colnames(layers_poly) <- paste0(colnames(layers), "_poly")
                        # to_pred <- cbind(layers, layers_poly)
                        p <- predict(x$model, to_pred, s = "lambda.1se", type = "response", ...)
                        # p <- predict(x$model, to_pred, type = "response", ...)
                        p[,1]
                      },
                      elasticnet = {
                        forms <- as.formula(paste("~ 1",
                                                  paste(colnames(layers), collapse = "+"),
                                                  paste(paste("I(", colnames(layers), "^2", ")", sep = ""), 
                                                        collapse = "+"), sep = "+"))
                        
                        to_pred <- model.matrix(forms, data = as.data.frame(layers)) 
                        to_pred <- to_pred[,-1]
                        # layers_poly <- apply(layers, 2, function(x) x^2)
                        # colnames(layers_poly) <- paste0(colnames(layers), "_poly")
                        # to_pred <- cbind(layers, layers_poly)
                        p <- predict(x$model, to_pred, s = "lambda.1se", type = "response", ...)
                        # p <- predict(x$model, to_pred, type = "response", ...)
                        p[,1]
                      },
                      ridge = {
                        forms <- as.formula(paste("~ 1",
                                                  paste(colnames(layers), collapse = "+"),
                                                  paste(paste("I(", colnames(layers), "^2", ")", sep = ""), 
                                                        collapse = "+"), sep = "+"))
                        
                        to_pred <- model.matrix(forms, data = as.data.frame(layers)) 
                        to_pred <- to_pred[,-1]
                        # layers_poly <- apply(layers, 2, function(x) x^2)
                        # colnames(layers_poly) <- paste0(colnames(layers), "_poly")
                        # to_pred <- cbind(layers, layers_poly)
                        p <- predict(x$model, to_pred, s = "lambda.1se", type = "response", ...)
                        # p <- predict(x$model, to_pred, type = "response", ...)
                        p[,1]
                      },
                      lgbm = predict(x$model, as.matrix(layers))
              )
            }
          }
          )


# predict.sdm_result <- function(x, layers, ...) {
#   if (class(layers)[1] == "SpatRaster") {
#     switch (x$name,
#             maxnet = predict(layers, x$model, type = "cloglog", na.rm = T, ...),
#             brt = predict(layers, x$model, type = "response", na.rm = T, ...),
#             rf_classification_ds = {
#               p <- predict(layers, x$model, type = "prob")
#               p[[2]]
#             },
#             gam_iwlr = predict(layers, x$model, type = "response", ...),
#             gam_dwpr = predict(layers, x$model, type = "response", ...),
#             gam_normal = predict(layers, x$model, type = "response", ...),
#             glm_iwlr = predict(layers, x$model, type = "response", ...),
#             glm_dwpr = predict(layers, x$model, type = "response", ...),
#             glm_normal = predict(layers, x$model, type = "response", ...),
#             lasso = {
#               lay_vals <- values(layers)
#               forms <- as.formula(paste("~ 1",
#                                         paste(colnames(lay_vals), collapse = "+"),
#                                         paste(paste("I(", colnames(lay_vals), "^2", ")", sep = ""), 
#                                               collapse = "+"), sep = "+"))
#               layers_poly <- model.matrix(forms,
#                                           model.frame(forms, as.data.frame(lay_vals),
#                                                       na.action=na.pass)) 
#               layers_poly <- layers_poly[,-1]
#               p <- predict(x$model, as.matrix(layers_poly), type = "response", ...)
#               f <- layers[[1]]
#               values(f) <- p
#               f
#             },
#             lgbm = {
#               vals <- values(layers)
#               pred <- predict(x$model, as.matrix(vals))
#               p <- layers[[1]]
#               values(p) <- pred
#               p <- mask(p, layers[[1]])
#               p
#             }
#     )
#   } else {
#     switch (x$name,
#             maxnet = predict(x$model, layers, type = "cloglog", na.rm = T, ...),
#             brt = predict(x$model, layers, type = "response", na.rm = T, ...),
#             rf_classification_ds = {
#               p <- predict(x$model, layers, type = "prob")
#               p[,2]
#             },
#             gam_iwlr = predict(x$model, layers, type = "response", ...),
#             gam_dwpr = predict(x$model, layers, type = "response", ...),
#             gam_normal = predict(x$model, layers, type = "response", ...),
#             glm_iwlr = predict(x$model, layers, type = "response", ...),
#             glm_dwpr = predict(x$model, layers, type = "response", ...),
#             glm_normal = predict(x$model, layers, type = "response", ...),
#             lasso = {
#               forms <- as.formula(paste("~ 1",
#                                         paste(colnames(layers), collapse = "+"),
#                                         paste(paste("I(", colnames(layers), "^2", ")", sep = ""), 
#                                               collapse = "+"), sep = "+"))
#               
#               to_pred <- model.matrix(forms, data = as.data.frame(layers)) 
#               to_pred <- to_pred[,-1]
#               # layers_poly <- apply(layers, 2, function(x) x^2)
#               # colnames(layers_poly) <- paste0(colnames(layers), "_poly")
#               # to_pred <- cbind(layers, layers_poly)
#               p <- predict(x$model, to_pred, type = "response", ...)
#               p[,1]
#             },
#             lgbm = predict(x$model, as.matrix(layers))
#     )
#   }
# }




#' Generate response curves for the SDM models
#'
#' @param model an object of class \code{sdm_result} returned by one of the SDM
#'   modules
#' @param layers environmental layers used to fit the model (SpatRaster format)
#' @param vars optional, a vector with names of the variables for which the
#'   response curves should be obtained.
#'   
#' @details
#' A plot method is available for the object generated by this function
#' 
#'
#' @return a data.frame containing predictions for the response curves
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
resp_curves <- function(model, layers, vars = NULL) {
  
  if (is.null(vars)) {
    vars <- names(layers)
  }
  
  env_data <- layers[[vars]]
  min_vals <- global(env_data, "min", na.rm=T)[,1]
  max_vals <- global(env_data, "max", na.rm=T)[,1]
  mea_vals <- global(env_data, "mean", na.rm=T)[,1]
  
  pred_list <- list()
  
  resp_data <- matrix(nrow = 100, ncol = length(vars))
  
  for (z in 1:length(vars)) {
    
    resp_data[,] <- rep(mea_vals, each = 100)
    
    resp_data[,z] <- seq(min_vals[z], max_vals[z], length.out = 100)
    
    resp_data <- as.data.frame(resp_data)
    names(resp_data) <- vars
    
    pred_resp <- predict(model, resp_data)
      
    pred_resp <- data.frame(response = pred_resp)
    
    pred_resp$base <- resp_data[,z]
    
    pred_list[[z]] <- pred_resp
    
  }
  
  names(pred_list) <- vars
  
  pred_list <- dplyr::bind_rows(pred_list, .id = "variable")
  
  class(pred_list) <- c("sdm_respcur", class(pred_list))
  
  return(pred_list)
}

# Plot method for the response curve generated by resp_curves
#' @export
plot.sdm_respcur <- function(response_curves) {
  
  if (!any(colnames(response_curves) %in% c("variable", "response", "base"))) {
    stop("response_curves should be an object with columns 'variable', 'response' and 'base'")
  }
  
  library(ggplot2)
  
  ggplot(response_curves) +
    geom_line(aes(x = base, y = response), color = "darkblue") +
    theme_light() +
    xlab("Value") + ylab("Response") +
    facet_wrap(~ variable, scales = "free_x")
  
}



#' Save the SDM model/metrics returned by the SDM modules
#'
#' @param model an sdm_result object returned by one of the SDM modules 
#' @param what what to save (e.g. model, cv_metrics)
#' @param where where to save (path)
#'
#' @return nothing, just saved objects
#' @export
#' 
#' @examples
#' \dontrun{
#' save_sdm(my_model, "model", "models_folder/folder")
#' }
save_sdm <- function(model, what, where) {
  to_save <- x[[what]]
  if (class(to_save)[1] != "data.frame") {
    saveRDS(to_save, file = paste0(where, "/model_", x$name, ".rds"))
  } else {
    write.csv(to_save,  paste0(where, "/", what, "_", x$name, ".csv"),
              row.names = F)
  }
  return(invisible(NULL))
}
