############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
######################## SDM modules accessory functions #######################


# Internal functions for SDM fitting ----

# Check object type
#' @export
.check_type <- function(x) {
  if (class(x)[1] != "sdm_dat") {
    cli::cli_abort("sdm_data should be of type sdm_dat (generated with {.fun mp_prepare_data})")
  } else {
    return(invisible(NULL))
  }
}

# Print messages
#' @export
.cat_sdm <- function(verbosity, to_print, bg = FALSE) {
  
  if (verbosity) {
    if (bg) {
      cli::cli_alert(cli::bg_cyan(to_print))
    } else {
      cli::cli_alert(to_print)
    }
  }
  
  return(invisible(NULL))
}

# Get time of proccess
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


# Support functions ----

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

#' Generate response curves for the SDM models
#'
#' @param model an object of class \code{sdm_result} returned by one of the SDM
#'   modules
#' @param layers environmental layers used to fit the model (SpatRaster format)
#' @param sdm_data optional, an `sdm_dat` object containing the fitting data.
#'   If supplied, a column called `in_range` will be added, with `0` for values
#'   out of the training range and `1` for those within the range.
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
resp_curves <- function(model, layers, sdm_data = NULL, vars = NULL) {
  
  if (is.null(vars)) {
    vars <- names(layers)
  }
  
  if (!is.null(sdm_data)) {
    input_data <- sdm_data$training[,vars]
  } else {
    input_data <- NULL
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
    
    if (!is.null(input_data)) {
      inp_range <- range(input_data[,z])
      pred_resp$in_range <- 0
      pred_resp$in_range[pred_resp$base >= inp_range[1] & pred_resp$base <= inp_range[2]] <- 1
    }
    
    pred_list[[z]] <- pred_resp
    
  }
  
  names(pred_list) <- vars
  
  pred_list <- dplyr::bind_rows(pred_list, .id = "variable")
  rownames(pred_list) <- 1:nrow(pred_list)
  
  class(pred_list) <- c("sdm_respcur", class(pred_list))
  
  return(pred_list)
}

# Plot method for the response curve generated by resp_curves
#' @export
plot.sdm_respcur <- function(response_curves) {
  
  if (!any(colnames(response_curves) %in% c("variable", "response", "base"))) {
    stop("response_curves should be an object with columns 'variable', 'response' and 'base'")
  }
  
  require(ggplot2)
  
  if ("in_range" %in% colnames(response_curves)) {
    # Prepare additional data
    fake <- data.frame(x = c(0, 0), y = c(0, 0), rg = factor(c(0, 1), levels = c(0,1)))
    rc_temp <- response_curves[response_curves$in_range == 1,]
    ab_dat <- tapply(rc_temp$base, rc_temp$variable, function(x) data.frame(min = range(x)[1], max = range(x)[2]))
    n_ab_dat <- names(ab_dat)
    ab_dat <- dplyr::bind_rows(ab_dat, .id = "variable")
    ab_dat$variable <- n_ab_dat
    
    ggplot(response_curves) +
      geom_line(aes(x = base, y = response, color = in_range), show_guide = FALSE) +
      scale_color_gradientn(colors = c("darkblue", "#7AAC1D")) +
      geom_vline(data = ab_dat, aes(xintercept = min), color = "darkblue",
                 linetype = 3, alpha = .5) +
      geom_vline(data = ab_dat, aes(xintercept = max), color = "darkblue",
                 linetype = 3, alpha = .5) +
      geom_point(data = fake, aes(x = x, y = y, fill = rg), shape = 22, alpha = 0) +
      scale_fill_manual(values = c("darkblue", "#7AAC1D"),
                        labels = c("Extrapolation", "Fit data"), drop=FALSE, name = "") +
      guides(fill = guide_legend(override.aes = list(alpha = 1, size = 4))) +
      theme_light() +
      xlab("Value") + ylab("Response") +
      facet_wrap(~ variable, scales = "free") #free_x
  } else {
    ggplot(response_curves) +
      geom_line(aes(x = base, y = response), color = "darkblue") +
      theme_light() +
      xlab("Value") + ylab("Response") +
      facet_wrap(~ variable, scales = "free") #free_x
  }
  
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



#' Assess the importance of variables for the model
#'
#' @param model 
#' @param sdm_data 
#' @param iterations 
#'
#' @return a data.frame with variables importance
#' @export
#' 
#' @details
#' This function calculates the importance of variables using the same principle
#' of [randomForest::randomForest()] and [biomod2::bm_VariablesImportance()].
#' 
#' Once each time, a variable is shuffled and a new prediction is made. The importance of the variable
#' is given as 1-(correlation between original and new prediction). Each variable
#' is shuffled several times (controled by the argument \code{iterations}) and the 
#' function returns the mean and standard deviation. The higher the score, higher is
#' the importance of the variable for the model.
#' 
#'
#' @examples
#' \dontrun{
#' variable_importance(model, sp_data)
#' }
variable_importance <- function(model, sdm_data, iterations = 10) {
  
  pred_data <- sdm_data$training[,2:ncol(sdm_data$training)]
  
  original_pred <- predict(model, pred_data)
  
  var_import <- data.frame(
    variable = names(pred_data),
    mean = NA,
    sd = NA
  )
  
  for (i in 1:ncol(pred_data)) {
    cor_val <- c()
    
    for (z in 1:iterations) {
      working_data <- pred_data
      
      working_data[,i] <- working_data[sample.int(nrow(working_data)),i]
      
      new_pred <- predict(model, working_data)
      
      cor_val <- c(cor_val, (1 - cor(original_pred, new_pred)))
    }
    
    var_import$mean[i] <- mean(cor_val, na.rm = T)
    var_import$sd[i] <- sd(cor_val, na.rm = T)
  }
  
  var_import$mean <- round(var_import$mean, 3)
  var_import$sd <- round(var_import$sd, 3)
  
  var_import <- var_import[order(var_import$mean, decreasing = T),]
  
  var_import
}



# Ensemble functions ----

#' Ensemble model predictions
#'
#' @param ens_method which method to use for ensemble (see details)
#' @param ... two or more SpatRaster objects, or a list of two or more SpatRaster objects
#'
#' @return ensemble in SpatRaster format, with the metric and error layers
#' @export
#' 
#' @details
#' Three methods are currently implemented:
#' - mean (average of model predictions, non-weighted)
#' - median (median of model predictions)
#' - committee averaging of binary predictions
#' 
#' Remember that for the last case, binary maps should be provided.
#'
#' @examples
#' \dontrun{
#' pred_ensemble <- ensemble_models("median", pred1, pred2)
#' }
ensemble_models <- function(ens_method, ...) {
  
  av_meth <- c("mean", "median", "comm_average")
  
  if (!ens_method %in% av_meth) {
    cli::cli_abort("{.var ens_method} should be one of {.var {av_meth}}, not {.var {ens_method}}.")
  }
  
  to_ensemble <- rast(list(...))
  
  if (nlyr(to_ensemble) < 2) {
    cli::cli_abort("At least 2 predictions should be supplied for ensembling.")
  }
  
  if (ens_method == "mean" | ens_method == "median") {
    ensemble <- eval(parse(text = paste0(ens_method, "(to_ensemble)")))
    sd_ensemble <- stdev(to_ensemble)
    
    ensemble <- c(ensemble, sd_ensemble)
  } else {
    ensemble <- sum(to_ensemble)
    ensemble <- ensemble/length(to_ensemble)
  }
  
  if (nlyr(ensemble) > 1) {
    names(ensemble) <- c(ens_method, "sd_ensemble")
  } else {
    names(ensemble) <- ens_method
  }
  
  return(ensemble)
}



#' Ensemble response curves
#'
#' @param ens_method which method to use for ensemble (see details)
#' @param normalize logical, if `TRUE`, 
#' then responses are normalized in the 0-1 scale before applying the ensemble
#' @param ... two or more response curves objects generated with [resp_curves()]
#'
#' @return data.frame
#' @export
#' 
#' @details
#' Three methods are currently implemented:
#' - mean (average of model predictions, non-weighted)
#' - median (median of model predictions)
#' - committee averaging of binary predictions
#' 
#' Remember that for the last case, binary maps should be provided.
#'
#' @examples
#' \dontrun{
#' ensresp <- ensemble_respcurves("median", resp1, resp2)
#' }
ensemble_respcurves <- function(ens_method, normalize, ...) {
  
  av_meth <- c("mean", "median", "comm_average")
  
  if (!ens_method %in% av_meth) {
    cli::cli_abort("{.var ens_method} should be one of {.var {av_meth}}, not {.var {ens_method}}.")
  }
  
  if (!is.logical(normalize)) {
    cli::cli_abort("{.var normalize} should be {.var TRUE} or {.var FALSE}.")
  }
  
  to_ensemble <- list(...)
  
  if (length(to_ensemble) < 2) {
    cli::cli_abort("At least 2 predictions should be supplied for ensembling.")
  }
  
  responses <- lapply(to_ensemble, function(x){
    resp <- x$response
    if (normalize) {
      resp <- (resp-min(resp))/(max(resp)-min(resp))
    }
    resp
  })
  
  responses <- do.call("cbind", responses)
  
  if (ens_method == "mean" | ens_method == "median") {
    ensemble <- apply(responses, 1, ens_method)
    sd_ensemble <- apply(responses, 1, sd)
  } else {
    ensemble <- apply(responses, 1, sum)
    ensemble <- ensemble/ncol(responses)
  }
  
  original <- to_ensemble[[1]]
  original$response <- ensemble
  
  if (ens_method == "mean" | ens_method == "median") {
    original <- cbind(original, sd_response = sd_ensemble)
  }
  
  class(original) <- c("sdm_respcur", class(original))
  
  return(original)
}


#' Evaluate ensemble
#'
#' @param ensemble the SpatRaster object containing the ensemble for the current
#'   period. If multiple layers are available, only the first one will be used.
#' @param sdm_data `sdm_dat` object containing coordinates of the occurrence
#'   records and optionally evaluation dataset
#' @param ... further arguments passed to [eval_metrics()]
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' ensemble_eval(ensemble_layer, sp_data)
#' }
ensemble_eval <- function(ensemble, sdm_data, ...) {
  
  .check_type(sdm_data)
  
  fit_pts <- sdm_data$coord_training
  eval_pts <- sdm_data$coord_eval
  
  fit_p <- sdm_data$training$presence
  eval_p <- sdm_data$eval_data$presence
  
  pred_fit <- terra::extract(emsemble[[1]], fit_pts, ID = F)
  
  evaluate_fit <- eval_metrics(fit_p, pred_fit[,1], ...)
  
  out_eval <- list(
    fit_metrics = evaluate_fit,
    eval_metrics = NULL
  )
  
  if (!is.null(eval_p)) {
    pred_eval <- terra::extract(emsemble[[1]], eval_pts, ID = F)
    
    evaluate_eval <- eval_metrics(eval_p, pred_eval[,1], ...)
    
    out_eval$eval_metrics <- evaluate_eval
  }
  
  return(out_eval)
  
}


# Print methods ----

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

# Print the SDM options object in a nicer view
#' @export
print.sdm_options <- function(x) {
  
  if (inherits(x, "sdm_sel_options")) {
    method <- attr(x, "method")
    x <- list(x)
    names(x) <- method
  }
  
  for (i in 1:length(x)) {
    cli::cli_h1(names(x)[i])
    for (z in 1:length(x[[i]])) {
      to_print <- glue::glue("{names(x[[i]])[z]} = {ifelse(is.null(x[[i]][[z]]), '`NULL`', paste(x[[i]][[z]], collapse = ', '))}")
      cli::cli_inform(to_print)
    }
  }
  
  cli::cat_line()
  cli::cat_line(cli::col_silver("Note: when more than 1 option is available, tuning is performed for that parameter."))
  
  return(invisible(NULL))
}

# Print the results of a multiple hypothesis SDM object
#' @export
print.sdm_mult_result <- function(x) {
  
  cli::cli_h1("Multiple hypothesis testing - {x$method}")
  
  cli::cat_line()
  cli::cat_line(glue::glue("Best model: {x$best_model}"))
  cli::cat_line(glue::glue("Best variables: {paste0(x$best_variables, collapse = ', ')}"))
  cli::cat_line("Metrics:")
  
  cli::cat_line()
  metrics <- paste(names(x$metrics), "=", round(x$metrics, 2))
  names(metrics) <- rep(" ", length(metrics))
  names(metrics)[which.max(x$metrics)[1]] <- "v"
  
  cli::cli_bullets(metrics)
  
  cli::cat_line()
  cli::cat_line(cli::col_silver("# You can access the model using `'object'$model`."))
}



# Predict methods for the models returned by the SDM modules ----
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
                        layers <- terra::subset(layers, x$variables)
                        vals <- terra::values(layers)
                        vals <- as.matrix(vals)
                        pred <- predict(x$model, vals)
                        p <- layers[[1]]
                        values(p) <- pred
                        p <- mask(p, layers[[1]])
                        p
                      },
                      xgboost = {
                        layers <- terra::subset(layers, x$variables)
                        vals <- terra::values(layers)
                        vals <- as.matrix(vals)
                        pred <- predict(x$model, vals)
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
                      lgbm = predict(x$model, as.matrix(layers)),
                      xgboost = {
                        layers <- layers[,x$variables]
                        predict(x$model, as.matrix(layers))
                      }
              )
            }
          }
)
