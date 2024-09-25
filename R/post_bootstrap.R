#' Re-fit model with best parameter for bootstrap
#' 
#' This function enable to re-fit any model using the best parameters and a new
#' dataset, intended mainly for bootstraping applications.
#'
#' @param sdm_data the `sdm_dat` object for model fitting
#' @param algo which algorithm to use for model fitting (see [obissdm::sdm_options()])
#' @param params the best parameters. Should be a named list, with each name
#'   equivalent to the parameters names from the [obissdm::sdm_options()]. Only
#'   one value for parameter is used, so passing multiple values for the same parameter
#'   is ignored.
#' @param verbose if `TRUE`, print messages
#'
#' @return fitted model of class `sdm_result` or `sdm_esm_result` (which is a
#'   list of `sdm_result` objects)
#' @export
#'
#' @examples
#' \dontrun{
#' bp <- log_object$model_bestparams$maxent
#' # For maxent, it is a list within a list, so one more level
#' bp <- bp[[1]]
#' model_bootstrap(sdm_data, "maxent", params = bp)
#' }
model_bootstrap <- function(sdm_data, algo, params, verbose = TRUE) {

    if (verbose) cli::cli_alert_info("Running simple model for {.var {algo}}")

    if (algo == "esm") {
        base_opt <- obissdm::sdm_options("maxent")

        esm_params_table <- params$individual_parameters
        esm_variables <- dplyr::bind_rows(
            lapply(1:length(params$variable_combinations), function(x) {
                cbind(data.frame(
                    variable_1 = params$variable_combinations[[x]][1],
                    variable_2 = params$variable_combinations[[x]][2]
                ), part = x)
            })
        )
        esm_scores <- params$scores
    } else {
        base_opt <- obissdm::sdm_options()
        base_opt <- base_opt[grepl(substr(algo, 1, 2), names(base_opt))][[1]]

        params_av <- names(params)

        for (i in 1:length(params_av)) {
            base_opt[params_av[i]] <- unlist(params[params_av[i]])[1]
        }
    }

    if (algo == "maxent") {
        p <- sdm_data$training$presence
        dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]

        model_fit <- maxnet::maxnet(
            p = p,
            data = dat,
            f = maxnet::maxnet.formula(
                p = p,
                data = dat,
                classes = base_opt$features
            ),
            regmult = base_opt$remult,
            addsamplestobackground = T
        )

        model_fit <- list(name = "maxent", model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "lasso") {
        method <- base_opt[["method"]]
        alpha_param <- base_opt[["alpha_param"]]
        total_area <- base_opt[["total_area"]]
        weight_resp <- base_opt[["weight_resp"]]

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
                wt[p == 0] <- total_area / sum(p == 0)
                resp_vec <- p / wt
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

        forms <- as.formula(paste("~ 1",
            paste(colnames(dat), collapse = "+"),
            paste(paste("I(", colnames(dat), "^2", ")", sep = ""),
                collapse = "+"
            ),
            sep = "+"
        ))

        training_poly <- model.matrix(forms, data = dat)
        training_poly <- training_poly[, -1]

        model_fit <- glmnet::cv.glmnet(
            x = training_poly,
            y = resp_vec,
            family = fam,
            alpha = alpha_param,
            weights = wt,
            nfolds = 5,
            foldid = NULL,
            type.measure = meas,
            keep = FALSE
        )

        if (final_alpha_param > 0 & final_alpha_param < 1) {
    model_type <- "elasticnet"
  } else if (final_alpha_param == 1) {
    model_type <- "lasso"
  } else {
    model_type <- "ridge"
  }

        model_fit <- list(name = model_type, model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "brt") {
        method <- base_opt[["method"]]
        n_trees <- base_opt[["n_trees"]]
        t_depth <- base_opt[["t_depth"]]
        lr <- base_opt[["lr"]]
        bag_fraction <- base_opt[["bag_fraction"]]
        weight_resp <- base_opt[["weight_resp"]]

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
                wt[p == 0] <- total_area / sum(p == 0)
                dat$presence <- p / wt
            } else {
                pres <- length(p[p == 1])
                bkg <- length(p[p == 0])
                wt <- ifelse(p == 1, 1, pres / bkg)
            }
        } else {
            wt <- NULL
        }

        model_fit <- gbm::gbm(presence ~ .,
            distribution = "poisson",
            data = dat,
            weights = wt,
            n.trees = n_trees,
            interaction.depth = t_depth,
            shrinkage = lr,
            bag.fraction = bag_fraction,
            cv.folds = 0,
            n.minobsinnode = 5,
            verbose = FALSE
        )
        model_fit <- list(name = "brt", model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "rf") {
        n_trees <- base_opt[["n_trees"]]
        method <- base_opt[["method"]]
        type <- base_opt[["type"]]
        mtry <- base_opt[["mtry"]]
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
            mtry <- unlist(lapply(mtry, function(x) {
                if (method == "classification") {
                    switch(x,
                        default = floor(sqrt(preds)),
                        double = floor(sqrt(preds)) * 2,
                        total = preds
                    )
                } else {
                    switch(x,
                        default = max(floor(ncol(preds) / 3), 1),
                        double = max(floor(ncol(preds) / 3), 1) * 2,
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
                mtry <- max(floor(ncol(preds) / 3), 1)
            }
        }

        if (type == "down-sampled") {
            pres <- sum(p)
            smpsize <- c("0" = pres, "1" = pres)
        } else {
            smpsize <- nrow(dat)
        }

        model_fit <- randomForest::randomForest(
            formula = presence ~ .,
            data = dat,
            ntree = n_trees,
            mtry = mtry,
            sampsize = smpsize,
            replace = TRUE
        )
        model_fit <- list(name = paste0("rf_",
                  method, "_",
                  ifelse(type == "down-sampled", "ds", "normal")), model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "glm") {
        method <- base_opt[["method"]]
        weight_resp <- base_opt[["weight_resp"]]
        quadratic <- base_opt[["quadratic"]]
        total_area <- base_opt[["total_area"]]

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
                wt[p == 0] <- total_area / sum(p == 0)
                dat$presence <- p / wt
            } else {
                pres <- length(p[p == 1])
                bkg <- length(p[p == 0])
                wt <- ifelse(p == 1, 1, pres / bkg)
            }
        } else {
            wt <- NULL
        }
        forms <- as.formula(
            paste(
                "presence ~",
                paste(
                    var_names,
                    collapse = "+"
                ),
                ifelse(quadratic, paste0("+ ", paste0(
                    "I(", var_names, "^2)",
                    collapse = " + "
                ), ""))
            )
        )

        if (method == "dwpr") {
            model_fit <- glm(forms, family = poisson(), data = dat, weights = wt)
        } else {
            model_fit <- glm(forms, family = binomial(), data = dat, weights = wt)
        }
        model_fit <- list(name = paste0("glm_", ifelse(method == "naive", "normal", method)), model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "gam") {
        method <- base_opt[["method"]]
        weight_resp <- base_opt[["weight_resp"]]
        k_val <- base_opt[["k_val"]]
        total_area <- base_opt[["total_area"]]
        select <- base_opt[["select"]]

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
                wt[p == 0] <- total_area / sum(p == 0)
                dat$presence <- p / wt
            } else {
                pres <- length(p[p == 1])
                bkg <- length(p[p == 0])
                wt <- ifelse(p == 1, 1, pres / bkg)
            }
        } else {
            wt <- NULL
        }

        var_names <- colnames(dat)[colnames(dat) != "presence"]
        forms <- as.formula(
            paste("presence ~", paste(
                "s(", var_names, ", k=", k_val, ", bs='cr')",
                collapse = "+"
            ))
        )

        model_fit <- mgcv::bam(forms,
            family = fam,
            data = dat,
            weights = wt,
            method = "fREML",
            select = select
        )
        model_fit <- list(name = paste0("gam_", ifelse(method == "naive", "normal", method)), model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "xgboost") {
        shrinkage <- base_opt[["shrinkage"]]
        gamma <- base_opt[["gamma"]]
        depth <- base_opt[["depth"]]
        rounds <- base_opt[["rounds"]]
        scale_pos_weight <- base_opt[["scale_pos_weight"]]

        # Separate data
        p <- sdm_data$training$presence
        dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]


        model_fit <- xgboost::xgboost(
            data = as.matrix(dat),
            label = p,
            max_depth = depth,
            gamma = gamma,
            scale_pos_weight = scale_pos_weight,
            nrounds = rounds,
            learning_rate = shrinkage,
            verbose = 0,
            nthread = 1,
            objective = "binary:logistic"
        )
        model_fit <- list(name = "xgboost", model = model_fit,
                          variables = colnames(dat)[colnames(dat) != "presence"])
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "lgbm") {
        objective <- base_opt[["objective"]]
        max_dep <- base_opt[["max_dep"]]
        lr <- base_opt[["lr"]]
        num_leaves <- base_opt[["num_leaves"]]
        feature_fraction <- base_opt[["feature_fraction"]]
        bagging_fraction <- base_opt[["bagging_fraction"]]
        early_stopping_round <- base_opt[["early_stopping_round"]]
        weight_resp <- base_opt[["weight_resp"]]

        # Separate data
        p <- sdm_data$training$presence
        dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]

        if (weight_resp) {
            pres <- length(p[p == 1])
            bkg <- length(p[p == 0])
            wt <- ifelse(p == 1, 1, pres / bkg)
        } else {
            wt <- NULL
        }

        to_include <- 1:nrow(dat)
        to_include <- to_include[-to_hold]

        model_fit <- lightgbm::lgb.train(
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
                # is_unbalance = T
                # feature_fraction = .8
            ),
            valids = list(test = lightgbm::lgb.Dataset(as.matrix(dat)[to_hold, ],
                label = p[to_hold], weight = wt[to_hold]
            )),
            data = lightgbm::lgb.Dataset(as.matrix(dat)[to_include, ],
                label = p[to_include], weight = wt[to_include]
            )
        )
        model_fit <- list(name = "lgbm", model = model_fit)
        class(model_fit) <- c("sdm_result", class(model_fit))
    } else if (algo == "esm") {
        bivar_df <- esm_variables

        n_vars <- nrow(bivar_df)

        models_list <- lapply(1:n_vars, function(x) NULL)

        for (va in 1:n_vars) {
            if (verbose) cli::cat_line(paste("Running ESM variable combination", va, "out of", n_vars))

            p <- sdm_data$training$presence
            dat <- sdm_data$training[, unlist(bivar_df[va, 1:2])]

            models_list[[va]] <- maxnet::maxnet(
                p = p,
                data = dat,
                f = maxnet::maxnet.formula(
                    p = p,
                    data = dat,
                    classes = esm_params_table$features[va]
                ),
                regmult = esm_params_table$remult[va],
                addsamplestobackground = T
            )

            models_list[[va]] <- list(name = "maxent", model = models_list[[va]])
            class(models_list[[va]]) <- c("sdm_result", class(models_list[[va]]))
        }

        model_fit <- models_list
        attr(model_fit, "scores") <- esm_scores
        class(model_fit) <- c("sdm_esm_result", class(model_fit))
    } else {
        cli::cli_abort("Algorithm {.arg {algo}} not implemented.")
    }

    return(model_fit)
}