#' Fast ensemble ecological niche model (ENM) with cross-validated weighting
#'
#' Fits an ensemble ecological niche model (ENM) using three algorithms:
#' Maxent (via \pkg{maxnet}), Random Forest (via \pkg{randomForest}), and
#' regularized logistic regression (via \pkg{glmnet}). Presence points are
#' calibrated against background (pseudo-absence) points sampled from valid
#' raster cells (or provided by the user). Models are evaluated via stratified
#' cross-validation and combined using performance-based weights (AUC or TSS).
#' The function returns an \code{"ensemble_enm"} object containing fitted models,
#' ensemble weights, an internal prediction function, evaluation results, and
#' calibration data.
#'
#' @param occurrences A data frame with occurrence coordinates. Required columns:
#'   \code{Lon} and \code{Lat} (decimal degrees). Rows with missing predictor values
#'   after extraction are removed.
#' @param predictors A \pkg{raster} predictor object (e.g., \code{RasterStack} or
#'   \code{RasterBrick}) containing environmental predictors used for calibration
#'   and prediction.
#' @param background_points Background specification. Default = \code{10000}.
#'   If numeric, the number of random background points sampled with
#'   \code{\link[dismo]{randomPoints}} (internally capped to 30\% of valid cells of
#'   the first predictor layer). If a data frame, it is interpreted as background
#'   coordinates and used directly. If \code{NULL}, defaults to sampling 10,000 points
#'   (also subject to the valid-cell cap).
#' @param cv_folds Integer. Default = \code{5}. Number of cross-validation folds.
#'   Folds are assigned separately to presences and background (stratified). If there
#'   are too few presences or background points, \code{cv_folds} is reduced automatically.
#' @param n_cores Integer. Default = \code{1}. Number of CPU cores to use for
#'   parallel cross-validation predictions. If \code{> 1}, folds are evaluated via
#'   \pkg{foreach}/\pkg{doParallel}. Model fitting on the full dataset remains
#'   sequential.
#'
#' @param maxent_classes Character. Default = \code{"default"}. Feature classes passed
#'   to \code{\link[maxnet]{maxnet.formula}} (e.g., \code{"l"}, \code{"lq"},
#'   \code{"lqh"}, \code{"default"}).
#' @param maxent_regmult Numeric. Default = \code{1}. Regularization multiplier for
#'   \pkg{maxnet} (\code{regmult}).
#'
#' @param rf_ntree Integer. Default = \code{500}. Number of trees for Random Forest.
#' @param rf_mtry Integer or \code{NULL}. Default = \code{NULL}. Number of variables
#'   randomly sampled as candidates at each split. If \code{NULL}, set to
#'   \code{floor(sqrt(p))} where \code{p} is the number of predictors.
#' @param rf_nodesize Integer. Default = \code{1}. Minimum size of terminal nodes for
#'   Random Forest.
#'
#' @param glmnet_alpha Numeric. Default = \code{1}. Elastic net mixing parameter
#'   (\code{1} = LASSO, \code{0} = ridge).
#' @param glmnet_lambda Character. Default = \code{"lambda.1se"}. Which lambda to use
#'   for prediction from \code{\link[glmnet]{cv.glmnet}} (commonly \code{"lambda.min"}
#'   or \code{"lambda.1se"}).
#'
#' @param weight_method Character. Default = \code{"auc"}. Method to derive ensemble
#'   weights from cross-validated predictions. One of \code{"auc"}, \code{"tss"}, or
#'   any other value (falls back to equal weights).
#'
#' @param thresh_stat Character. Default = \code{"no_omission"}. Threshold statistic
#'   used by \code{\link[dismo]{threshold}} for the reported optimal threshold.
#'   One of \code{"no_omission"}, \code{"spec_sens"}, \code{"equal_sens_spec"},
#'   \code{"prevalence"}, \code{"kappa"}, or \code{"sensitivity"}.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Extract predictor values for occurrence coordinates and remove rows with
#'   missing predictor data.
#'   \item Create background points: if numeric, sample from valid raster cells using
#'   \code{\link[dismo]{randomPoints}}; if a data frame, use provided background
#'   coordinates. The number of sampled points is capped at 30\% of the valid cells
#'   in the first predictor layer to avoid requesting more points than available.
#'   \item Fit three models on the full dataset (presence/background labels):
#'   \itemize{
#'     \item Maxent via \pkg{maxnet} with \code{maxnet.formula(..., classes = maxent_classes)}
#'       and \code{regmult = maxent_regmult}.
#'     \item Random Forest via \pkg{randomForest} with \code{ntree}, \code{mtry},
#'       and \code{nodesize}.
#'     \item Regularized logistic regression via \pkg{glmnet} using \code{cv.glmnet}
#'       (binomial family) with \code{alpha = glmnet_alpha}.
#'   }
#'   \item Compute out-of-fold predictions for each model using stratified folds.
#'   Fold-level model fits are repeated within each fold (and can be parallelized).
#'   \item Derive model weights from cross-validated performance using \code{weight_method}
#'   (AUC or TSS). Weights are normalized to sum to 1.
#'   \item Evaluate the ensemble using cross-validated ensemble predictions and compute
#'   threshold-based diagnostics using \code{thresh_stat}.
#' }
#'
#' The returned \code{ensemble_predict} function supports:
#' \itemize{
#'   \item a \pkg{raster} \code{Raster*} object as \code{newdata}, returning a single-layer
#'   suitability raster (with \code{NA} where predictors are incomplete), or
#'   \item a data frame/matrix of predictor values, returning a numeric vector of suitability.
#' }
#'
#' @return An object of class \code{"ensemble_enm"} (a list) with:
#' \describe{
#'   \item{models}{List with fitted models: \code{$maxent}, \code{$rf}, \code{$glmnet}.}
#'   \item{ensemble_weights}{Named numeric vector of length 3 with weights for
#'   \code{"Maxent"}, \code{"RandomForest"}, and \code{"GLMNET"}.}
#'   \item{ensemble_predict}{A function that generates weighted ensemble suitability
#'   predictions for new data (raster or tabular).}
#'   \item{evaluation}{List with ensemble evaluation statistics:
#'     \code{auc}, \code{cor}, \code{tss}, \code{optimal_threshold}, \code{threshold_stat},
#'     and \code{confusion_matrix} (using the chosen threshold).}
#'   \item{calibration_data}{List with calibration inputs:
#'     \code{occ_env}, \code{bg_env}, \code{occ_coords}, \code{bg_coords},
#'     and \code{n_unique_presences}.}
#'   \item{parameters}{List summarizing key parameter settings for Maxent, RF, GLMNET,
#'     and cross-validation/parallel options.}
#'   \item{summary_matrix}{A data frame summarizing the modeling run, including
#'     number of unique presences, number of background points, AUC, TSS, and
#'     optimal threshold.}
#' }
#'
#' @seealso \code{\link[dismo]{randomPoints}}, \code{\link[dismo]{evaluate}},
#'   \code{\link[dismo]{threshold}}, \code{\link[maxnet]{maxnet}},
#'   \code{\link[maxnet]{maxnet.formula}}, \code{\link[randomForest]{randomForest}},
#'   \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#' \dontrun{
#' fit <- fast.ensemble.enm(occurrences, predictors,
#'                          background_points = 1000,
#'                          cv_folds = 5,
#'                          n_cores = 2,
#'                          maxent_classes = "lq",
#'                          maxent_regmult = 1,
#'                          rf_ntree = 500,
#'                          glmnet_alpha = 1,
#'                          glmnet_lambda = "lambda.1se",
#'                          weight_method = "auc",
#'                          thresh_stat = "spec_sens")
#'
#' fit$evaluation
#' print(fit$summary_matrix)
#'
#' # Predict on a RasterStack/RasterBrick:
#' suit <- fit$ensemble_predict(predictors)
#' plot(suit)
#' }
#'
#' @export

fast.ensemble.enm <- function(occurrences, predictors, 
                              background_points = 10000, 
                              cv_folds = 5, 
                              n_cores = 1,
                              
                              # Maxent parameters
                              maxent_classes = "default",
                              maxent_regmult = 1,
                              
                              # Default values for other algorithms
                              rf_ntree = 500,
                              rf_mtry = NULL,
                              rf_nodesize = 1,
                              
                              glmnet_alpha = 1,
                              glmnet_lambda = "lambda.1se",
                              
                              # Ensemble options
                              weight_method = "auc",
                              
                              thresh_stat = c("no_omission", "spec_sens", "equal_sens_spec", 
                                             "prevalence", "kappa", "sensitivity")) {
  
  # Input validation
  if(!all(c("Lon", "Lat") %in% colnames(occurrences))) {
    stop("Occurrence data must have 'Lon' and 'Lat' columns")
  }
  
  thresh_stat <- match.arg(thresh_stat)
  
  # ============================================
  # 1. Prepare data - WITH DUPLICATE REMOVAL
  # ============================================
  cat("Preparing data...\n")
  
  # Get species name if present (for summary matrix)
  species_name <- if("Sp" %in% colnames(occurrences)) {
    unique(occurrences$Sp)[1]
  } else {
    "Unknown"
  }
  
  # Remove duplicate coordinates
  orig_n <- nrow(occurrences)
  occurrences <- occurrences[!duplicated(occurrences[, c("Lon", "Lat")]), ]
  dup_removed <- orig_n - nrow(occurrences)
  if(dup_removed > 0) {
    cat(paste0("  Removed ", dup_removed, " duplicate occurrence records\n"))
  }
  
  # Check if we have enough unique presences
  if(nrow(occurrences) < 3) {
    warning(paste0("Species ", species_name, " has fewer than 3 unique occurrences. Cannot fit model."))
    
    # Return a minimal object with summary matrix
    result <- list(
      models = NULL,
      ensemble_weights = NULL,
      ensemble_predict = NULL,
      evaluation = list(
        auc = NA,
        cor = NA,
        tss = NA,
        optimal_threshold = NA,
        threshold_stat = thresh_stat,
        confusion_matrix = NULL
      ),
      calibration_data = list(
        occ_env = NULL,
        bg_env = NULL,
        occ_coords = occurrences[, c("Lon", "Lat")],
        bg_coords = NULL,
        n_unique_presences = nrow(occurrences)
      ),
      parameters = list(
        maxent = list(classes = maxent_classes, regmult = maxent_regmult),
        rf = list(ntree = rf_ntree, mtry = rf_mtry, nodesize = rf_nodesize),
        glmnet = list(alpha = glmnet_alpha, lambda = glmnet_lambda),
        parallel = list(n_cores = n_cores, cv_folds = cv_folds)
      ),
      summary_matrix = data.frame(
        species = species_name,
        n_unique_presences = nrow(occurrences),
        n_background = NA,
        n_total_records = NA,
        auc = NA,
        tss = NA,
        optimal_threshold = NA,
        status = "FAILED - insufficient presences",
        stringsAsFactors = FALSE
      )
    )
    
    class(result) <- "ensemble_enm"
    return(result)
  }
  
  # Extract occurrence data
  occ_coords <- occurrences[, c("Lon", "Lat")]
  occ_pres <- extract(predictors, occ_coords)
  
  # Remove occurrences with NA values
  complete_cases <- complete.cases(occ_pres)
  occ_pres <- occ_pres[complete_cases, ]
  occ_coords <- occ_coords[complete_cases, ]
  
  na_removed <- sum(!complete_cases)
  if(na_removed > 0) {
    cat(paste0("  Removed ", na_removed, " occurrences with NA environmental data\n"))
  }
  
  # Check again after NA removal
  if(nrow(occ_pres) < 3) {
    warning(paste0("Species ", species_name, " has fewer than 3 unique occurrences with complete environmental data. Cannot fit model."))
    
    result <- list(
      models = NULL,
      ensemble_weights = NULL,
      ensemble_predict = NULL,
      evaluation = list(
        auc = NA,
        cor = NA,
        tss = NA,
        optimal_threshold = NA,
        threshold_stat = thresh_stat,
        confusion_matrix = NULL
      ),
      calibration_data = list(
        occ_env = NULL,
        bg_env = NULL,
        occ_coords = occ_coords,
        bg_coords = NULL,
        n_unique_presences = nrow(occ_pres)
      ),
      parameters = list(
        maxent = list(classes = maxent_classes, regmult = maxent_regmult),
        rf = list(ntree = rf_ntree, mtry = rf_mtry, nodesize = rf_nodesize),
        glmnet = list(alpha = glmnet_alpha, lambda = glmnet_lambda),
        parallel = list(n_cores = n_cores, cv_folds = cv_folds)
      ),
      summary_matrix = data.frame(
        species = species_name,
        n_unique_presences = nrow(occ_pres),
        n_background = NA,
        n_total_records = NA,
        auc = NA,
        tss = NA,
        optimal_threshold = NA,
        status = "FAILED - insufficient presences after NA removal",
        stringsAsFactors = FALSE
      )
    )
    
    class(result) <- "ensemble_enm"
    return(result)
  }
  
  # ============================================
  # Generate background points
  # ============================================
  cat("Generating background points...\n")
  
  # Calculate valid (non-NA) cells in raster
  valid_cells <- sum(!is.na(getValues(predictors[[1]])))
  
  # Adjust number of background points if necessary
  if(is.null(background_points)) {
    bg_n <- 10000
  } else if(is.data.frame(background_points)) {
    bg_coords <- background_points
    bg_n <- NULL
  } else {
    bg_n <- background_points
  }
  
  if(!is.null(bg_n)) {
    # Ensure we don't request more points than exist
    bg_n <- min(bg_n, round(valid_cells * 0.3))
    
    # Generate points with warning handling
    bg_coords <- tryCatch({
      randomPoints(predictors, n = bg_n)
    }, warning = function(w) {
      # If warning, reduce further
      bg_n_adj <- round(bg_n * 0.5)
      cat(paste0("  Adjusting to ", bg_n_adj, " points\n"))
      randomPoints(predictors, n = bg_n_adj)
    })
  }
  
  bg_env <- extract(predictors, bg_coords)
  bg_coords <- bg_coords[complete.cases(bg_env), ]
  bg_env <- bg_env[complete.cases(bg_env), ]
  
  # ============================================
  # Convert to data frames and ensure numeric
  # ============================================
  occ_pres <- as.data.frame(occ_pres)
  bg_env <- as.data.frame(bg_env)
  
  occ_pres <- data.frame(lapply(occ_pres, as.numeric))
  bg_env <- data.frame(lapply(bg_env, as.numeric))
  
  # Combine presence and background
  all_data <- rbind(occ_pres, bg_env)
  all_labels <- c(rep(1, nrow(occ_pres)), rep(0, nrow(bg_env)))
  all_data <- data.frame(lapply(all_data, as.numeric))
  
  n_pres <- nrow(occ_pres)
  n_bg <- nrow(bg_env)
  n_total <- nrow(all_data)
  
  cat(paste0("  Total records: ", n_total, " (", 
             n_pres, " unique presences, ", n_bg, " background)\n"))
  
  # ============================================
  # 2. Cross-validation setup
  # ============================================
  set.seed(42)
  
  # Adjust cv_folds if too few points
  min_samples <- min(n_pres, n_bg)
  if(min_samples < cv_folds) {
    cv_folds <- max(3, min_samples)
    cat(paste0("  Adjusting cv_folds to ", cv_folds, " due to few points\n"))
  }
  
  pres_folds <- sample(rep(1:cv_folds, length.out = n_pres))
  bg_folds <- sample(rep(1:cv_folds, length.out = n_bg))
  all_folds <- c(pres_folds, bg_folds)
  
  # ============================================
  # 3. Train models
  # ============================================
  cat("\nTraining ensemble models...\n")
  
  models <- list()
  cv_predictions <- matrix(NA, nrow = nrow(all_data), ncol = 3)
  
  # ============================================
  # MODEL 1: MAXENT
  # ============================================
  cat("\n--- Training Maxent ---\n")
  cat(paste0("  Using feature classes: ", maxent_classes, "\n"))
  
  presence_data <- all_data[all_labels == 1, , drop = FALSE]
  presence_labels <- all_labels[all_labels == 1]
  
  feature_formula <- tryCatch({
    maxnet.formula(presence_data, presence_labels, classes = maxent_classes)
  }, error = function(e) {
    cat("  Warning: maxnet.formula failed, using linear formula\n")
    as.formula(paste("~", paste(colnames(all_data), collapse = " + ")))
  })
  
  models$maxent <- tryCatch({
    maxnet(all_labels, all_data, 
                   formula = feature_formula,
                   regmult = maxent_regmult)
  }, error = function(e) {
    cat("  Warning: maxnet training failed, using simplified model\n")
    # Fallback to simple logistic regression
    glm(all_labels ~ ., data = all_data, family = binomial)
  })
  
  if(inherits(models$maxent, "maxnet")) {
    cat("Maxent trained successfully\n")
  } else {
    cat("Using fallback GLM model\n")
  }
  
  # ============================================
  # MODEL 2: RANDOM FOREST
  # ============================================
  cat("\n--- Training Random Forest ---\n")
  
  if(is.null(rf_mtry)) {
    rf_mtry <- max(1, floor(sqrt(ncol(all_data))))
  }
  
  # Adjust RF for very small samples
  rf_ntree_adj <- ifelse(n_pres < 10, min(100, rf_ntree), rf_ntree)
  
  models$rf <- randomForest(all_data, as.factor(all_labels),
                                          ntree = rf_ntree_adj,
                                          mtry = rf_mtry,
                                          nodesize = min(rf_nodesize, floor(n_pres/2)),
                                          importance = TRUE)
  # FIXED: changed rf_mty to rf_mtry
  cat(paste0("  Using: ntree=", rf_ntree_adj, ", mtry=", rf_mtry, 
             ", nodesize=", min(rf_nodesize, floor(n_pres/2)), "\n"))
  
  # ============================================
  # MODEL 3: GLMNET - MODIFIED FOR SMALL SAMPLES
  # ============================================
  cat("\n--- Training GLMNET ---\n")
  
  # Check if we have enough presences for cross-validation
  if(n_pres < 10) {
    cat("WARNING: Very few unique presence points (", n_pres, "). Using glmnet without cross-validation with fixed lambda.\n")
    
    # Use glmnet without CV (fixed lambda)
    models$glmnet <- glmnet(as.matrix(all_data), all_labels,
                                   family = "binomial",
                                   alpha = glmnet_alpha,
                                   lambda = 0.01)  # Small fixed regularization
    
    cat(paste0("  Using: alpha=", glmnet_alpha, ", lambda=0.01 (fixed)\n"))
    
  } else {
    # Normal CV for GLMNET
    glm_cv_folds <- min(cv_folds, n_pres)  # Use presence count for folds
    glm_cv_folds <- max(3, glm_cv_folds)
    
    cat(paste0("  Using nfolds = ", glm_cv_folds, "\n"))
    
    models$glmnet <- cv.glmnet(as.matrix(all_data), all_labels,
                                      family = "binomial",
                                      alpha = glmnet_alpha,
                                      nfolds = glm_cv_folds)
    cat(paste0("  Using: alpha=", glmnet_alpha, ", lambda=", glmnet_lambda, "\n"))
  }
  
  # ============================================
  # Cross-validation predictions - PARALLELIZED
  # ============================================
  cat("\n--- Cross-validation (parallel) ---\n")
  
  if(cv_folds >= 3 && nrow(all_data) > cv_folds) {
    
    # Setup parallelization
    if(n_cores > 1) {
      cat(paste0("  Using ", n_cores, " cores\n"))
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      
      # Export necessary variables to workers
      clusterExport(cl, c("all_data", "all_labels", "all_folds", 
                         "maxent_classes", "maxent_regmult",
                         "rf_ntree_adj", "rf_mtry", "rf_nodesize",
                         "glmnet_alpha", "glmnet_lambda", "n_pres"),
                   envir = environment())
      
      # Load libraries in workers
      clusterEvalQ(cl, {
        library(maxnet)
        library(randomForest)
        library(glmnet)
      })
    }
    
    # Function to process a single fold
    process_fold <- function(fold) {
      fold_idx <- which(all_folds != fold)
      valid_idx <- which(all_folds == fold)
      
      fold_predictions <- matrix(NA, nrow = length(valid_idx), ncol = 3)
      
      # Maxent CV
      presence_idx_fold <- which(all_labels[fold_idx] == 1)
      if(length(presence_idx_fold) > 0) {
        presence_data_fold <- all_data[fold_idx, ][presence_idx_fold, , drop = FALSE]
        presence_labels_fold <- all_labels[fold_idx][presence_idx_fold]
        
        fold_formula <- tryCatch({
          maxnet.formula(presence_data_fold, presence_labels_fold, 
                                 classes = maxent_classes)
        }, error = function(e) {
          as.formula(paste("~", paste(colnames(all_data), collapse = " + ")))
        })
        
        m_cv <- tryCatch({
          maxnet(all_labels[fold_idx], all_data[fold_idx, ], 
                         formula = fold_formula,
                         regmult = maxent_regmult)
        }, error = function(e) {
          NULL
        })
        
        if(!is.null(m_cv)) {
          fold_predictions[, 1] <- predict(m_cv, all_data[valid_idx, ], type = "logistic")
        }
      }
      
      # RF CV
      rf_cv <- randomForest(all_data[fold_idx, ], 
                                          as.factor(all_labels[fold_idx]),
                                          ntree = rf_ntree_adj,
                                          mtry = rf_mtry,
                                          nodesize = min(rf_nodesize, floor(sum(all_labels[fold_idx] == 1)/2)))
      fold_predictions[, 2] <- predict(rf_cv, all_data[valid_idx, ], type = "prob")[,2]
      
      # GLMNET CV - handle differently for small samples
      if(n_pres < 10) {
        # For small samples, use glmnet without CV in each fold too
        glm_cv <- glmnet(as.matrix(all_data[fold_idx, ]), all_labels[fold_idx],
                                family = "binomial", 
                                alpha = glmnet_alpha,
                                lambda = 0.01)
        fold_predictions[, 3] <- predict(glm_cv, as.matrix(all_data[valid_idx, ]),
                                        type = "response")[,1]
      } else {
        # Normal CV for GLMNET in each fold
        glm_cv_folds_fold <- min(3, length(unique(all_labels[fold_idx])))
        glm_cv_folds_fold <- max(3, glm_cv_folds_fold)
        
        glm_cv <- cv.glmnet(as.matrix(all_data[fold_idx, ]), all_labels[fold_idx],
                                   family = "binomial", 
                                   alpha = glmnet_alpha, 
                                   nfolds = glm_cv_folds_fold)
        fold_predictions[, 3] <- predict(glm_cv, as.matrix(all_data[valid_idx, ]),
                                        s = glmnet_lambda, type = "response")[,1]
      }
      
      return(list(valid_idx = valid_idx, predictions = fold_predictions))
    }
    
    # Execute in parallel or sequential
    if(n_cores > 1) {
      results <- foreach(fold = 1:cv_folds, .combine = 'c') %dopar% {
        list(process_fold(fold))
      }
      stopCluster(cl)
    } else {
      results <- list()
      for(fold in 1:cv_folds) {
        results[[fold]] <- process_fold(fold)
      }
    }
    
    # Combine results
    for(res in results) {
      cv_predictions[res$valid_idx, ] <- res$predictions
    }
    
  } else {
    cat("  Not enough data for cross-validation, using training data only\n")
    cv_predictions[, 1] <- predict(models$maxent, all_data, type = "logistic")
    cv_predictions[, 2] <- predict(models$rf, all_data, type = "prob")[,2]
    
    if(n_pres < 10) {
      cv_predictions[, 3] <- predict(models$glmnet, as.matrix(all_data), 
                                     type = "response")[,1]
    } else {
      cv_predictions[, 3] <- predict(models$glmnet, as.matrix(all_data), 
                                     s = glmnet_lambda, type = "response")[,1]
    }
  }
  
  # ============================================
  # Calculate ensemble weights
  # ============================================
  cat("\nCalculating ensemble weights...\n")
  
  if(weight_method == "auc") {
    model_aucs <- apply(cv_predictions, 2, function(pred) {
      tryCatch({
        if(all(is.na(pred))) return(0)
        evaluate(p = pred[all_labels == 1], 
                        a = pred[all_labels == 0])@auc
      }, error = function(e) 0)
    })
    weights <- model_aucs / sum(model_aucs)
    
  } else if(weight_method == "tss") {
    model_tss <- apply(cv_predictions, 2, function(pred) {
      tryCatch({
        if(all(is.na(pred))) return(0)
        eval <- evaluate(p = pred[all_labels == 1], 
                                a = pred[all_labels == 0])
        thresh <- threshold(eval, stat = "spec_sens")
        eval@TPR + eval@TNR - 1
      }, error = function(e) 0)
    })
    weights <- model_tss / sum(model_tss)
    
  } else {
    weights <- rep(1/3, 3)
  }
  
  weights[is.na(weights)] <- 0
  names(weights) <- c("Maxent", "RandomForest", "GLMNET")
  
  cat("  Final weights:\n")
  cat(paste0("    Maxent: ", round(weights[1], 3), "\n"))
  cat(paste0("    RF:     ", round(weights[2], 3), "\n"))
  cat(paste0("    GLMNET: ", round(weights[3], 3), "\n"))
  
  # ============================================
  # Create ensemble prediction function
  # ============================================
  ensemble_predict <- function(newdata) {
    if(inherits(newdata, "Raster")) {
      raster_values <- getValues(newdata)
      pred_data <- as.data.frame(raster_values)
      colnames(pred_data) <- names(newdata)
      
      complete_idx <- complete.cases(pred_data)
      ensemble_pred <- rep(NA, nrow(pred_data))
      
      if(sum(complete_idx) > 0) {
        pred_data_complete <- pred_data[complete_idx, , drop = FALSE]
        
        pred_max <- tryCatch({
          if(inherits(models$maxent, "maxnet")) {
            predict(models$maxent, pred_data_complete, type = "logistic")
          } else {
            predict(models$maxent, newdata = pred_data_complete, type = "response")
          }
        }, error = function(e) rep(NA, nrow(pred_data_complete)))
        
        pred_rf <- tryCatch({
          predict(models$rf, pred_data_complete, type = "prob")[,2]
        }, error = function(e) rep(NA, nrow(pred_data_complete)))
        
        if(n_pres < 10) {
          pred_glm <- tryCatch({
            predict(models$glmnet, as.matrix(pred_data_complete), 
                    type = "response")[,1]
          }, error = function(e) rep(NA, nrow(pred_data_complete)))
        } else {
          pred_glm <- tryCatch({
            predict(models$glmnet, as.matrix(pred_data_complete), 
                    s = glmnet_lambda, type = "response")[,1]
          }, error = function(e) rep(NA, nrow(pred_data_complete)))
        }
        
        ensemble_pred[complete_idx] <- weights[1] * pred_max + 
                                        weights[2] * pred_rf + 
                                        weights[3] * pred_glm
      }
      
      result_raster <- raster(newdata)
      result_raster[] <- ensemble_pred
      names(result_raster) <- "suitability"
      return(result_raster)
      
    } else {
      pred_data <- as.data.frame(newdata)
      pred_data <- data.frame(lapply(pred_data, as.numeric))
      
      pred_max <- tryCatch({
        if(inherits(models$maxent, "maxnet")) {
          predict(models$maxent, pred_data, type = "logistic")
        } else {
          predict(models$maxent, newdata = pred_data, type = "response")
        }
      }, error = function(e) rep(NA, nrow(pred_data)))
      
      pred_rf <- tryCatch({
        predict(models$rf, pred_data, type = "prob")[,2]
      }, error = function(e) rep(NA, nrow(pred_data)))
      
      if(n_pres < 10) {
        pred_glm <- tryCatch({
          predict(models$glmnet, as.matrix(pred_data), 
                  type = "response")[,1]
        }, error = function(e) rep(NA, nrow(pred_data)))
      } else {
        pred_glm <- tryCatch({
          predict(models$glmnet, as.matrix(pred_data), 
                  s = glmnet_lambda, type = "response")[,1]
        }, error = function(e) rep(NA, nrow(pred_data)))
      }
      
      return(weights[1] * pred_max + weights[2] * pred_rf + weights[3] * pred_glm)
    }
  }
  
  # ============================================
  # Final evaluation
  # ============================================
  cat("\nCalculating final evaluation metrics...\n")
  
  ensemble_cv_pred <- cv_predictions %*% weights
  
  evaluation <- tryCatch({
    evaluate(p = ensemble_cv_pred[all_labels == 1],
                    a = ensemble_cv_pred[all_labels == 0])
  }, error = function(e) NULL)
  
  if(!is.null(evaluation)) {
    thresh <- threshold(evaluation, stat = thresh_stat)
    tss <- evaluation@TPR + evaluation@TNR - 1
    
    cat(paste0("  AUC: ", round(evaluation@auc, 3), "\n"))
    cat(paste0("  TSS: ", round(max(tss), 3), "\n"))
    cat(paste0("  Optimal threshold (", thresh_stat, "): ", round(thresh, 3), "\n"))
  } else {
    thresh <- 0.5
    tss <- NA
    cat("  Warning: Could not calculate evaluation metrics\n")
  }
  
  # ============================================
  # Create summary matrix
  # ============================================
  summary_matrix <- data.frame(
    species = species_name,
    n_unique_presences = n_pres,
    n_background = n_bg,
    n_total_records = n_total,
    auc = if(!is.null(evaluation)) evaluation@auc else NA,
    tss = if(!is.null(evaluation)) max(tss) else NA,
    optimal_threshold = thresh,
    status = "SUCCESS",
    stringsAsFactors = FALSE
  )
  
  # ============================================
  # Return results
  # ============================================
  cat("\nDone!\n")
  
  result <- list(
    models = models,
    ensemble_weights = weights,
    ensemble_predict = ensemble_predict,
    evaluation = list(
      auc = if(!is.null(evaluation)) evaluation@auc else NA,
      cor = if(!is.null(evaluation)) evaluation@cor else NA,
      tss = if(!is.null(evaluation)) max(tss) else NA,
      optimal_threshold = thresh,
      threshold_stat = thresh_stat,
      confusion_matrix = if(!is.null(evaluation)) {
        table(pred = ensemble_cv_pred > thresh, obs = all_labels)
      } else {
        NULL
      }
    ),
    calibration_data = list(
      occ_env = occ_pres,
      bg_env = bg_env,
      occ_coords = occ_coords,
      bg_coords = bg_coords,
      n_unique_presences = n_pres
    ),
    parameters = list(
      maxent = list(classes = maxent_classes, regmult = maxent_regmult),
      rf = list(ntree = rf_ntree_adj, mtry = rf_mtry, nodesize = min(rf_nodesize, floor(n_pres/2))),
      glmnet = list(alpha = glmnet_alpha, lambda = ifelse(n_pres < 10, "fixed_0.01", glmnet_lambda)),
      parallel = list(n_cores = n_cores, cv_folds = cv_folds)
    ),
    summary_matrix = summary_matrix
  )
  
  class(result) <- "ensemble_enm"
  return(result)
}
