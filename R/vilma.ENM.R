#' Fit Vilma ensemble ENMs for multiple species with optional M-area definition
#'
#' Fits ensemble ecological niche models (ENMs) for each species in a multi-species
#' occurrence table using \code{\link{fast.ensemble.enm}}. For each species, the
#' function can define an accessible area (M; sensu BAM) from the full study area,
#' from all occurrences pooled, or separately per taxon using a minimum convex
#' polygon (MCP) or buffer workflow. Optionally, collinear predictors within M are
#' removed using a VIF-based stepwise procedure. Species with fewer than 3 unique
#' occurrences are automatically skipped. Finally, an estimated species richness
#' raster is computed by summing the binary predictions across all successfully
#' modeled species.
#'
#' @param occ A data frame of occurrences with exactly three columns. The function
#'   expects either \code{Sp}, \code{Lon}, \code{Lat} or \code{Sp}, \code{Longitude},
#'   \code{Latitude}. Coordinates must be in decimal degrees. Columns are internally
#'   standardized to \code{Sp}, \code{Lon}, \code{Lat}. Duplicate coordinates are
#'   handled within \code{\link{fast.ensemble.enm}}.
#' @param envs Environmental predictors. Can be a \pkg{raster} \code{RasterStack},
#'   \code{RasterBrick}, or a \pkg{terra} \code{SpatRaster}. If a \code{SpatRaster}
#'   is provided, it is coerced to a \code{RasterStack} via \code{raster::stack()}
#'   for compatibility with downstream functions.
#' @param Marea Character. Default = \code{"byTaxon"}. Defines how the accessible area
#'   (M) is constructed:
#'   \describe{
#'     \item{\code{"byTaxon"}}{Build M separately for each species. For species with
#'       \code{>5} occurrences, uses an MCP buffered by \code{buffer.dist} meters.
#'       For species with \code{≤5} occurrences, uses a circular buffer around all
#'       points (also of radius \code{buffer.dist}) to avoid unstable MCPs.}
#'     \item{\code{"byArea"}}{Build a single M from all occurrences pooled and apply
#'       it to all species. Currently implemented as a placeholder; requires
#'       additional development.}
#'     \item{\code{"user"}}{Use \code{envs} as provided (no M-area cropping).}
#'   }
#' @param background_points Integer. Default = \code{10000}. Number of background
#'   points passed to \code{\link{fast.ensemble.enm}}. Automatically capped at 30\%
#'   of valid raster cells.
#' @param cv_folds Integer. Default = \code{5}. Number of cross-validation folds
#'   passed to \code{\link{fast.ensemble.enm}}. Automatically reduced for small
#'   sample sizes.
#' @param n_cores Integer. Default = \code{1}. Number of CPU cores for parallel
#'   cross-validation inside \code{\link{fast.ensemble.enm}}.
#' @param test.VIF Logical. Default = \code{TRUE}. If \code{TRUE}, remove collinear
#'   predictors inside each M-area. For species with \code{>5} occurrences, applies
#'   \code{\link[usdm]{vifstep}} to environmental values extracted at occurrence
#'   locations, then uses \code{\link[usdm]{exclude}} to remove correlated variables.
#'   For species with \code{≤5} occurrences, VIF is skipped due to insufficient
#'   sample size.
#' @param buffer.dist Numeric. Default = \code{10000}. Buffer distance in meters
#'   for MCP expansion or point buffers when constructing M areas. Used in
#'   \code{.envs_m_area} (for MCPs) and \code{.envs_m_buffer} (for point buffers).
#' @param maxent_classes Character. Default = \code{"default"}. Maxent feature classes
#'   passed to \code{\link{fast.ensemble.enm}}. Can be \code{"default"}, \code{"l"},
#'   \code{"lq"}, \code{"lh"}, \code{"lqhp"}, etc.
#' @param maxent_regmult Numeric. Default = \code{1}. Regularization multiplier for
#'   Maxent passed to \code{\link{fast.ensemble.enm}}.
#' @param rf_ntree Integer. Default = \code{500}. Number of trees in Random Forest
#'   passed to \code{\link{fast.ensemble.enm}}.
#' @param rf_mtry Integer or \code{NULL}. Default = \code{NULL}. Number of candidate
#'   variables at each split. If \code{NULL}, determined internally by
#'   \code{\link{fast.ensemble.enm}} as \code{floor(sqrt(p))}.
#' @param rf_nodesize Integer. Default = \code{1}. Minimum node size in Random Forest
#'   passed to \code{\link{fast.ensemble.enm}}.
#' @param glmnet_alpha Numeric. Default = \code{1}. Elastic net mixing parameter
#'   (\code{1}=LASSO, \code{0}=ridge) passed to \code{\link{fast.ensemble.enm}}.
#' @param glmnet_lambda Character. Default = \code{"lambda.1se"}. Lambda choice for
#'   predictions from \code{\link[glmnet]{cv.glmnet}} passed to
#'   \code{\link{fast.ensemble.enm}}. For species with \code{<10} unique presences,
#'   a fixed lambda of 0.01 is used without cross-validation.
#' @param weight_method Character. Default = \code{"auc"}. Weighting method used inside
#'   \code{\link{fast.ensemble.enm}} to combine models. Options: \code{"auc"},
#'   \code{"tss"}, or any other value (falls back to equal weights).
#' @param thresh_stat Character. Threshold statistic passed to
#'   \code{\link{fast.ensemble.enm}} and used to select an optimal threshold for
#'   binarizing the ensemble prediction. One of \code{"no_omission"},
#'   \code{"spec_sens"}, \code{"equal_sens_spec"}, \code{"prevalence"},
#'   \code{"kappa"}, or \code{"sensitivity"}.
#'
#' @details
#' For each species, \code{vilma.ENM}:
#' \enumerate{
#'   \item Subsets \code{occ} to the focal species.
#'   \item Checks if the species has at least 3 unique occurrences. If not, skips it
#'     with a message and continues to the next species.
#'   \item Constructs the accessible area (M) according to \code{Marea}:
#'     \itemize{
#'       \item \code{"byTaxon"}: Uses \code{.envs_m_area} (MCP + buffer) for species
#'         with >5 occurrences, and \code{.envs_m_buffer} (circular buffers around
#'         points) for species with ≤5 occurrences to avoid unstable MCPs.
#'       \item \code{"byArea"}: Placeholder for pooled M area.
#'       \item \code{"user"}: Uses full \code{envs} without cropping.
#'     }
#'   \item Optionally removes collinear predictors in M using a VIF-based filter
#'     (if \code{test.VIF = TRUE}). For species with >5 occurrences, applies
#'     \code{vifstep} to environmental values at occurrence points; for ≤5
#'     occurrences, VIF is skipped due to insufficient sample size.
#'   \item Fits an ensemble ENM using \code{\link{fast.ensemble.enm}} with all
#'     specified parameters.
#'   \item Produces two outputs: (i) a continuous suitability raster (ensemble
#'     prediction), and (ii) a binary prediction raster using the model's selected
#'     optimal threshold.
#' }
#' 
#' After processing all species, the function:
#' \itemize{
#'   \item Aligns all binary rasters to a common extent (the maximum extent
#'     covered by any successfully modeled species).
#'   \item Sums the binary rasters to generate an estimated species richness map.
#'   \item Sets pixels with zero summed value to NA (indicating areas not predicted
#'     suitable for any species).
#' }
#'
#' @return An object of class \code{"vilma.ENM"} containing a named list with four elements:
#' \describe{
#'   \item{Prediction}{A list of continuous suitability rasters (one per successfully
#'     modeled species) from the ensemble prediction. Names correspond to species names.}
#'   \item{Boolean}{A list of binary rasters (one per successfully modeled species)
#'     indicating suitability above the optimal threshold (based on \code{thresh_stat}).}
#'   \item{Est_richness}{A single raster layer of class \code{RasterLayer} representing 
#'     estimated species richness, calculated as the sum of all binary predictions across 
#'     successfully modeled species, aligned to a common extent. Areas with no predicted 
#'     suitability for any species are set to NA.}
#'   \item{Skipped}{A character vector with names of species that were skipped due to
#'     insufficient occurrences (<3 unique points).}
#' }
#'
#' @note
#' \itemize{
#'   \item The \code{"byArea"} option is currently a placeholder and requires
#'     additional implementation to function properly.
#'   \item For species with very few occurrences (<10 unique presences),
#'     \code{\link{fast.ensemble.enm}} issues warnings and adjusts modeling
#'     parameters (e.g., uses GLMNET without CV). Results for such species should
#'     be interpreted with caution.
#'   \item VIF-based variable selection requires the \pkg{usdm} package. If not
#'     available, the function will attempt to install it.
#'   \item The estimated richness raster is based solely on the binary predictions
#'     of successfully modeled species; skipped species are not included.
#'   \item All binary rasters are extended to a common extent before summation,
#'     ensuring that the richness estimate covers the maximum geographic range of
#'     any modeled species.
#'   \item The function returns an S3 object of class \code{"vilma.ENM"} to enable
#'     custom printing and plotting methods in future package versions.
#' }
#'
#' @seealso \code{\link{fast.ensemble.enm}}, \code{\link[maxnet]{maxnet}},
#'   \code{\link[randomForest]{randomForest}}, \code{\link[glmnet]{cv.glmnet}},
#'   \code{\link[usdm]{vifstep}}, \code{\link[usdm]{exclude}}
#'
#' @examples
#' \dontrun{
#' # Basic usage with per-species M areas
#' out <- vilma.ENM(
#'   occ = occ_data,
#'   envs = envs,
#'   Marea = "byTaxon",
#'   background_points = 5000,
#'   cv_folds = 5,
#'   n_cores = 2,
#'   test.VIF = TRUE,
#'   buffer.dist = 10000,
#'   maxent_classes = "lq",
#'   maxent_regmult = 1,
#'   rf_ntree = 500,
#'   glmnet_alpha = 1,
#'   glmnet_lambda = "lambda.1se",
#'   weight_method = "auc",
#'   thresh_stat = "spec_sens"
#' )
#'
#' # Check object class and structure
#' class(out)  # "vilma.ENM"
#' names(out)
#'
#' # Check results
#' names(out$Prediction)       # Successfully modeled species
#' out$Skipped                 # Species skipped due to insufficient data
#' 
#' # Visualize outputs
#' plot(out$Prediction[[1]])    # Continuous suitability for first species
#' plot(out$Boolean[[1]])       # Binary prediction for first species
#' plot(out$Est_richness)       # Estimated species richness map
#' 
#' # Access specific species results
#' sp_name <- names(out$Prediction)[1]
#' sp_continuous <- out$Prediction[[sp_name]]
#' sp_binary <- out$Boolean[[sp_name]]
#' 
#' # Summary of modeling results
#' cat("Successfully modeled:", length(out$Prediction), "species\n")
#' cat("Skipped:", length(out$Skipped), "species\n")
#' }
#'
#' @export

vilma.ENM <- function(occ, envs, Marea = c("byTaxon", "byArea", "user"),
                      background_points = 10000, cv_folds = 5, n_cores = 1,
                      test.VIF = TRUE, buffer.dist = 10000,
                      maxent_classes = "default", maxent_regmult = 1,
                      rf_ntree = 500, rf_mtry = NULL, rf_nodesize = 1,
                      glmnet_alpha = 1, glmnet_lambda = "lambda.1se",
                      weight_method = "auc",
                      thresh_stat = c("no_omission", "spec_sens", "equal_sens_spec",
                                      "prevalence", "kappa", "sensitivity")) {

  # ----------------------------
  # Safety checks
  # ---------------------------

  if(!all(c("Lon", "Lat") %in% colnames(occ))){
    if(!all(c("Longitude", "Latitude") %in% colnames(occ))){ 
      stop("Occurrence data must have 'lon' and 'lat', or 'Longitude' and 'Latitude' columns") 
    } 
  } 
  
  if(length(colnames(occ))!=3){ 
    stop("Occurrences data must have three columns: 'Sp', 'Lon', and 'Lat'") 
  } 
  
  if(!inherits(envs,"RasterStack")){ 
    if(inherits(envs, "SpatRaster")){ 
    envs <- stack(envs) 
    }else{ 
      stop("'envs' object must be either 'RasterBrick' or 'SpatRaster' class")
    } 
  }

  Marea <- match.arg(Marea)
  
  .envs_m_area <- function(occ, envs, dist = 10000){
  
                           m_area <- SpatialPoints(coords = occ[,c("Lon","Lat")], proj4string = CRS(proj4string(envs)))
                           sp_mcp <- mcp(m_area, percent = 100) %>% 
                           st_as_sf() %>%
                           st_make_valid() %>%
                           st_buffer(dist = dist)
  
                           m_envs <- crop(envs, y = sp_mcp) %>% mask(mask = sp_mcp)
                           return(m_envs)
   }

  .envs_m_buffer <- function(occ, envs, dist = 10000){
  
                             pts_sf <- st_as_sf(occ, coords = c("Lon", "Lat"), crs = 4326)
                             sp_mcp <- st_buffer(pts_sf, dist = dist)
  
                              m_envs <- crop(envs, y = sp_mcp) %>% mask(mask = sp_mcp)
                              return(m_envs)
   }

  # ----------------------------
  # Initial M area selection
  # ----------------------------

  if (Marea == "user") {
    m_envs <- envs
  }

  if (Marea == "byArea") {
    m_envs <- .envs_m_area(occ = occ, envs = envs)
  }

  # ----------------------------
  # Iterating over species
  # ----------------------------

  sp.list <- unique(occ$Sp)
  models <- vector("list", length(sp.list))
  preds  <- vector("list", length(sp.list))
  bools  <- vector("list", length(sp.list))
  skipped <- character(0)

  for (i in seq_along(sp.list)) {

    sp.occ <- occ[occ$Sp %in% sp.list[i], , drop = FALSE]
    
    # ============================================
    # Check for minimum occurrences (NUEVO)
    # ============================================
    # Remove duplicates for accurate count
    sp.occ.unique <- sp.occ[!duplicated(sp.occ[, c("Lon", "Lat")]), ]
    
    if(nrow(sp.occ.unique) < 3) {
      cat(paste0("  Species ", sp.list[i], " has fewer than 3 unique occurrences. Skipping.\n"))
      skipped <- c(skipped, sp.list[i])
      next
    }

    if (Marea == "byTaxon") {
      if(length(sp.occ$Sp)>5){
        m_envs <- .envs_m_area(occ = sp.occ, envs = envs, dist = buffer.dist)
      }else{
        m_envs <- .envs_m_buffer(occ = sp.occ, envs = envs, dist = buffer.dist)
      } 
    }

    if (isTRUE(test.VIF)) {
      if(length(sp.occ$Sp)>5){
        m_envs <- suppressWarnings(
          exclude(m_envs, vifstep(extract(m_envs, sp.occ[, c("Lon", "Lat")]))))
      }else{
        m_envs <- m_envs
      }
    }

    cat(paste0("Running model: ", sp.list[i], "\n"))
    cat(paste0("  - Species ", i, " of ", length(sp.list), "\n\n"))

    model <- fast.ensemble.enm(
      occurrences = sp.occ,
      predictors = m_envs,
      n_cores = n_cores,
      background_points = background_points,
      cv_folds = cv_folds,
      maxent_classes = maxent_classes,
      maxent_regmult = maxent_regmult,
      rf_ntree = rf_ntree,
      rf_mtry = rf_mtry,
      rf_nodesize = rf_nodesize,
      glmnet_alpha = glmnet_alpha,
      glmnet_lambda = glmnet_lambda,
      weight_method = weight_method,
      thresh_stat = thresh_stat
    )

    models[[i]] <- model$ensemble_predict(m_envs)
    bools[[i]]  <- models[[i]] >= model$evaluation$optimal_threshold
  }

  # Remove NULL entries from lists
  names(models) <- sp.list
  names(bools) <- sp.list
  
  models <- models[!sapply(models, is.null)]
  bools <- bools[!sapply(bools, is.null)]
  
  # Richness
  
  ### Fin common extent
  
  preds <- bools
  
  xmin_comun <- min(sapply(preds, function(x) xmin(x)))
  xmax_comun <- max(sapply(preds, function(x) xmax(x)))
  ymin_comun <- min(sapply(preds, function(x) ymin(x)))
  ymax_comun <- max(sapply(preds, function(x) ymax(x)))
  
  # Common extent
  ext_total <- extent(xmin_comun, xmax_comun, ymin_comun, ymax_comun)
  
  # Raster template
  template <- raster(ext_total, resolution = res(preds[[1]])[1], crs = crs(preds[[1]]))
  
  #Richness out raster
  sp_rich <- raster(template)
  values(sp_rich) <- 0
  
  # Sumar todos los rasters
  for(i in seq_along(preds)){
  
    # Extender cada raster al template (rellena con NA donde no hay datos)
    rtmp <- extend(preds[[i]], template, value = NA)
  
    # Reemplazar NA con 0 para la suma
    rtmp[is.na(rtmp)] <- 0
  
    # Sumar
    sp_rich <- sp_rich + rtmp
  }

  values(sp_rich)[which(values(sp_rich)==0)] <- NA

 return(
   structure(
     list(Prediction = models,
                Boolean = bools,
                Est_richness = sp_rich,
                Skipped = skipped),
     class = "vilma.ENM")
   
   )
              
}
