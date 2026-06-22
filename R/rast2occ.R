#' Convert binary ENM rasters to occurrence points
#'
#' Converts binary ensemble niche model (ENM) rasters (values 1 and 0) into 
#' occurrence points by extracting one point for every cell with a value of 1. 
#' The function handles various input types and behaves differently depending 
#' on whether \code{enmRasters} is a list or a single raster object.
#'
#' @param enmRasters Binary ENM rasters to convert. Can be one of the following:
#'   \itemize{
#'     \item A \code{vilma.ENM} object (output from \code{\link{vilma.ENM}})
#'     \item A \code{RasterStack} (or \code{RasterBrick}) containing binary layers
#'     \item A \code{SpatRaster} (from the \pkg{terra} package)
#'     \item A \code{list} of raster objects (each element must be a \code{Raster*} 
#'           or \code{SpatRaster})
#'   }
#'   \strong{All layers or list elements must have names} – these are used as 
#'   species identifiers. For non‑list inputs, names must match (at least partially) 
#'   those in the \code{occ} data frame when \code{includeSkipped = TRUE}.
#' @param occ A data frame of original occurrence records with exactly three columns.
#'   The function expects either \code{Sp}, \code{Lon}, \code{Lat} or 
#'   \code{Sp}, \code{Longitude}, \code{Latitude}. Coordinates must be in decimal degrees.
#'   \strong{Note:} This argument is \emph{required} as a formal parameter, but is 
#'   \strong{only used} when \code{enmRasters} is a single raster object (not a list). 
#'   For list inputs, \code{occ} is ignored.
#' @param includeSkipped Logical. Default = \code{TRUE}. If \code{TRUE}, the function
#'   appends original occurrence records for species that are present in \code{occ}
#'   but absent from \code{enmRasters} (e.g., species skipped during ENM fitting).
#'   This option only applies when \code{enmRasters} is \strong{not} a list; for list
#'   inputs, this parameter has no effect.
#'
#' @details
#' The function operates in two distinct modes:
#' \enumerate{
#'   \item \strong{List input:} When \code{enmRasters} is a list, the function:
#'     \itemize{
#'       \item Validates that every element is a raster object.
#'       \item Iterates over each element, extracts points from cells with value 1
#'         using \code{\link[raster]{rasterToPoints}}.
#'       \item Returns a data frame with fixed column names \code{Sp}, \code{Lon}, 
#'         \code{Lat} (these are \emph{not} standardised to the column names of \code{occ}).
#'       \item Does \emph{not} use \code{occ} or \code{includeSkipped} – skipped species
#'         are not appended.
#'     }
#'   \item \strong{Other input:} When \code{enmRasters} is a \code{RasterStack},
#'     \code{RasterBrick}, \code{SpatRaster}, or \code{vilma.ENM} object (the latter is 
#'     treated as a \code{RasterStack}), the function:
#'     \itemize{
#'       \item Validates the column names of \code{occ}.
#'       \item Identifies common species between raster layers and \code{occ}, and
#'         species present only in \code{occ} (if \code{includeSkipped = TRUE}).
#'       \item Extracts points from each layer (each layer is treated as a separate species).
#'       \item Returns a data frame with columns renamed to match the column names of \code{occ}
#'         (i.e., \code{Sp}, \code{Lon}, \code{Lat} or \code{Sp}, \code{Longitude}, \code{Latitude}).
#'       \item If \code{includeSkipped = TRUE}, appends the original occurrence records
#'         for species not found in the rasters.
#'     }
#' }
#'
#' @return A data frame containing occurrence points. The column names depend on the input:
#'   \itemize{
#'     \item If \code{enmRasters} is a \strong{list}: always columns \code{Sp}, \code{Lon}, \code{Lat}.
#'     \item If \code{enmRasters} is a \strong{Other object}: columns are renamed
#'       to match the column names of \code{occ} (either \code{Sp}, \code{Lon}, \code{Lat} or
#'       \code{Sp}, \code{Longitude}, \code{Latitude}).
#'   }
#'   In both cases, the data frame contains one row per cell with value 1, plus optionally
#'   original records for skipped species (only in the single‑raster mode).
#'
#' @note
#' \itemize{
#'   \item The function assumes that binary rasters contain only values 0 and 1,
#'     where 1 indicates predicted presence.
#'   \item For large rasters with many presence cells, the output can be very large.
#'   \item Names are required for \code{enmRasters} – if missing, the function stops.
#'   \item When \code{enmRasters} is a list, elements may have different extents,
#'     projections, or origins; they are not stacked, so no alignment is forced.
#'   \item A progress bar is shown during point extraction for both modes.
#' }
#'
#' @seealso \code{\link{vilma.ENM}} for generating binary ENM rasters,
#'   \code{\link[raster]{rasterToPoints}} for the underlying extraction,
#'   \code{\link[raster]{xyFromCell}} for alternative coordinate extraction.
#'
#' @examples
#' \dontrun{
#' # Example 1: Convert a RasterStack (non‑list) with skipped species
#' out <- vilma.ENM(occ = occ_data, envs = envs, ...)
#' new_occ <- rast2occ(enmRasters = out, 
#'                     occ = occ_data, 
#'                     includeSkipped = TRUE)
#' head(new_occ)  # Column names match occ_data
#'
#' # Example 2: Convert a list of rasters (occ is ignored)
#' binary_list <- list(Species_A = raster_a, Species_B = raster_b)
#' # occ must still be supplied (even though it is not used)
#' dummy_occ <- data.frame(Sp = "dummy", Lon = 0, Lat = 0)
#' points_list <- rast2occ(enmRasters = binary_list, 
#'                         occ = dummy_occ, 
#'                         includeSkipped = FALSE)
#' # Output columns are always Sp, Lon, Lat
#' colnames(points_list)  # "Sp" "Lon" "Lat"
#'
#' # Example 3: Visualise extracted points
#' library(sp)
#' new_points <- rast2occ(out, occ_data)
#' coordinates(new_points) <- ~ Lon + Lat
#' proj4string(new_points) <- CRS("+proj=longlat +datum=WGS84")
#' plot(envs[[1]], main = "Original (red) vs Extracted (blue)")
#' points(occ_data[, c("Lon", "Lat")], col = "red", cex = 0.7)
#' points(new_points, col = "blue", cex = 0.3, pch = 3)
#' }
#'
#' @export

rast2occ <- function(enmRasters, occ, includeSkipped = TRUE){
  
  if (!inherits(enmRasters, "vilma.ENM")) {
    if (inherits(enmRasters, "RasterStack")) {
      # Already a RasterStack – leave as is
      
    } else if (inherits(enmRasters, "SpatRaster")) {
      # Convert single SpatRaster to RasterStack
      enmRasters <- stack(enmRasters)
      
    } else if (inherits(enmRasters, "list")) {
      # Validate that every element in the list is a raster object
      are_rasters <- sapply(enmRasters, function(x) {
        inherits(x, "Raster") || inherits(x, "SpatRaster")
      })
      
      if (!all(are_rasters)) {
        stop("If enmRasters is a list, all elements must be raster objects (Raster* or SpatRaster).")
      }
      # If validation passes, keep as list – stacking is not attempted
      # because extents/projections may differ.
      
    } else {
      stop("enmRaster must be either a 'vilma.ENM', 'RasterStack', 'SpatRaster', or a 'list' of Rasters")
    }
  }
  
  if(is.null(names(enmRasters))){
    stop("enmRasters' must have names")
  }
  
  ##########################################
  
  if(inherits(enmRasters, "list")){
    
    message(paste0("\n","Extracting points ...", "\n"))
    
    pb <- txtProgressBar(min = 0, max = length(enmRasters) , style = 3)
    
    for(i in seq_along(enmRasters)){
      
      spID <- names(enmRasters)[i]
      
      if(i == 1){
        occ_out <- rasterToPoints(raster(enmRasters[[i]]), fun = function(x){x==1})[, c(1,2)]
        sp_out <- rep(spID, length(occ_out[,1]))
      }else{
        occ_t <- rasterToPoints(raster(enmRasters[[i]]), fun = function(x){x==1})[, c(1,2)]
        occ_out <- rbind(occ_out, occ_t)
        sp_out <- c(sp_out, rep(spID, length(occ_t[,1])))
      }
      
      setTxtProgressBar(pb, i)
      
    }
    
    out <- data.frame(Sp = sp_out,
                      Lon = occ_out[ , 1],
                      Lat = occ_out[ , 2])
    
    return(out)
    
  }else{
    
    #Safety checks
    
    if(!all(c("Lon", "Lat") %in% colnames(occ))){
      if(!all(c("Longitude", "Latitude") %in% colnames(occ))){ 
        stop("Occurrence data must have 'lon' and 'lat', or 'Longitude' and 'Latitude' columns") 
      } 
    } 
    
    
    common_sp <- intersect(names(enmRasters), occ[ , 1])
    only_occ <- unique(occ[-which(occ[,1]%in%names(enmRasters)) ,1])
    
    if(includeSkipped == TRUE){
      if(length(common_sp)==0){
        stop("No common species found between enmRaster and occ: Species names must be the same in both")
      }
    }
    
    message(paste0("\n","Extracting points ...", "\n"))
    
    pb <- txtProgressBar(min = 0, max = length(enmRasters) , style = 3)
    
    for(i in seq_along(enmRasters)){
      
      spID <- names(enmRasters)[i]
      
      if(i == 1){
        occ_out <- rasterToPoints(raster(enmRasters[[i]]), fun = function(x){x==1})[, c(1,2)]
        sp_out <- rep(spID, length(occ_out[,1]))
      }else{
        occ_t <- rasterToPoints(raster(enmRasters[[i]]), fun = function(x){x==1})[, c(1,2)]
        occ_out <- rbind(occ_out, occ_t)
        sp_out <- c(sp_out, rep(spID, length(occ_t[,1])))
      }
      
      setTxtProgressBar(pb, i)
      
    }
    
    out <- data.frame(Sp = sp_out,
                      Lon = occ_out[ , 1],
                      Lat = occ_out[ , 2])
    
    if(includeSkipped == TRUE){
      colnames(out) <- colnames(occ)
      
      if(length(only_occ) == 0){
        message("No skipped species detected")
      }else{
        message(paste0("\n Added original occurences for ", length(only_occ), " skipped species"))
        out <- rbind(out, occ[which(occ$Sp%in%only_occ),]) 
      }
      return(out)
      
    }else{
      colnames(out) <- colnames(occ)
      return(out)
    }
  }
}


