#' Convert binary ENM rasters to occurrence points
#'
#' Converts binary ensemble niche model (ENM) rasters (with values 1 and 0) into 
#' occurrence points by extracting one point for every cell with a value of 1. 
#' The function handles various input types including \code{vilma.ENM} objects, 
#' RasterStacks, SpatRasters, or lists of rasters. Optionally, it can include 
#' original occurrence records for species that were skipped during modeling.
#'
#' @param enmRasters Binary ENM rasters to convert. Can be one of the following:
#'   \itemize{
#'     \item A \code{vilma.ENM} object (output from \code{\link{vilma.ENM}})
#'     \item A \code{RasterStack} or \code{RasterBrick} containing binary layers
#'     \item A \code{SpatRaster} (from the \pkg{terra} package)
#'     \item A named \code{list} of \code{RasterLayer} objects
#'   }
#'   All layers must be named with species identifiers that match (at least partially)
#'   those in the \code{occ} data frame.
#' @param occ A data frame of original occurrence records with exactly three columns.
#'   The function expects either \code{Sp}, \code{Lon}, \code{Lat} or 
#'   \code{Sp}, \code{Longitude}, \code{Latitude}. Coordinates must be in decimal degrees.
#'   This is used primarily when \code{includeSkipped = TRUE} to append records for
#'   species not represented in the binary rasters.
#' @param includeSkipped Logical. Default = \code{TRUE}. If \code{TRUE}, the function
#'   appends original occurrence records for species that are present in \code{occ}
#'   but absent from \code{enmRasters} (e.g., species skipped during ENM fitting due
#'   to insufficient data). If \code{FALSE}, only the points derived from the binary
#'   rasters are returned.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input formats and column names.
#'   \item Identifies common species between \code{enmRasters} and \code{occ}, as well
#'     as species present only in \code{occ} (if \code{includeSkipped = TRUE}).
#'   \item For each raster layer in \code{enmRasters}, extracts the coordinates of
#'     all cells with value 1 using \code{\link[raster]{rasterToPoints}}.
#'   \item Combines all extracted points into a single data frame with three columns:
#'     species name, longitude, and latitude.
#'   \item If \code{includeSkipped = TRUE}, appends the original occurrence records
#'     for species not represented in the binary rasters.
#' }
#'
#' @return A data frame with three columns (standardized to \code{Sp}, \code{Lon}, 
#'   \code{Lat}) containing:
#'   \itemize{
#'     \item Points extracted from binary rasters (one point per cell with value 1)
#'     \item Optionally, original occurrence records for skipped species
#'   }
#'   Column names are standardized to match the input \code{occ} format.
#'
#' @note
#' \itemize{
#'   \item The function assumes that binary rasters contain only values 0 and 1,
#'     where 1 indicates predicted presence and 0 indicates absence.
#'   \item For large rasters with many presence cells, the output data frame can
#'     become very large (potentially millions of rows). Ensure sufficient memory
#'     is available.
#'   \item When \code{includeSkipped = TRUE}, the function checks for species
#'     present in \code{occ} but missing from \code{enmRasters}. If none are found,
#'     a message informs the user.
#'   \item A progress bar shows extraction status when processing multiple rasters.
#' }
#'
#' @seealso \code{\link{vilma.ENM}} for generating binary ENM rasters,
#'   \code{\link[raster]{rasterToPoints}} for the underlying point extraction,
#'   \code{\link[raster]{xyFromCell}} for alternative coordinate extraction.
#'
#' @examples
#' \dontrun{
#' # Example 1: Convert vilma.ENM output to occurrence points
#' # First, run vilma.ENM to get binary predictions
#' out <- vilma.ENM(occ = occ_data, envs = envs, ...)
#'
#' # Convert binary rasters to points, including skipped species
#' new_occ <- rast2occ(enmRasters = out, 
#'                     occ = occ_data, 
#'                     includeSkipped = TRUE)
#'
#' # Check results
#' head(new_occ)
#' cat("Total points extracted:", nrow(new_occ))
#' cat("Species included:", unique(new_occ$Sp))
#'
#' # Example 2: Convert a list of binary rasters without original data
#' # Suppose you have a named list of binary rasters
#' binary_list <- list(Species_A = raster_a, Species_B = raster_b)
#'
#' # Convert to points (skipped species not included as occ is NULL/missing)
#' # Note: occ is still required for column name standardization
#' dummy_occ <- data.frame(Sp = "dummy", Lon = 0, Lat = 0)
#' points_only <- rast2occ(enmRasters = binary_list, 
#'                         occ = dummy_occ, 
#'                         includeSkipped = FALSE)
#'
#' # Example 3: Extract points and visualize
#' library(sp)
#' library(raster)
#'
#' new_points <- rast2occ(out, occ_data)
#'
#' # Convert to SpatialPoints for mapping
#' coordinates(new_points) <- ~ Lon + Lat
#' proj4string(new_points) <- CRS("+proj=longlat +datum=WGS84")
#'
#' # Plot original vs extracted points
#' plot(envs[[1]], main = "Original occurrences (red) vs Extracted points (blue)")
#' points(occ_data[, c("Lon", "Lat")], col = "red", cex = 0.7)
#' points(new_points, col = "blue", cex = 0.3, pch = 3)
#' }
#'
#' @export

rast2occ <- function(enmRasters, occ, includeSkipped = TRUE){
  
  #Safety checks
  
  if(!all(c("Lon", "Lat") %in% colnames(occ))){
    if(!all(c("Longitude", "Latitude") %in% colnames(occ))){ 
      stop("Occurrence data must have 'lon' and 'lat', or 'Longitude' and 'Latitude' columns") 
    } 
  } 
  
  if(!inherits(enmRasters, "vilma.ENM")){
    if(!inherits(enmRasters, "RasterStack")){
      if(inherits(enmRasters, "SpatRaster")){
        enmRasters <- stack(enmRasters)
      }else{
        if(!inherits(enmRasters, "list")){
          stop("enmRaster must be either a 'vilma.ENM', 'RasterStack', 'SpatRaster' or a 'list' with Rasters")
        }
      }
    }
  }
  
  if(is.null(names(enmRasters))){
    stop("enmRasters' must have names")
  }
  
  common_sp <- intersect(names(enmRasters), occ[ , 1])
  only_occ <- unique(occ[-which(occ[,1]%in%names(enmRasters)) ,1])
  
  if(includeSkipped == TRUE){
    if(length(common_sp)==0){
      stop("No common species found between enmRaster and occ: Species names must be the same in both")
    }
  }
  
  pb <- txtProgressBar(min = 0, max = length(enmRasters) , style = 3)
  
  for(i in seq_along(enmRasters)){
    
    spID <- names(enmRasters)[i]
    
    if(i == 1){
      occ_out <- rasterToPoints(enmRasters[[i]], fun = function(x){x==1})[, c(1,2)]
      sp_out <- rep(spID, length(occ_out[,1]))
    }else{
      occ_t <- rasterToPoints(enmRasters[[i]], fun = function(x){x==1})[, c(1,2)]
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


