#' Convert polygon ranges to occurrence records and optional raster layers
#'
#' Converts polygon range maps (e.g., IUCN-style distributions) into a standard
#' occurrence table (\code{Sp}, \code{Lon}, \code{Lat}) for use in the *vilma*
#' workflow. Polygons are rasterized on a global template in geographic coordinates
#' (EPSG:4326) at a user-selected resolution. Optionally, the function also returns
#' the per-species rasterized layers used to generate occurrences.
#'
#' @param pols An \code{sf} object containing polygon geometries and an attribute
#'   column with species identifiers.
#' @param sp_id Character. Name of the column in \code{pols} that contains species
#'   names/IDs. This is used to label occurrences and to name any returned raster
#'   layers.
#' @param res Character. Default = \code{"1m"}. Output raster resolution. One of:
#'   \code{"30s"} (30 arc-seconds), \code{"1m"} (1 arc-minute), \code{"2.5m"}
#'   (2.5 arc-minutes), or \code{"5m"} (5 arc-minutes). If multiple values are
#'   provided, the function defaults to \code{"1m"}.
#' @param return_raster Logical. Default = \code{TRUE}. If \code{TRUE}, return both
#'   the occurrence table and the list of rasterized polygon layers. If \code{FALSE},
#'   return only the occurrence table.
#'
#' @details
#' The function builds a global template raster spanning \code{-180:180} longitude
#' and \code{-90:90} latitude with CRS \code{EPSG:4326}. Each polygon geometry is
#' cropped to the template extent and rasterized to a binary layer.
#'
#' If rasterization yields an all-\code{NA} raster (which can occur when polygons are
#' very small relative to the chosen resolution), the function iteratively increases
#' raster resolution (halving cell size) within the polygon extent until at least one
#' raster cell is assigned a value, then resamples back to the original cropped
#' template. A warning is issued when the resolution must be refined.
#'
#' Occurrence records are extracted by converting each rasterized polygon layer to
#' points (one point per occupied raster cell) and combining them into a single
#' three-column data frame: \code{Sp} (species identifier), \code{Lon}, and \code{Lat}.
#'
#' @return If \code{return_raster = TRUE}, a list with:
#' \describe{
#'   \item{occ}{A data frame with columns \code{Sp}, \code{Lon}, \code{Lat} containing
#'   one occurrence per occupied raster cell.}
#'   \item{rasters}{A named list of \code{\link[terra]{SpatRaster}} objects (one per
#'   polygon), where list names correspond to species identifiers in \code{pols[[sp_id]]}.}
#' }
#' If \code{return_raster = FALSE}, only the occurrence data frame (\code{Sp}, \code{Lon}, \code{Lat})
#' is returned.
#'
#' @seealso \code{\link[sf]{st_geometry}}, \code{\link[terra]{rast}},
#'   \code{\link[terra]{crop}}, \code{\link[terra]{rasterize}},
#'   \code{\link[terra]{resample}}, \code{\link[raster]{rasterToPoints}},
#'   \code{\link{rast2occ}}
#'
#' @examples
#' \dontrun{
#' # pols is an sf object with polygons and a species column, e.g. "species"
#' x <- pol2occ(pols, sp_id = "species", res = "1m", return_raster = TRUE)
#' head(x$occ)
#' names(x$rasters)
#'
#' # Return only occurrences
#' occ <- pol2occ(pols, sp_id = "species", res = "2.5m", return_raster = FALSE)
#' head(occ)
#' }
#'
#' @export

pol2occ <- function(pols, sp_id, res = c("30s","1m","2.5m","5m"), return_raster = TRUE){
  
  if(!inherits(sp_id, "character")){
    stop("'sp_id' must be a character object")
  }
  
  if(any(class(pols)%in%c("sf"))){
    if(any(class(pols)%in%c("data.frame"))){
      sp_names <- pols[[sp_id]]
      pols_geom <- st_geometry(pols)
    }
  }
  
  #### Set default resolution
  
  if(length(res)>1){
    res <- "1m"
  }
  
  ### Create resolution
  
  if(res == "30s"){
    res <- 30/3600   # 0.008333333
  }
  if(res == "1m"){
    res <- 1/60     # 0.01666667
  }
  if(res == "2.5m"){
    res <- 2.5/60     # 0.04166667
  }
  if(res == "5m"){
    res <- 5/60     # 0.08333333
  }
  
  # Create template ras
  r <- rast(
    xmin = -180, xmax = 180,
    ymin = -90,  ymax = 90,
    resolution = res,
    crs = "EPSG:4326"
  )
  
  values(r) <- NA
  
  out_r <- vector("list", length = length(sp_names))
  sp_out <- c()
  occ_out <- c()
  
  pb <- txtProgressBar(min = 0, max = length(pols_geom) , style = 3)
  
  for(i in seq_along(pols_geom)){
    
    pol_i <- pols_geom[i]
    r_tmp <- crop(r, y=pol_i)
    r_i <- rasterize(vect(pol_i), r_tmp)
    
    valNA <- which(is.na(values(r_i)))
    
    # Some polygons could be to small and had not values
    # Then, the resolution must be changed
    while(length(valNA) == length(values(r_i))){
      
      warning(paste0("Polygon of: ", sp_names[i], " is to small for the selected resolution. Resolution reduced to: ", res(r_i)/2))
      
      r_agg <- rast(ext(r_i), 
                    res = res(r_i)/2,
                    crs = crs(r_i))
      r_i <- rasterize(vect(pol_i), r_agg)
      r_i <- resample(r_i, r_tmp, method = "bilinear")
      valNA <- which(is.na(values(r_i)))
      
    }
    
    out_r[[i]] <- r_i
    
    ### Occ extraction
    
    if(i == 1){
      occ_out <- rasterToPoints(raster(r_i), fun = function(x){x==1})[, c(1,2)]
      
      if(inherits(occ_out,"numeric")){occ_out <- t(as.matrix(occ_out))} # Check proper format for single occ
      
      sp_out <- rep(sp_names[i], length(occ_out[,1]))
    }else{
      occ_t <- rasterToPoints(raster(r_i), fun = function(x){x==1})[, c(1,2)]
      
      if(inherits(occ_t,"numeric")){occ_t <- t(as.matrix(occ_t))} # Check proper format for single occ
      
      occ_out <- rbind(occ_out, occ_t)
      sp_out <- c(sp_out, rep(sp_names[i], length(occ_t[,1])))
    }
    
    
    setTxtProgressBar(pb, i)
    
  }
  
  out <- data.frame(Sp = sp_out,
                    Lon = occ_out[ , 1],
                    Lat = occ_out[ , 2])
  
  names(out_r) <- sp_names
  
  if(return_raster == TRUE){
    
    out <- list(occ = out,
                rasters= out_r)
    
    return(out)
    
  }else{
    
    return(out)
    
  }
  
}
