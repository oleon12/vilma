#' Convert Species Occurrence Points to Raster Distribution
#' @description
#' Converts species occurrence points (Species, Longitude, Latitude) to a
#' \code{vilma.dist} object that includes a template grid and (optionally)
#' raster layers of species richness and abundance per cell.
#'
#' @param points A \code{matrix} or \code{data.frame} with exactly three columns,
#'   in the order: Species, Longitude, Latitude.
#' @param crs Coordinate Reference System. Either an EPSG code as integer
#'   (e.g., \code{4326} for WGS84) or a PROJ string. Defaults to EPSG:4326.
#' @param ext A \code{terra::ext} extent defining output raster bounds. If \code{NULL}
#'   (default), the extent is computed from \code{points} and buffered slightly.
#' @param res Numeric resolution of the output grid in CRS units. Either a single
#'   value (square cells; \code{xres = yres}) or a length-2 numeric vector
#'   \code{c(xres, yres)}. Defaults to \code{1} (degrees if CRS is geographic).
#' @param doRast Logical; if \code{TRUE} (default) raster layers are created
#'   (richness and abundance). If \code{FALSE}, returns only the distribution
#'   table and the grid definition.
#' @param symmetrical Logical; if \code{TRUE}, forces square pixels by setting
#'   \code{yres <- xres} and snaps the extent so that its width/height are exact
#'   multiples of the (possibly forced) resolution. Default \code{FALSE}.
#'
#' @details
#' \strong{Resolution and symmetry:} \code{res} can be a scalar or \code{c(xres, yres)}.
#' When \code{symmetrical = TRUE}, the function enforces \code{yres = xres} and
#' snaps the extent to whole-cell multiples, ensuring symmetric (square) pixels.
#'
#' \strong{On geographic CRSs:} If the CRS is geographic (lon/lat), square cells
#' are \emph{square degrees}, not equal-area squares. For equal-area, project to
#' an appropriate (equal-area) CRS before calling this function.
#'
#' \strong{Outputs:} Cell IDs are assigned via \code{terra::cellFromXY}. If
#' \code{doRast = TRUE}, abundance is computed as record counts per cell and
#' richness as the number of unique species per cell.
#'
#' @return An object of class \code{vilma.dist} with:
#' \itemize{
#'   \item \code{distribution}: \code{data.frame} with original points and their \code{Cell} IDs.
#'   \item \code{grid}: \code{terra::SpatRaster} template grid.
#'   \item \code{r.raster}: richness raster (\verb{#} unique species per cell) or \code{NULL}.
#'   \item \code{ab.raster}: abundance raster (\verb{#} records per cell) or \code{NULL}.
#' }
#'
#' @importFrom terra ext rast cellFromXY vect rasterize is.lonlat
#' @importFrom sf st_as_sf st_coordinates
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @examples
#' \dontrun{
#' # Sample occurrence data
#' tree <- examplae_tree()
#' dist <- example_dist()
#'
#' # Rectangular cells (0.5 x 0.25 degrees)
#' d2 <- points_to_raster(dist, res = c(0.5, 0.25))
#'
#' # Force square cells with snapping (use xres = 0.25)
#' d3 <- points_to_raster(dist, res = 0.25, symmetrical = TRUE)
#'
#' # Return only table + grid (no rasters)
#' d4 <- points_to_raster(dist, doRast = FALSE)
#' }
#' @export

points_to_raster <- function(points, crs = 4326, ext = NULL, res = 1, doRast = TRUE, symmetrical = FALSE) {

  # 1. Input validation - Enhanced checks
  if (!inherits(points, c("matrix", "data.frame"))) {
    stop("'points' must be a matrix or data.frame")
  }

  if (ncol(points) != 3) {
    stop("'points' must have exactly three columns: Species, Longitude, Latitude")
  }

  # Check for missing values in coordinates
  if (any(is.na(points[, 2:3]))) {
    stop("Longitude and Latitude columns cannot contain missing values (NA)")
  }
  
  if (any(is.na(points[, 1]))) {
    stop("Species column contains missing values (NA). Please clean the input data.")
  }
  

  # 2. Set column names more robustly
  colnames(points) <- c("Sp", "Lon", "Lat")
  
  # Convert to data.frame for safer handling
  points_df <- as.data.frame(points)

  # 3. Process CRS parameter
  if (is.numeric(crs)) {
    crs_string <- paste0("EPSG:", crs)
  } else {
    crs_string <- crs
  }

  # --- resolution handling --------------------------------------------------
  # Allow res to be a scalar (x=y) or a vector c(xres, yres)
  if (length(res) == 1L) {
    xres <- res
    yres <- res
  } else if (length(res) == 2L) {
    xres <- res[1]
    yres <- res[2]
  } else {
    stop("'res' must be length 1 (square) or length 2 (xres, yres).")
  }

  # 4. Create spatial points and calculate extent
  occ_sf <- st_as_sf(points_df, coords = c("Lon", "Lat"), crs = crs_string)
  
  if (is.null(ext)) {
    ext <- ext(occ_sf)
    # Add small buffer to ensure all points are within the extent
    ext <- ext + (max(xres, yres) * 0.1)   ### NEW: use max res for buffer
  }

  # 4.1 Force symmetrical (square) pixels if requested
  if (isTRUE(symmetrical)) {
    # Make yres match xres (keep user's xres)
    if (!isTRUE(all.equal(xres, yres))) {
      message("Forcing symmetrical pixels: setting yres = xres = ", xres)
    }
    yres <- xres

    # Heads-up for geographic CRS: square degrees != square meters
    # We detect after grid creation below; message again there.
  }

  # 4.2 Snap extent to exact multiples of (xres, yres) so cells are full-sized
  xmin <- ext[1]; xmax <- ext[2]; ymin <- ext[3]; ymax <- ext[4]
  nx <- ceiling((xmax - xmin) / xres)
  ny <- ceiling((ymax - ymin) / yres)
  xmax_new <- xmin + nx * xres
  ymax_new <- ymin + ny * yres
  ext <- ext(xmin, xmax_new, ymin, ymax_new)   ### NEW: snapped extent

  # 5. Create raster grid and calculate cell IDs
  rgrid <- rast(ext, resolution = c(xres, yres), crs = crs_string)

  # If symmetrical=TRUE and grid is lon/lat, warn that "square" is in degrees
  if (isTRUE(symmetrical) && isTRUE(is.lonlat(rgrid))) {
    message("Note: symmetrical=TRUE in a geographic CRS yields square degrees, not square meters.")
  }

  cellR <- cellFromXY(rgrid, st_coordinates(occ_sf))
  outDist <- cbind(points_df, Cell = cellR)

  # 6. Create raster outputs if requested
  if (doRast) {
    occ_vect <- vect(points_df, geom = c("Lon", "Lat"), crs = crs_string)
    
    # Abundance raster (count of all records)
    ab_rast <- rasterize(occ_vect, rgrid, field = "Sp", fun = "count")
    
    # Richness raster (count of unique species)
    ri_rast <- rasterize(occ_vect, rgrid, field = "Sp", 
                                fun = function(x) length(unique(na.omit(x))))
    
    # Return vilma.dist object with raster data
    structure(
      list(
        distribution = outDist,
        grid = rgrid,
        r.raster = ri_rast,
        ab.raster = ab_rast
      ),
      class = "vilma.dist"
    )
  } else {
    # Return vilma.dist object without raster data
    structure(
      list(
        distribution = outDist,
        grid = rgrid,
        r.raster = NULL,
        ab.raster = NULL
      ),
      class = "vilma.dist"
    )
  }
}
