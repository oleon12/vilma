#' Write a \code{vilma.region} Object to Disk
#'
#' Exports a \code{vilma.region} object produced by \code{\link{vilma.regionalize}}
#' by writing its tabular outputs (\code{cell.info} and \code{group.info}) and its
#' raster outputs (region IDs and mean within-region dissimilarity) to disk. File
#' names are generated from the user-provided \code{file} prefix.
#'
#' @param vilma.region A \code{vilma.region} object (see \code{\link{vilma.regionalize}}).
#' @param file Character string giving the prefix for output file names.
#'   The function appends appropriate suffixes and extensions automatically.
#' @param raster.format Character; output format for raster layers.
#'   One of \code{"tif"}, \code{"grd"}, or \code{"img"}.
#' @param overwrite Logical; whether to overwrite existing files. Default = \code{TRUE}.
#'
#' @details
#' The function writes four files:
#' \itemize{
#'   \item \code{<file>cell_info.csv} — per-cell assignments, including the cell index,
#'   region (group) ID, and mean within-region dissimilarity for each cell.
#'   \item \code{<file>group_info.csv} — per-group summary, including group size and
#'   mean within-region dissimilarity.
#'   \item \code{<file>_groups.<format>} — raster of integer region (group) IDs.
#'   \item \code{<file>_mean_groups.<format>} — raster of mean within-region
#'   dissimilarity per cell (singleton groups may appear as \code{-1} if encoded that way
#'   in \code{vilma.region} for visualization).
#' }
#'
#' Output paths are printed to the console, and the function invisibly returns
#' normalized (absolute) paths.
#'
#' @return Invisibly returns a named list with normalized file paths:
#' \describe{
#'   \item{cell.info.csv}{Path to the per-cell CSV file.}
#'   \item{group.info.csv}{Path to the per-group CSV file.}
#'   \item{group.raster}{Path to the region-ID raster.}
#'   \item{group.mean}{Path to the mean within-region dissimilarity raster.}
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{\link{vilma.regionalize}}
#'
#' @examples
#' \dontrun{
#' reg <- vilma.regionalize(beta)
#' write.vilma.region(reg, file = "results/region_export", raster.format = "tif")
#' }
#'
#' @export

write.vilma.region <- function(vilma.region, file, raster.format = c("tif","grd","img"), overwrite = TRUE){


  if(!inherits(vilma.region, "vilma.region")){
    stop("Input is not a vilma.region object. See vilma.regionalize().")
  }
  
  if(!inherits(file,"character")){
    stop("File name must be a character object.")
  }
  
  if(length(raster.format) > 1){
    raster.format <- "tif"
    cat("\n")
    message("TIF format selected")
    cat("\n")
  }
  
  ###########################################################
  #                  Output file names                      #
  ###########################################################
  
  cell.info.name <- paste0(file, "_cell_info.csv")
  group.info.name <- paste0(file, "_group_info.csv")
  group.name <- paste0(file, "_groups.", raster.format)
  mean.name <- paste0(file, "_mean_groups.", raster.format)
  
  ###########################################################
  #                     Save files                          #
  ###########################################################
  
  write.csv(x = vilma.region$cell.info, file = cell.info.name,
            row.names = FALSE)
  write.csv(x = vilma.region$group.info, file = group.info.name,
            row.names = FALSE)
            
  suppressMessages(writeRaster(x = vilma.region$raster$group.raster,
                               filename = group.name, overwrite = overwrite))
                               
  suppressMessages(writeRaster(x = vilma.region$raster$mean.group.raster,
                               filename = mean.name, overwrite = overwrite))
  
  ###########################################################
  #                  Output messages                        #
  ###########################################################
  
  cat("Saved files: \n\n")
  cat(paste0("Cell groups: ", getwd(), "/", cell.info.name,"\n"))
  cat(paste0("Groups info: ", getwd(), "/", group.info.name,"\n"))
  cat(paste0("Groups raster: ", getwd(), "/", group.name,"\n"))
  cat(paste0("Groups' mean raster: ", getwd(), "/", mean.name,"\n"))                             
  
  ###########################################################
  #                 Return absolute paths                   #
  ###########################################################
  
  out_paths <- list(
    cell.info.csv = normalizePath(cell.info.name, mustWork = FALSE),
    group.info.csv = normalizePath(group.info.name, mustWork = FALSE),
    group.raster = normalizePath(group.name, mustWork = FALSE),
    group.mean = normalizePath(mean.name, mustWork = FALSE)
  )
  
  invisible(out_paths)    
}
