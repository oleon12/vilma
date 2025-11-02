write.vilma.dist <- function(vilma.dist, file, raster.format = c("tif","grd","img"), overwrite = TRUE){
  
  if (class(vilma.dist) != "vilma.dist") {
    stop("Input is not a vilma.dist object. See points.to.raster().")
  }
  
  if (class(file) != "character") {
    stop("File name must be a character object.")
  }
  
  if (length(raster.format) > 1) {
    raster.format <- "tif"
    cat("\n")
    message("TIF format selected")
    cat("\n")
  }
  
  ###########################################################
  #                  Output file names                      #
  ###########################################################
  
  dist.name <- paste0(file, ".csv")
  r.name    <- paste0(file, "_richness.", raster.format)
  ab.name   <- paste0(file, "_abundance.", raster.format)
  log.name  <- paste0(file, "_log.txt")
  
  ###########################################################
  #                     Save files                          #
  ###########################################################
  
  write.csv(x = vilma.dist$distribution, file = dist.name, quote = FALSE, row.names = FALSE)
  suppressMessages(terra::writeRaster(x = vilma.dist$r.raster,  filename = r.name, overwrite = overwrite))
  suppressMessages(terra::writeRaster(x = vilma.dist$ab.raster, filename = ab.name, overwrite = overwrite))
  capture.output(print(vilma.dist), file = log.name)
  
  ###########################################################
  #                  Output messages                        #
  ###########################################################
  
  cat("Saved files:\n\n")
  cat(paste0("Distribution: ", getwd(), "/", dist.name, "\n"))
  cat(paste0("Richness raster: ", getwd(), "/", r.name, "\n"))
  cat(paste0("Abundance raster: ", getwd(), "/", ab.name, "\n"))
  cat(paste0("Summary log: ", getwd(), "/", log.name, "\n"))
  
  ###########################################################
  #                 Return absolute paths                   #
  ###########################################################
  
  out_paths <- list(
    csv       = normalizePath(dist.name, mustWork = FALSE),
    richness  = normalizePath(r.name, mustWork = FALSE),
    abundance = normalizePath(ab.name, mustWork = FALSE),
    log       = normalizePath(log.name, mustWork = FALSE)
  )
  
  invisible(out_paths)
}
