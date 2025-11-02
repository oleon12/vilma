#' Write a \code{vilma.dist} Object to Disk
#' @description
#' Exports a \code{vilma.dist} object to disk by writing its core outputs —
#' the species distribution table and the corresponding richness and abundance
#' raster layers — along with a text summary generated from \code{print.vilma.dist()}.
#' The function automatically names and saves these files based on the
#' user-provided \code{file} prefix.
#'
#' @param vilma.dist A \code{vilma.dist} object (see \code{\link{points_to_raster}}).
#' @param file Character string giving the prefix for output file names.
#'   The function appends appropriate extensions automatically.
#' @param raster.format Character; output format for raster layers.
#'   One of \code{"tif"}, \code{"grd"}, or \code{"img"}.
#'   Defaults to \code{"tif"} if multiple or invalid options are provided.
#' @param overwrite Logical; whether to overwrite existing files.
#'   Defaults to \code{TRUE}.
#'
#' @details
#' The function generates four files in the current working directory:
#' \itemize{
#'   \item \code{<file>.csv} — species-by-cell distribution table.
#'   \item \code{<file>_richness.<format>} — raster map of species richness.
#'   \item \code{<file>_abundance.<format>} — raster map of abundance.
#'   \item \code{<file>_log.txt} — textual summary produced by \code{print.vilma.dist()}.
#' }
#'
#' All files are written to the current working directory; absolute paths are
#' printed to the console upon completion.
#'
#' @return
#' Invisibly returns a named list with absolute file paths:
#' \itemize{
#'   \item \code{$csv} — path to the distribution CSV table.
#'   \item \code{$richness} — path to the richness raster file.
#'   \item \code{$abundance} — path to the abundance raster file.
#'   \item \code{$log} — path to the summary text file.
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{points_to_raster}}, \code{\link{print.vilma.dist}},
#' \code{\link{faith.pd}}, \code{\link{write.vilma.pd}}
#'
#' @examples
#' \dontrun{
#' data("example_vilma_dist")
#' write.vilma.dist(example_vilma_dist, file = "results/dist_export")
#' }
#'
#' @export

write.vilma.dist <- function(vilma.dist, file, raster.format = c("tif","grd","img"), overwrite = TRUE){
  
  if (!inherits(vilma.dist, "vilma.dist")) {
    stop("Input is not a vilma.dist object. See points_to_raster().")
  }
  
  if (!inherits(file, "character")) {
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
  suppressMessages(writeRaster(x = vilma.dist$r.raster,  filename = r.name, overwrite = overwrite))
  suppressMessages(writeRaster(x = vilma.dist$ab.raster, filename = ab.name, overwrite = overwrite))
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
