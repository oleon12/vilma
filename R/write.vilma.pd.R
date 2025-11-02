#' Write a \code{vilma.pd} Object to Disk
#' @description
#' Exports a \code{vilma.pd} object to disk by writing its core outputs —
#' the phylogenetic diversity table and corresponding raster layer —
#' along with a text summary file generated from \code{print.vilma.pd()}.
#' The function automatically names and saves these files based on the
#' user-provided \code{file} prefix.
#'
#' @param vilma.pd A \code{vilma.pd} object, typically created using one of
#'   the alpha-diversity functions such as \code{\link{faith.pd}},
#'   \code{\link{mpd.calc}}, or \code{\link{mntd.calc}}.
#' @param file Character string giving the prefix for output file names.
#'   The function appends appropriate extensions automatically.
#' @param raster.format Character; output format for raster layers.
#'   One of \code{"tif"}, \code{"grd"}, or \code{"img"}.
#'   Defaults to \code{"tif"} if multiple or invalid options are provided.
#' @param overwrite Logical; whether to overwrite existing files.
#'   Defaults to \code{TRUE}.
#'
#' @details
#' The function generates three files in the current working directory:
#' \itemize{
#'   \item \code{<file>.csv} — cell-wise phylogenetic diversity values.
#'   \item \code{<file>_pd_raster.<format>} — raster map of PD values.
#'   \item \code{<file>_log.txt} — textual summary produced by \code{print.vilma.pd()}.
#' }
#'
#' All files are written to the current working directory, and absolute paths
#' are displayed in the console upon completion.
#'
#' @return
#' Invisibly returns a named list containing the absolute file paths of all
#' written outputs:
#' \itemize{
#'   \item \code{$pd.table} — path to the CSV table of PD values.
#'   \item \code{$pd.raster} — path to the raster file of PD.
#'   \item \code{$log} — path to the summary text file.
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{faith.pd}}, \code{\link{mpd.calc}}, \code{\link{mntd.calc}},
#' \code{\link{print.vilma.pd}}
#'
#' @examples
#' \dontrun{
#' data("example_vilma_pd")
#' write.vilma.pd(example_vilma_pd, file = "results/pd_export")
#' }
#'
#'
#' @export

write.vilma.pd <- function(vilma.pd, file, raster.format = c("tif","grd","img"), overwrite = TRUE) {
  
  if (!inherits(vilma.pd, "vilma.pd")) {
    stop("Input is not a vilma.pd object.")
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
  
  pd.name   <- paste0(file, ".csv")
  rpd.name  <- paste0(file, "_pd_raster.", raster.format)
  log.name  <- paste0(file, "_log.txt")
  
  ###########################################################
  #                     Save files                          #
  ###########################################################
  
  write.csv(x = vilma.pd$pd.table, file = pd.name, quote = FALSE, row.names = FALSE)
  suppressMessages(writeRaster(x = vilma.pd$rasters$pd.raster, filename = rpd.name, overwrite = overwrite))
  capture.output(print(vilma.pd), file = log.name)
  
  ###########################################################
  #                  Output messages                        #
  ###########################################################
  
  cat("Saved files:\n\n")
  cat(paste0("PD table: ", getwd(), "/", pd.name, "\n"))
  cat(paste0("PD raster: ", getwd(), "/", rpd.name, "\n"))
  cat(paste0("Summary log: ", getwd(), "/", log.name, "\n"))
  
  ###########################################################
  #                 Return absolute paths                   #
  ###########################################################
  
  out_paths <- list(
    pd.table  = normalizePath(pd.name, mustWork = FALSE),
    pd.raster = normalizePath(rpd.name, mustWork = FALSE),
    log       = normalizePath(log.name, mustWork = FALSE)
  )
  
  invisible(out_paths)
}
