#' Write a \code{vilma.null} Object to Disk
#'
#' Exports a \code{vilma.null} object to disk. For cell-based null models
#' (\code{Method == "cell"}), the function writes a CSV table of cell values,
#' a SES raster, and a text summary produced by \code{print.vilma.null()}.
#' For global null models (\code{Method != "cell"}), the function saves a PNG
#' histogram via \code{plot(vilma.null)} and a text summary.
#'
#' @param vilma.null A \code{vilma.null} object returned by one of the null
#'   model functions (e.g., \code{\link{faith.pd.null}}, \code{\link{mpd.calc.null}},
#'   \code{\link{mntd.calc.null}}, \code{\link{pe.calc.null}}, \code{\link{rao.calc.null}}).
#' @param file Character string giving the prefix for output file names.
#'   The function appends appropriate extensions automatically.
#' @param raster.format Character; output format for raster layers (cell mode).
#'   One of \code{"tif"}, \code{"grd"}, or \code{"img"}.
#'   Defaults to \code{"tif"} if multiple or invalid options are provided.
#' @param overwrite Logical; whether to overwrite existing files (cell mode rasters).
#'   Defaults to \code{TRUE}.
#'
#' @details
#' Output depends on \code{vilma.null$Method}:
#'
#' \strong{Cell method} (\code{Method == "cell"}):
#' \itemize{
#'   \item \code{<file>.csv} — cell-wise null results (e.g., SES, p-values), from \code{vilma.null$CellValues}.
#'   \item \code{<file>_ses_raster.<format>} — SES raster from \code{vilma.null$Raster}.
#'   \item \code{<file>_log.txt} — textual summary via \code{print.vilma.null()}.
#' }
#'
#' \strong{Global method} (\code{Method != "cell"}):
#' \itemize{
#'   \item \code{<file>null_hist.png} — histogram saved from \code{plot(vilma.null)}.
#'   \item \code{<file>_log.txt} — textual summary via \code{print.vilma.null()}.
#' }
#'
#' All files are written to the current working directory, and absolute paths
#' are displayed in the console upon completion.
#'
#' @return
#' Invisibly returns a named list of absolute file paths:
#' \itemize{
#'   \item \emph{Cell method}: \code{$pd.table}, \code{$pd.raster}, \code{$log}.
#'   \item \emph{Global method}: \code{$png.img}, \code{$log}.
#' }
#'
#' @examples
#' \dontrun{
#' # For a cell-based null result:
#' write.vilma.null(vilma_null_cell, file = "results/mpd_null_cell")
#'
#' # For a global null result:
#' write.vilma.null(vilma_null_global, file = "results/mpd_null_global")
#' }
#'
#' @seealso
#' \code{\link{faith.pd.null}}, \code{\link{mpd.calc.null}}, \code{\link{mntd.calc.null}},
#' \code{\link{pe.calc.null}}, \code{\link{rao.calc.null}},
#' \code{\link{print.vilma.null}}, \code{\link{plot.vilma.null}}
#'
#' @export

write.vilma.null <- function(vilma.null, file, raster.format = c("tif","grd","img"), overwrite = TRUE) {
  
  if (class(vilma.null) != "vilma.null") {
    stop("Input is not a vilma.null object.")
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
  
  
  if(vilma.null$Method == "cell"){
    ###########################################################
    #                  Output file names                      #
    ###########################################################
    
    pd.name   <- paste0(file, ".csv")
    rpd.name  <- paste0(file, "_ses_raster.", raster.format)
    log.name  <- paste0(file, "_log.txt")
    
    ###########################################################
    #                     Save files                          #
    ###########################################################
    
    write.csv(x = vilma.null$CellValues, file = pd.name, quote = FALSE, row.names = FALSE)
    suppressMessages(terra::writeRaster(x = vilma.null$Raster, filename = rpd.name, overwrite = overwrite))
    capture.output(print(vilma.null), file = log.name)
    
    ###########################################################
    #                  Output messages                        #
    ###########################################################
    
    cat("Saved files:\n\n")
    cat(paste0("Null table: ", getwd(), "/", pd.name, "\n"))
    cat(paste0("SES raster: ", getwd(), "/", rpd.name, "\n"))
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
  }else{
    
    ###########################################################
    #                  Output file names                      #
    ###########################################################
    
    png.name <- paste0(file, "null_hist.png")
    log.name <- paste0(file, "_log.txt")
    
    ###########################################################
    #                     Save files                          #
    ###########################################################
    
    png(filename = png.name)
      plot(vilma.null)
    dev.off()
    
    capture.output(print(vilma.null), file = log.name)
    
    ###########################################################
    #                  Output messages                        #
    ###########################################################
    
    cat("Saved files:\n\n")
    cat(paste0("Histogram : ", getwd(), "/", png.name, "\n"))
    cat(paste0("Summary log: ", getwd(), "/", log.name, "\n"))
    
    ###########################################################
    #                 Return absolute paths                   #
    ###########################################################
    
    out_paths <- list(
      png.img  = normalizePath(png.name, mustWork = FALSE),
      log       = normalizePath(log.name, mustWork = FALSE)
    )
    
    invisible(out_paths)
    
  }
  
  
  
}
