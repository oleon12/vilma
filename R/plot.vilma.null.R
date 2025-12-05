#' Plot Results of a VILMA Null Model
#' @description
#' Provides a quick visualization of null model results from \code{vilma.null} objects.  
#' Depending on the null model method, it plots either:
#' \itemize{
#'   \item \strong{Global} – A histogram of randomized mean PD values with the observed value indicated.
#'   \item \strong{Cell} – A raster map of SES values, with numerical SES displayed per cell.
#' }
#'
#' @param x An object of class \code{vilma.null}, typically the output of a null model function.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @details
#' This function is designed for quick inspection of null model outputs, 
#' either summarizing global PD deviations or mapping SES values per cell.  
#' For more advanced plotting and customization, users may directly access 
#' the raster or data stored in the \code{vilma.null} object.
#'
#' @return Produces a plot (histogram or raster map) depending on the method used 
#' to generate the null model. No values are returned.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
#' @method plot vilma.null

plot.vilma.null <- function(x, ...){
  
  if (!inherits(x, "vilma.null")) {
    stop("Input must be an object of class 'vilma.null'.")
  }
  
  if (identical(x$Method, "global")) {
    z <- suppressWarnings(as.numeric(x$null.pd))
    z <- z[is.finite(z)]
    if (length(z) == 0L) stop("No finite null values to plot for global method.")
    obs <- suppressWarnings(as.numeric(x$pd.obs))
    
    graphics::hist(z,
                   main = "Global Null Model: Mean PD",
                   xlab = "Mean PD",
                   col = "lightgray", border = "white")
    if (length(obs) == 1L && is.finite(obs)) {
      graphics::abline(v = obs, col = "red", lwd = 2)
      graphics::legend("topright", legend = c("Observed mean PD"), col = "red", lwd = 2, bty = "n")
    }
    return(invisible())
  }
  
  if (identical(x$Method, "cell")) {
    if (is.null(x$Raster)) stop("Missing SES raster in 'vilma.null' (x$Raster).")
    terra::plot(x$Raster, main = "Null model: SES per Cell")
    
    # Safe cell labeling: use SES values and valid cell ids
    #if (!is.null(x$CellValues) && all(c("Cell","SES") %in% names(x$CellValues))) {
    #  cells <- suppressWarnings(as.integer(x$CellValues$Cell))
    #  ses   <- suppressWarnings(as.numeric(x$CellValues$SES))
    #  ok    <- is.finite(ses) & is.finite(cells) &
    #           cells >= 1L & cells <= terra::ncell(x$Raster)
    #  if (any(ok)) {
    #    crds <- try(terra::xyFromCell(x$Raster, cells[ok]), silent = TRUE)
    #    if (!inherits(crds, "try-error")) {
    #      graphics::text(crds[,1], crds[,2], labels = round(ses[ok], 2),
    #                     cex = 0.6, col = "white")
    #    }
    #  }
    #}
    return(invisible())
  }
  
  stop("Unknown 'Method' in vilma.null. Expected 'global' or 'cell'.")
}

