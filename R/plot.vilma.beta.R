#' Plot method for \code{vilma.beta} objects
#' @description
#' Produces quick-look maps of mean PhyloSor similarity and mean PhyloSor
#' dissimilarity per cell from a \code{vilma.beta} object (output of
#' \code{\link{phylosor.calc}}). Values are drawn using the underlying
#' \code{terra} raster and annotated with the cell-wise means.
#'
#' @param x An object of class \code{vilma.beta}, typically returned by
#'   \code{\link{phylosor.calc}} and containing the rasters
#'   \code{$rasters$mean.similarity} and \code{$rasters$mean.dissimilarity}.
#' @param ... Additional graphical arguments (currently ignored).
#'
#' @details
#' The function opens a 2-panel plotting layout (\code{par(mfrow=c(2,1))})
#' and renders:
#' \enumerate{
#'   \item Mean PhyloSor similarity per cell (\code{$rasters$mean.similarity}).
#'   \item Mean PhyloSor dissimilarity per cell (\code{$rasters$mean.dissimilarity}).
#' }
#' For each panel, numeric values are overlaid at cell centroids (rounded to 2
#' decimals). The plotting is performed with \pkg{terra} (i.e., \code{terra::plot}).
#'
#' Expect \code{NA} cells when communities were filtered (e.g., \code{method = "exclude"})
#' or when no value could be computed. The function temporarily modifies
#' graphical parameters and restores them on exit.
#'
#' @return
#' Invisibly returns \code{x} (the input object), allowing further use in pipelines.
#'
#' @seealso
#' \code{\link{phylosor.calc}}, \code{\link{print.vilma.beta}}, \code{\link{faith.pd}}
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
#' @method plot vilma.beta

plot.vilma.beta <- function(x, ...){
  
  if (!inherits(x, "vilma.beta")) {
    stop("Input must be an object of class 'vilma.beta'.")
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # small helper to safely annotate numbers on a SpatRaster
  .label_cells <- function(r){
    vals  <- try(terra::values(r), silent = TRUE)
    if (inherits(vals, "try-error") || length(vals) == 0) return(invisible())
    # finite indices
    idx   <- which(is.finite(vals))
    if (length(idx) == 0) return(invisible())
    crds  <- try(terra::xyFromCell(r, idx), silent = TRUE)
    if (inherits(crds, "try-error")) return(invisible())
    graphics::text(crds[,1], crds[,2], labels = round(vals[idx], 2),
                   cex = 0.6, col = "white")
    invisible()
  }
  
  if (x$algorithm == "PhyloSor"){
    
    have_sim  <- !is.null(x$rasters$mean.similarity)
    have_diss <- !is.null(x$rasters$mean.dissimilarity)
    
    if (have_sim && have_diss) {
      par(mfrow = c(2, 1), mar = c(3, 3, 3, 5))
    } else {
      par(mfrow = c(1, 1), mar = c(3, 3, 3, 5))
    }
    
    if (have_sim) {
      terra::plot(x$rasters$mean.similarity, main = "Mean similarity per cell")
      #.label_cells(x$rasters$mean.similarity)
    }
    
    if (have_diss) {
      terra::plot(x$rasters$mean.dissimilarity, main = "Mean dissimilarity per cell")
      #.label_cells(x$rasters$mean.dissimilarity)
    }
  }
  
  if (x$algorithm == "UniFrac"){
    if (!is.null(x$rasters$mean.unifrac)) {
      par(mar = c(3, 3, 3, 5))
      terra::plot(x$rasters$mean.unifrac, main = "Mean UniFrac per cell")
      #.label_cells(x$rasters$mean.unifrac)
    }
  }
  
  if (x$algorithm == "beta.MPD"){
    if (!is.null(x$rasters$mean.bMPD)) {
      par(mar = c(3, 3, 3, 5))
      terra::plot(x$rasters$mean.bMPD, main = "Mean beta-MPD per cell")
      #.label_cells(x$rasters$mean.bMPD)
    }
  }
  
  if (x$algorithm == "beta.MNTD"){
    if (!is.null(x$rasters$mean.bMNTD)) {
      par(mar = c(3, 3, 3, 5))
      terra::plot(x$rasters$mean.bMNTD, main = "Mean beta-MNTD per cell")
      #.label_cells(x$rasters$mean.bMNTD)
    }
  }
  
  if (x$algorithm == "PhyloBeta"){
    if (!is.null(x$rasters$mean.total)) {
      par(mar = c(3, 3, 3, 5))
      terra::plot(x$rasters$mean.total, main = "Mean Phylogenetic beta-Diversity per cell")
      #.label_cells(x$rasters$mean.total)
    }
  }
  
  if (x$algorithm == "RaoBeta"){
    if (!is.null(x$rasters$mean.dissimilarity)) {
      par(mar = c(3, 3, 3, 5))
      terra::plot(x$rasters$mean.dissimilarity, main = "Mean Rao beta-Phylogenetic Diversity per cell")
      #.label_cells(x$rasters$mean.dissimilarity)
    }
  }
  
  invisible(x)
}

