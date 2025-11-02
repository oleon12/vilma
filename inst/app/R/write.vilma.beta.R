#' Write a \code{vilma.beta} Object to Disk
#'
#' Exports a \code{vilma.beta} object to disk. The set of files written depends
#' on the beta-diversity algorithm used and stored in \code{vilma.beta$algorithm}:
#' \itemize{
#'   \item \strong{PhyloSor}: writes similarity and dissimilarity CSV matrices, a set of
#'         rasters (one per entry in \code{vilma.beta$rasters}), and a text summary via
#'         \code{print.vilma.beta()}.
#'   \item \strong{PhyloBeta}: writes total dissimilarity, turnover, and nestedness CSV tables,
#'         a set of rasters (one per entry in \code{vilma.beta$rasters}), and a text summary.
#'   \item \strong{Other algorithms} (e.g., UniFrac, Rao β): writes a single CSV named after the
#'         algorithm, a set of rasters (one per entry in \code{vilma.beta$rasters}), and a text summary.
#' }
#'
#' @param vilma.beta A \code{vilma.beta} object returned by one of the beta-diversity
#'   functions (e.g., \code{\link{phylo.beta}}, \code{\link{phylosor.calc}},
#'   \code{\link{rao.beta}}, \code{\link{unifrac.calc}}).
#' @param file Character string giving the prefix for output file names.
#'   The function appends appropriate extensions automatically.
#' @param raster.format Character; output format for raster layers.
#'   One of \code{"tif"}, \code{"grd"}, or \code{"img"}.
#'   Defaults to \code{"tif"} if multiple or invalid options are provided.
#' @param overwrite Logical; whether to overwrite existing raster files.
#'   Defaults to \code{TRUE}.
#'
#' @details
#' The following files are created in the current working directory:
#'
#' \strong{If \code{vilma.beta$algorithm == "PhyloSor"}:}
#' \itemize{
#'   \item \code{<file>_dissimilarity.csv} — dissimilarity matrix.
#'   \item \code{<file>_similarity.csv} — similarity matrix.
#'   \item \code{<file>_<raster_name>.<format>} — one raster per entry in \code{vilma.beta$rasters}.
#'   \item \code{<file>_log.txt} — textual summary from \code{print.vilma.beta()}.
#' }
#'
#' \strong{If \code{vilma.beta$algorithm == "PhyloBeta"}:}
#' \itemize{
#'   \item \code{<file>_dissimilarity.csv} — total dissimilarity matrix.
#'   \item \code{<file>_turnover.csv} — turnover component.
#'   \item \code{<file>_nestedness.csv} — nestedness component.
#'   \item \code{<file>_<raster_name>.<format>} — one raster per entry in \code{vilma.beta$rasters}.
#'   \item \code{<file>_log.txt} — textual summary from \code{print.vilma.beta()}.
#' }
#'
#' \strong{Otherwise} (e.g., UniFrac, Rao β):
#' \itemize{
#'   \item \code{<file>_<algorithm>.csv} — dissimilarity matrix for the algorithm.
#'   \item \code{<file>_<raster_name>.<format>} — one raster per entry in \code{vilma.beta$rasters}.
#'   \item \code{<file>_log.txt} — textual summary from \code{print.vilma.beta()}.
#' }
#'
#' Absolute paths to the written files are printed to the console on completion.
#'
#' @return
#' Invisibly returns \code{NULL}. Files are written to disk.
#'
#' @examples
#' \dontrun{
#' # PhyloSor output
#' write.vilma.beta(vilma_beta_phylosor, file = "results/phylosor")
#'
#' # PhyloBeta output
#' write.vilma.beta(vilma_beta_phylobeta, file = "results/phylobeta")
#'
#' # Other algorithm (e.g., UniFrac / Rao beta)
#' write.vilma.beta(vilma_beta_unifrac, file = "results/unifrac")
#' }
#'
#' @seealso
#' \code{\link{phylo.beta}}, \code{\link{phylosor.calc}}, \code{\link{rao.beta}},
#' \code{\link{unifrac.calc}}, \code{\link{print.vilma.beta}}, \code{\link{view.vilma.beta}}
#'
#' @export

write.vilma.beta <- function(vilma.beta, file, raster.format = c("tif","grd","img"), overwrite = TRUE){
  
  if (class(vilma.beta) != "vilma.beta") {
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
  
  
  if(vilma.beta$algorithm == "PhyloSor"){
    
    ###########################################################
    #                  Output file names                      #
    ###########################################################
    
    dis.name <- paste0(file, "_dissimilarity.csv")
    sim.name <- paste0(file, "_similarity.csv")
    log.name  <- paste0(file, "_log.txt")
    r.names <- list()
    
    for(i in seq_along(vilma.beta$rasters)){
      r.names[[i]] <- paste0(file,"_",names(vilma.beta$rasters)[i],".",raster.format)
    }
    
    names(r.names) <- names(vilma.beta$rasters)
    
    ###########################################################
    #                     Save files                          #
    ###########################################################
    
    write.csv(x = as.matrix(vilma.beta$similarity), file = sim.name, quote = FALSE, row.names = FALSE)
    write.csv(x = as.matrix(vilma.beta$dissimilarity), file = dis.name, quote = FALSE, row.names = FALSE)
    capture.output(print(vilma.beta), file = log.name)
    
    for(i in seq_along(r.names)){
      suppressMessages(terra::writeRaster(x = vilma.beta$rasters[[i]], filename = r.names[[i]], overwrite = overwrite))
    }
    
    ###########################################################
    #                  Output messages                        #
    ###########################################################
    
    cat("Saved files:\n\n")
    cat(paste0("Dissimilarity table: ", getwd(), "/", dis.name, "\n"))
    cat(paste0("Similarity table: ", getwd(), "/", sim.name, "\n"))
    
    for(i in seq_along(r.names)){
      cat(paste0("Raster: ", getwd(), "/", r.names[[i]],"\n"))
    }
    
    cat(paste0("Summary log: ", getwd(), "/", log.name, "\n"))
    
  }
  
  if(vilma.beta$algorithm == "PhyloBeta"){
    
    ###########################################################
    #                  Output file names                      #
    ###########################################################
    
    dis.name <- paste0(file, "_dissimilarity.csv")
    turn.name <- paste0(file, "_turnover.csv")
    nest.name <- paste0(file, "_nestedness.csv")
    log.name  <- paste0(file, "_log.txt")
    r.names <- list()
    
    for(i in seq_along(vilma.beta$rasters)){
      r.names[[i]] <- paste0(file,"_",names(vilma.beta$rasters)[i],".",raster.format)
    }
    
    ###########################################################
    #                     Save files                          #
    ###########################################################
    
    write.csv(x = as.matrix(vilma.beta$total.dissimilarity), file = dis.name, quote = FALSE, row.names = FALSE)
    write.csv(x = as.matrix(vilma.beta$turnover), file = turn.name, quote = FALSE, row.names = FALSE)
    write.csv(x = as.matrix(vilma.beta$nestedness), file = nest.name, quote = FALSE, row.names = FALSE)
    
    capture.output(print(vilma.beta), file = log.name)
    
    for(i in seq_along(r.names)){
      suppressMessages(terra::writeRaster(x = vilma.beta$rasters[[i]], filename = r.names[[i]], overwrite = overwrite))
    }
    
    ###########################################################
    #                  Output messages                        #
    ###########################################################
    
    cat("Saved files:\n\n")
    cat(paste0("Dissimilarity table: ", getwd(), "/", dis.name, "\n"))
    cat(paste0("Turnover table: ", getwd(), "/", turn.name, "\n"))
    cat(paste0("Nestedness table: ", getwd(), "/", nest.name, "\n"))
    
    
    for(i in seq_along(r.names)){
      cat(paste0("Raster: ", getwd(), "/", r.names[[i]],"\n"))
    }
    
    cat(paste0("Summary log: ", getwd(), "/", log.name, "\n"))
    
  }
  
  if(vilma.beta$algorithm != "PhyloBeta" && vilma.beta$algorithm != "PhyloSor"){
    
    ###########################################################
    #                  Output file names                      #
    ##########################################################
    
    dis.name <- paste0(file,"_",vilma.beta$algorithm, ".csv")
    log.name  <- paste0(file, "_log.txt")
    r.names <- list()
    
    for(i in seq_along(vilma.beta$rasters)){
      r.names[[i]] <- paste0(file,"_",names(vilma.beta$rasters)[i],".",raster.format)
    }
    
    ###########################################################
    #                     Save files                          #
    ###########################################################
    
    write.csv(x = as.matrix(vilma.beta[[2]]), file = dis.name, quote = FALSE, row.names = FALSE)
    capture.output(print(vilma.beta), file = log.name)
    
    for(i in seq_along(r.names)){
      suppressMessages(terra::writeRaster(x = vilma.beta$rasters[[i]], filename = r.names[[i]], overwrite = overwrite))
    }
    
    ###########################################################
    #                  Output messages                        #
    ###########################################################
    
    cat("Saved files:\n\n")
    cat(paste0(vilma.beta$algorithm," table: ", getwd(), "/", dis.name, "\n"))

    for(i in seq_along(r.names)){
      cat(paste0("Raster: ", getwd(), "/", r.names[[i]],"\n"))
    }
    
    cat(paste0("Summary log: ", getwd(), "/", log.name, "\n"))
    
  }
  
}
