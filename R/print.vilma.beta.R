#' Print summary for a \code{vilma.beta} object
#' @description
#' Displays a concise summary of the results from \code{\link{phylosor.calc}} objects,
#' including general dataset information, summary statistics for the PhyloSor
#' similarity and dissimilarity matrices, and ordination diagnostics (PCoA and NMDS).
#'
#' @param x An object of class \code{vilma.beta}, typically produced by
#'   \code{\link{phylosor.calc}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The print method provides an overview of the main components contained in a
#' \code{vilma.beta} object:
#' \itemize{
#'   \item The number of analyzed communities (cells) and pairwise comparisons.
#'   \item Summary statistics (mean +/- SD) for PhyloSor similarity and dissimilarity.
#'   \item The proportion of variance explained by the first two PCoA axes (based on
#'         the dissimilarity matrix).
#'   \item The NMDS stress value, indicating ordination fit quality.
#'   \item A list of available raster outputs (mean similarity/dissimilarity per cell,
#'         PCoA axes 1-2, and NMDS axes 1-2).
#' }
#'
#' This function is automatically invoked when printing a \code{vilma.beta} object
#' to the console.
#'
#' @return
#' Invisibly returns the input object \code{x}, allowing further use in pipelines.
#'
#' @seealso
#' \code{\link{phylosor.calc}}, \code{\link{faith.pd}}, \code{\link{pd.taxon}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' phs <- phylosor.calc(tree, dist)
#' print(phs)
#'
#' # The summary will display overall similarity, dissimilarity,
#' # and ordination statistics.
#' }
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @export
#' @method print vilma.beta

print.vilma.beta <- function(x, ...){
  # 1. Input validation
  if (!inherits(x, "vilma.beta")) {
    stop("Object must be of class 'vilma.beta'")
  }
  
  if (is.null(x$distribution)) {
    stop("'distribution' component is missing from this vilma.pd object")
  }
  
  if(x$algorithm == "PhyloSor"){
    #Extract information for cleaner code
  
    n_cells <- length(unique(x$distribution$Cell))
  
    sim_vals <- as.numeric(x$similarity)
    dis_vals <- as.numeric(x$dissimilarity)
  
    # 3. Formatted output
    cat("\n")
    cat("vilma.beta Object Summary")
    cat("\n")
    cat("-------------------------")
    cat("\n\n")
  
    cat(sprintf("Communities (cells): %d",n_cells))
    cat("\n")
    cat("Pairwise comparisons:", n_cells * (n_cells - 1) / 2, "\n")
  
    cat("\n")
  
    cat("Similarity:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(sim_vals, na.rm =T), sd(sim_vals, na.rm =T)))
    cat("\n")
  
    cat("\nDissimilarity:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(dis_vals, na.rm =T), sd(dis_vals, na.rm =T)))
    cat("\n")
  
    eig <- sum(x$pcoa.eig) * 100
    cat("\nPCoA (dissimilarity)\n")
    cat(sprintf("  - First two axes explain ~ %.2f%% of total variance\n", eig))
  
    cat("\n")
  
    cat("NMDS \n")
    cat(sprintf("  - Stress: %.3f", x$ndms.stress))
  
    cat("\n")
  
    cat("\nAvailable rasters:\n")
    cat("  - Mean similarity and dissimilarity per cell\n")
    cat("  - PCoA axes 1-2 and NMDS axes 1-2 (spatially mapped)\n")
    cat("\n")  
  
  }
  
  if(x$algorithm == "UniFrac"){
  
    #Extract information for cleaner code
    n_cells <- length(unique(x$distribution$Cell))
    dis_vals <- as.numeric(x$UniFrac)
    
    # 3. Formatted output
    cat("\n")
    cat("vilma.beta Object Summary")
    cat("\n")
    cat("-------------------------")
    cat("\n\n")
  
    cat(sprintf("Communities (cells): %d",n_cells))
    cat("\n")
    cat("Pairwise comparisons:", n_cells * (n_cells - 1) / 2, "\n")
    
    cat("\nUniFrac:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(dis_vals, na.rm =T), sd(dis_vals, na.rm =T)))
    cat("\n")
    
    eig <- sum(x$pcoa.eig) * 100
    cat("\nPCoA\n")
    cat(sprintf("  - First two axes explain ~ %.2f%% of total variance\n", eig))
  
    cat("\n")
  
    cat("NMDS \n")
    cat(sprintf("  - Stress: %.3f", x$ndms.stress))
  
    cat("\n")
  
    cat("\nAvailable rasters:\n")
    cat("  - Mean similarity and dissimilarity per cell\n")
    cat("  - PCoA axes 1-2 and NMDS axes 1-2 (spatially mapped)\n")
    cat("\n")
  
  }
  
  if(x$algorithm == "beta.MPD"){
  
    #Extract information for cleaner code
    n_cells <- length(unique(x$distribution$Cell))
    dis_vals <- as.numeric(x$bMPD)
    
    # 3. Formatted output
    cat("\n")
    cat("vilma.beta Object Summary")
    cat("\n")
    cat("-------------------------")
    cat("\n\n")
  
    cat(sprintf("Communities (cells): %d",n_cells))
    cat("\n")
    cat("Pairwise comparisons:", n_cells * (n_cells - 1) / 2, "\n")
    
    cat("\nbetaMPD:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(dis_vals, na.rm =T), sd(dis_vals, na.rm =T)))
    cat("\n")
    
    eig <- sum(x$pcoa.eig) * 100
    cat("\nPCoA\n")
    cat(sprintf("  - First two axes explain ~ %.2f%% of total variance\n", eig))
  
    cat("\n")
  
    cat("NMDS \n")
    cat(sprintf("  - Stress: %.3f", x$ndms.stress))
  
    cat("\n")
  
    cat("\nAvailable rasters:\n")
    cat("  - Mean similarity and dissimilarity per cell\n")
    cat("  - PCoA axes 1-2 and NMDS axes 1-2 (spatially mapped)\n")
    cat("\n")
  
  }
  
  if(x$algorithm == "beta.MNTD"){
  
    #Extract information for cleaner code
    n_cells <- length(unique(x$distribution$Cell))
    dis_vals <- as.numeric(x$bMNTD)
    
    # 3. Formatted output
    cat("\n")
    cat("vilma.beta Object Summary")
    cat("\n")
    cat("-------------------------")
    cat("\n\n")
  
    cat(sprintf("Communities (cells): %d",n_cells))
    cat("\n")
    cat("Pairwise comparisons:", n_cells * (n_cells - 1) / 2, "\n")
    
    cat("\nbetaMNTD:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(dis_vals, na.rm =T), sd(dis_vals, na.rm =T)))
    cat("\n")
    
    eig <- sum(x$pcoa.eig) * 100
    cat("\nPCoA\n")
    cat(sprintf("  - First two axes explain ~ %.2f%% of total variance\n", eig))
  
    cat("\n")
  
    cat("NMDS \n")
    cat(sprintf("  - Stress: %.3f", x$ndms.stress))
  
    cat("\n")
  
    cat("\nAvailable rasters:\n")
    cat("  - Mean similarity and dissimilarity per cell\n")
    cat("  - PCoA axes 1-2 and NMDS axes 1-2 (spatially mapped)\n")
    cat("\n")
  
  }
  
  if(x$algorithm == "PhyloBeta"){
    
    #Extract information for cleaner code
    n_cells <- length(unique(x$distribution$Cell))
    total <- as.numeric(x$total.dissmilarity)
    turnover <- as.numeric(x$turnover)
    nes <- as.numeric(x$nestedness)
    
    # 3. Formatted output
    cat("\n")
    cat("vilma.beta Object Summary")
    cat("\n")
    cat("-------------------------")
    cat("\n\n")
  
    cat(sprintf("Communities (cells): %d",n_cells))
    cat("\n")
    cat("Pairwise comparisons:", n_cells * (n_cells - 1) / 2, "\n")
    
    cat("\nPhylo Beta Diversity:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(total, na.rm =T), sd(total, na.rm =T)))
    cat("\n")
    
    cat("\nTurnover:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(turnover, na.rm =T), sd(turnover, na.rm =T)))
    cat("\n")
    
    cat("\nNestedness:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(nes, na.rm =T), sd(nes, na.rm =T)))
    cat("\n")
    
    eig <- apply(x$stats[,1:2], 1, sum)
    eig <- eig * 100
    
    cat("\nPCoA\n")
    cat(sprintf("  - Total: First two axes explain ~ %.2f%% of total variance\n", eig[1]))
    cat(sprintf("  - Turnover: First two axes explain ~ %.2f%% of total variance\n", eig[2]))
    cat(sprintf("  - Nestedness: First two axes explain ~ %.2f%% of total variance\n", eig[3]))
    cat("\n")
  
    cat("NMDS \n")
    cat(sprintf("  - Total Stress: %.3f\n", x$stats[1,3]))
    cat(sprintf("  - Turnover Stress: %.3f\n", x$stats[2,3]))
    cat(sprintf("  - Netedness Stress: %.3f\n", x$stats[3,3]))
  
    cat("\n")
  
    cat("\nAvailable rasters:\n")
    cat("  - Mean Total, Turnover, and Nestedness per cell\n")
    cat("  - PCoA axes 1-2 and NMDS axes 1-2 (spatially mapped)\n")
    cat("\n")
  
  }
  
  if(x$algorithm == "RaoBeta"){
    
    #Extract information for cleaner code
    n_cells <- length(unique(x$distribution$Cell))
    dis_vals <- as.numeric(x$Rao.Beta)
    
    # 3. Formatted output
    cat("\n")
    cat("vilma.beta Object Summary")
    cat("\n")
    cat("-------------------------")
    cat("\n\n")
  
    cat(sprintf("Communities (cells): %d",n_cells))
    cat("\n")
    cat("Pairwise comparisons:", n_cells * (n_cells - 1) / 2, "\n")
    
    cat("\nRao Beta Diversity:\n")
    cat(sprintf("  - Mean +/- SD: %.1f +/- %.1f" , 
              mean(dis_vals, na.rm =T), sd(dis_vals, na.rm =T)))
    cat("\n")
    
    eig <- sum(x$pcoa.eig) * 100
    cat("\nPCoA\n")
    cat(sprintf("  - First two axes explain ~ %.2f%% of total variance\n", eig))
  
    cat("\n")
  
    cat("NMDS \n")
    cat(sprintf("  - Stress: %.3f", x$ndms.stress))
  
    cat("\n")
  
    cat("\nAvailable rasters:\n")
    cat("  - Mean Rao Beta per cell\n")
    cat("  - PCoA axes 1-2 and NMDS axes 1-2 (spatially mapped)\n")
    cat("\n")
    
  }
  
}

