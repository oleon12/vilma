#' Net Relatedness Index (NRI)
#' @description
#' Computes the Net Relatedness Index (NRI) for communities (grid cells) based on 
#' phylogenetic structure. NRI is derived from the standardized effect size (SES) 
#' of the Mean Pairwise Distance (MPD) under a null model, with the sign inverted so 
#' that positive values indicate phylogenetic clustering and negative values indicate 
#' phylogenetic overdispersion.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo}, with branch lengths.
#' @param dist An object of class \code{vilma.dist}, typically produced by 
#'   \code{points.to.raster()}, containing species distribution data.
#' @param mpd.method Character string specifying how to handle single-species cells in MPD. 
#'   Options are:
#'   \itemize{
#'     \item \code{"root"} – For cells with single species, uses distance from root to taxon.
#'     \item \code{"node"} – For cells with single species, uses branch length to nearest ancestral node.
#'     \item \code{"exclude"} – Excludes single-species cells from calculations (default).
#'   }
#' @param abundance Logical. If \code{TRUE}, weights MPD by relative species abundances 
#'   within each cell. Default is \code{FALSE}.
#' @param iterations Integer. Number of randomizations for the null model (default = 999).
#' @param sampling Character. Randomization strategy for null model:
#'   \itemize{
#'     \item \code{"taxa.label"} – Permutes species labels on the tree.
#'     \item \code{"range"} – Randomizes species ranges using swap algorithms.
#'     \item \code{"neighbor"} – Swaps species occurrences between adjacent cells.
#'     \item \code{"regional"} – Assembles communities by drawing from a regional pool.
#'   }
#' @param n.directions Character. Neighborhood adjacency definition for 
#'   \code{"neighbor"} method. Options: \code{"rook"}, \code{"bishop"}, or \code{"queen"}
#'   (default = \code{"queen"}).
#' @param regional.weight Character. Weighting scheme for species in the regional null model:
#'   \itemize{
#'     \item \code{"uniform"} – All species have equal probability.
#'     \item \code{"frequency"} – Weighted by the number of records per species.
#'     \item \code{"range"} – Weighted by the number of unique cells occupied.
#'   }
#'
#' @details
#' The Net Relatedness Index (NRI) is defined as:
#' \deqn{NRI = - \frac{MPD_{obs} - \mu(MPD_{null})}{\sigma(MPD_{null})}}
#' Positive NRI values indicate that species within communities are more 
#' closely related (phylogenetic clustering) than expected under the null model. 
#' Negative NRI values indicate that species are less related (phylogenetic overdispersion).
#'
#' NRI is computed at the \code{"cell"} level, returning results per grid cell, 
#' including a raster of NRI values for visualization.
#'
#' @return An object of class \code{vilma.pd}, which is a list containing:
#' \itemize{
#'   \item \code{distribution} – Original species distribution.
#'   \item \code{grid} – The reference raster grid.
#'   \item \code{pd.table} – Data frame with cells, species richness, and NRI values.
#'   \item \code{rasters} – List with abundance, richness, and NRI rasters.
#'   \item \code{calculation.method} – Method used for single-species cells in MPD.
#'   \item \code{abundance} – Logical, whether abundance-weighting was applied.
#'   \item \code{index} – The string \code{"nri.calc"}.
#' }
#'
#' @references
#' Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002). 
#' Phylogenies and community ecology. \emph{Annual Review of Ecology and Systematics}, 
#' 33, 475–505. \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' Gotelli, N. J., & Graves, G. R. (1996). \emph{Null Models in Ecology}. 
#' Smithsonian Institution Press.
#'
#' Gotelli, N. J. (2000). Null model analysis of species co‐occurrence patterns. 
#' \emph{Ecology}, 81(9), 2606–2621. \doi{10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2}
#'
#' @author 
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{\link{mpd.calc}}, \code{\link{mpd.calc.null}}, \code{\link{faith.pd}}, \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' # Calculate NRI using default null model (taxa.label)
#' nri_result <- nri.calc(tree = tree, dist = dist, iterations = 999)
#' 
#' print(nri_result)
#' plot(nri_result)
#' view.vilma(nri_result)
#'
#' # Calculate NRI using neighbor null model
#' nri_neighbor <- nri.calc(tree = tree, dist = dist,
#'                          iterations = 999, sampling = "neighbor",
#'                          n.directions = "queen")
#' }
#'
#'
#'
#' @export

nri.calc <- function(tree, dist,
                     mpd.method = c("root","node","exclude"),
                     abundance = FALSE, 
                     iterations = 999,
                     sampling = c("taxa.label","range","neigbor","regional"),
                     n.directions = c("rook","bishop","queen"),
                     regional.weight = c("uniform","frequency","range")){
  
  ##############################################################
  #                        VERIFICATION                        #
  ##############################################################
  
  if(!inherits(tree, "phylo")){
    stop("Input 'tree' must be an objecto of class 'phylo'.")
  }
  
  if (is.null(tree$edge.length)) {
    stop("The phylogenetic tree must have branch lengths.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(sampling) != 1){
    sampling <- "taxa.label"
  }
  
  if(length(n.directions) != 1){
    n.directions <- "queen"
  }
  
  if(length(regional.weight) != 1){
    regional.weight <- "uniform"
  }
  
  if(length(mpd.method)>2){
    method <- "exclude"
    message("Using method 'exclude' \n")
  }
  
  ##############################################################
  #                      END VERIFICATION                      #
  ##############################################################
  
  #####################################################################
  
  #########################################################################
  
  mpd1 <- mpd.calc(tree = tree, dist = dist, method = mpd.method, abundance = abundance)
  
  mpd.null <- suppressMessages(mpd.calc.null(pd = mpd1, tree = tree, dist = dist, iterations = iterations,
                                             method = "cell", sampling = sampling, n.directions = n.directions,
                                             regional.weight = regional.weight))
  
  nri.vals <- (-1 * values(mpd.null$Raster))
  mpd1$pd.table$PD <- nri.vals[!is.na(nri.vals)]
  
  grid0 <- dist$grid   # instead of mpd1$rasters$pd.raster
  grid0[] <- NA
  grid0[as.numeric(mpd1$pd.table$Cell)] <- mpd1$pd.table$PD

  structure(
    list(
      distribution = dist$distribution,
      grid = dist$grid,
      pd.table = mpd1$pd.table,
      rasters = list(ab.raster = dist$ab.raster,
                     r.raster = dist$r.raster,
                     pd.raster = grid0),
      calculation.method = mpd.method,
      abundance = abundance,
      index = "nri.calc"
    ),
    class = "vilma.pd"
  )
  
  
}
