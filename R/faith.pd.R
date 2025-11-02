#' Faith's Phylogenetic Diversity (PD)
#' @description
#' Calculates Faith’s Phylogenetic Diversity (PD) for species assemblages across spatial cells
#' using a rooted phylogenetic tree and species distributions. 
#' The function supports different treatments of single-species cells through the 
#' \code{method} argument.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist}, typically produced by 
#'   \code{points.to.raster()}, containing the species distribution data.
#' @param method Character string specifying how to handle single-species cells.
#'   Options are:
#'   \itemize{
#'     \item \code{"root"} – For cells with single species, PD is calculated considering the root path to the species.
#'     \item \code{"node"} – For cells with single species, PD considers the closest ancestral node length.
#'     \item \code{"exclude"} – Single-species cells are excluded from the calculation (default).
#'   }
#'
#' @details
#' Faith’s PD is calculated as the sum of the branch lengths of the minimum spanning
#' path connecting all species in a cell. Species not shared between \code{tree} and
#' \code{dist} are ignored, and informative messages are returned about missing species.
#'
#' @return An object of class \code{vilma.pd}, which is a list containing:
#' \itemize{
#'   \item \code{distribution} – Original species distribution.
#'   \item \code{grid} – The reference raster grid.
#'   \item \code{pd.table} – Data frame with cells, species richness, and PD values.
#'   \item \code{rasters} – List with abundance, richness, and PD rasters.
#'   \item \code{calculation.method} – Method used for single-species cells.
#'   \item \code{index} – The string \code{"faith.pd"}.
#' }
#'
#' @references
#' Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity. 
#' \emph{Biological Conservation}, 61(1), 1–10. 
#' \doi{10.1016/0006-3207(92)91201-3}
#'
#' Webb, C.O., Ackerly, D.D., McPeek, M.A., & Donoghue, M.J. (2002). 
#' Phylogenies and community ecology. 
#' \emph{Annual Review of Ecology and Systematics}, 33, 475–505. 
#' \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' # Calculate PD excluding single-species cells
#' fpd <- faith.pd(tree = tree, dist = dist, method = "exclude")
#' 
#' print(fpd)
#' plot(fpd)
#' view.vilma(fpd)
#' }
#' @export

faith.pd <- function(tree, dist, method = c("root", "node", "exclude")) {
  
  ##############################################################
  #                        VERIFICATION                        #
  ##############################################################
  
  if (!inherits(tree, "phylo")) {
    stop("Input 'tree' must be an object of class 'phylo'.")
  }
  
  if (is.null(tree$edge.length)) {
    stop("The phylogenetic tree must have branch lengths.")
  }
  
  if (!is.rooted(tree)) {
    stop("Faith's PD requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(method)>2){
    method <- "exclude"
    message("Using method 'exclude' \n")
  }
  
  
  ##############################################################
  #                     DATA PREPARATION                       #
  ##############################################################
  
  ##### Dist Processing #####
  distM <- dist$distribution
  
  # 3. Find and use only species that match between tree and dist
  tree.species <- tree$tip.label
  dist.species <- unique(distM$Sp)
  
  missing.in.tree <- setdiff(dist.species, tree.species)
  missing.in.dist <- setdiff(tree.species, dist.species)
  
  if (length(missing.in.tree) > 0) {
    message("Note: ", length(missing.in.tree), 
            " species from 'dist' not found in 'tree': ", 
            paste(head(missing.in.tree, 5), collapse = ", "),
            ifelse(length(missing.in.tree) > 5, "...", ""))
  }
  
  if (length(missing.in.dist) > 0) {
    message("Note: ", length(missing.in.dist),
            " species from 'tree' not found in 'dist': ",
            paste(head(missing.in.dist, 5), collapse = ", "),
            ifelse(length(missing.in.dist) > 5, "...", ""))
  }
  
  
  ############################################################### 
 
  # Use only intersecting species
  common.species <- intersect(tree.species, dist.species)
  
  if (length(common.species) == 0) {
    stop("No species in common between the tree and distribution data.")
  }
  
  message("Using ", length(common.species), " species in common between tree and distribution.")
  distM <- distM[distM$Sp %in% common.species, ]
  
  ###############################################################

  
  #Get unique cells and calculte Richness for each on
  
  presabs <- table(distM$Cell, distM$Sp)
  presabs <- (presabs > 0) *1
  presabs <- as.matrix(presabs)
  presabs <- as.matrix(apply(presabs, 1, sum))
  
  Cell <- rownames(presabs)
  SR <- as.vector(presabs)
  
  pd.res <- data.frame(Cell = Cell,
                       SR = SR)
                       
  ##############################################################
  #                      PD CALCULATION                        #
  ##############################################################
  
  # Handle single-species cells based on method
  if (method == "exclude") {
    message("Excluding ", sum(pd.res$SR == 1), " single-species cells.")
    pd.res <- pd.res[pd.res$SR > 1, ]
  }
  
  PD <- c()  # Pre-allocate for efficiency
  
  # Precompute root distances once
  depths <- node.depth.edgelength(tree)
  
  set.seed(123)
  
  for(i in 1:length(rownames(pd.res))) {
    
    cell.id <- pd.res$Cell[i]
    sp.list <- unique(distM$Sp[distM$Cell == cell.id])
    
    if(length(sp.list) != 1){
      # Multiple species: calculate full PD
      tree.subset <- drop.tip(tree, 
                                   tip = tree$tip.label[!tree$tip.label %in% sp.list])
      PDcore <- sum(tree.subset$edge.length)
      
      #To be consistent with method root for single species cells
      if(method == "root"){
        mrca_node <- getMRCA(tree, sp.list)
        root_to_mrca <- depths[mrca_node]
        PD <- c(PD,(PDcore + root_to_mrca))
      }else{
        PD <- c(PD, PDcore)
      }
      
    }else{
      PD <- c(PD, pd.taxon(tree, sp.list, method))
    }
}
  
  pd.res$PD <- PD
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################

  grid0 <- dist$grid
  
  values(grid0) <- NA
  
  #### Raster ####
  
  set.values(grid0, as.numeric(as.character(pd.res$Cell)),  as.numeric(pd.res$PD))
  
  
  ##############################################################
  #                        VILMA OBJECT                        #
  ##############################################################
  
  
  structure(
    list(
      distribution = dist$distribution,
      grid = dist$grid,
      pd.table = pd.res,
      rasters = list(
        ab.raster = dist$ab.raster,
        r.raster = dist$r.raster,
        pd.raster = grid0
      ),
      calculation.method = method,
      index = "faith.pd"
    ),
    class = "vilma.pd"
  )
}
