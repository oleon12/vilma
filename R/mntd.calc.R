#' Mean Nearest Taxon Distance (MNTD) Phylogenetic Diversity Index
#' @description
#' Computes the Mean Nearest Taxon Distance (MNTD) among species within spatial cells
#' based on a phylogenetic tree and species distribution data.  
#' The function allows different treatments of single-species cells and can optionally
#' incorporate species abundances.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#'   For methods \code{"root"} and \code{"node"}, the tree must be rooted.
#' @param dist An object of class \code{vilma.dist}, typically produced by 
#'   \code{points.to.raster()}, containing the species distribution data.
#' @param method Character string specifying how to handle single-species cells. Options are:
#'   \itemize{
#'     \item \code{"root"} – For cells with single species, uses distance from root to taxon for single-species cells.
#'     \item \code{"node"} – For cells with single species, uses branch length to the nearest ancestral node.
#'     \item \code{"exclude"} – Excludes single-species cells from calculations (default).
#'   }
#' @param abundance Logical. If \code{TRUE}, weights pairwise distances by relative species abundances 
#'   within each cell. Default is \code{FALSE}.
#'
#' @details
#' MNTD is calculated as the average phylogenetic distance from each species in a community (cell) 
#' to its closest relative (nearest neighbor) within the same community. 
#' When \code{abundance = TRUE}, MNTD is weighted by the relative abundance of each species, 
#' such that dominant taxa contribute more to the index than rare taxa.
#'
#' If a cell contains only one species, treatment depends on the \code{method} argument.
#'
#' @return An object of class \code{vilma.pd}, which is a list containing:
#' \itemize{
#'   \item \code{distribution} – Original species distribution.
#'   \item \code{grid} – The reference raster grid.
#'   \item \code{pd.table} – Data frame with cells, species richness, and MPD values.
#'   \item \code{rasters} – List with abundance, richness, and MNTD rasters.
#'   \item \code{calculation.method} – Method used for single-species cells.
#'   \item \code{abundance} – Logical, whether abundance-weighting was applied.
#'   \item \code{index} – The string \code{"mntd.calc"}.
#' }
#'
#' @references
#' Webb, C.O., Ackerly, D.D., McPeek, M.A., & Donoghue, M.J. (2002). 
#' Phylogenies and community ecology. 
#' \emph{Annual Review of Ecology and Systematics}, 33, 475–505. 
#' \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' Hardy, O.J. (2008). Testing the spatial phylogenetic structure of local communities: 
#' statistical performances of different null models and test statistics on a locally neutral community. 
#' \emph{Journal of Ecology}, 96(5), 914–926. 
#' \doi{10.1111/j.1365-2745.2008.01421.x}
#'
#' Swenson, N.G. (2014). \emph{Functional and phylogenetic ecology in R}. 
#' Springer, New York. 
#' \doi{10.1007/978-1-4614-9542-0}
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com/}
#'
#' @seealso \code{\link{mpd.calc}}, \code{\link{mntd.calc.null}}, \code{\link{points_to_raster}}
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' # Calculate abundance-unweighted MNTD
#' result <- mntd.calc(tree = tree, dist = dist, method = "exclude")
#' 
#' print(result)
#' plot(result)
#' view.vilma(result)
#'
#' # Calculate abundance-weighted MNTD
#' result_ab <- mntd.calc(tree = tree, dist = dist, 
#'                       method = "exclude", abundance = TRUE)
#' }
#' @export

mntd.calc <- function(tree, dist, method = c("root","node","exclude"), abundance = FALSE){
  ##############################################################
  #                        VERIFICATION                        #
  ##############################################################
  
  if(!inherits(tree, "phylo")){
    stop("Input 'tree' must be an object of class 'phylo'.")	
  }
  
  if(is.null(tree$edge.length)){
    stop("The phylogenetic tree must have branch lengths.")
  }
  
  if(!is.rooted(tree)){
    if(method %in% c("node","root")){
      stop("MPD with methods node/root requires a rooted tree.")
    }	
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
  
  distM <- dist$distribution
  
  tree.species <- tree$tip.label
  dist.species <- unique(distM$Sp)
  
  missing.in.tree <- setdiff(dist.species, tree.species)
  missing.in.dist <- setdiff(tree.species, dist.species)
  
  if(length(missing.in.tree) > 0){
    message("Note: ", length(missing.in.tree),
            "species from 'dist' not found in 'tree': ",
            paste(head(missing.in.tree, 5), collapse = ","),
            ifelse(length(missing.in.tree) > 5, "..." , "") )
  }
  
  common.species <- intersect(tree.species, dist.species)
  
  if(length(common.species) == 0){
    stop("No species in common between the tree and distribution data.")
  }
  
  message("Using ", length(common.species), " species in common between tree and distribution.")
  distM <- distM[distM$Sp %in% common.species, ]
  
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
  
  if(method == "exclude"){
    message("Excluding ", sum(pd.res$SR == 1), " single-species cells.")
    pd.res <- pd.res[pd.res$SR > 1, ]
  }
  
  if(abundance == TRUE){
    AbCellList <- split(distM$Sp, distM$Cell)
    AbCellList <- lapply(AbCellList, function(x) table(factor(x, levels = unique(distM$Sp))))
  }
  
  PD <- c() 
  set.seed(123)
  
  for(i in 1:length(rownames(pd.res))){
    cell.id <- pd.res$Cell[i]
    sp.list <- unique(distM$Sp[distM$Cell == cell.id])
    
    if(length(sp.list) > 1){
      if(abundance == TRUE){
        # Extract relative abundances for the species in this cell
        ab.list <- AbCellList[[cell.id]][match(sp.list, names(AbCellList[[cell.id]]))]
        rel_ab <- ab.list / sum(ab.list)
        rel_ab <- rel_ab[!duplicated(names(rel_ab))]
        
        # Trim tree to species in the cell
        tree.subset <- keep.tip(tree,
                                     tip = tree$tip.label[tree$tip.label %in% sp.list])
        # Cophenetic distance matrix
        co.dist <- cophenetic.phylo(tree.subset)
        diag(co.dist) <- NA
        
        # Nearest neighbor distance for each species
        minD <- apply(co.dist, 1, min, na.rm = TRUE)
        
        rel_ab <- rel_ab[names(minD)]   # align weights to minD species order
        rel_ab[is.na(rel_ab)] <- 0
        rel_ab <- rel_ab / sum(rel_ab)  # re-normalize (safe if any NAs/zeros)
        PD <- c(PD, sum(minD * rel_ab))
        
      } else {
        # Abundance = FALSE, unweighted MNTD
        tree.subset <- keep.tip(tree,
                                     tip = tree$tip.label[tree$tip.label %in% sp.list])
        co.dist <- cophenetic.phylo(tree.subset)
        diag(co.dist) <- NA
        minD <- apply(co.dist, 1, min, na.rm = TRUE)
        
        PD <- c(PD, mean(minD))
      }
    }else{
      
      pd_single <- pd.taxon(tree, sp.list, method)
      
      if(abundance == TRUE){
        
        rel_ab <- AbCellList[[cell.id]][sp.list] / sum(AbCellList[[cell.id]])
        PD <- c(PD, (pd_single * rel_ab))
        
      }else{
        
        PD <- c(PD, pd_single)
        
      }
    }
  }
  
  pd.res$PD <- PD
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################	
  
  grid0 <- dist$grid
  values(grid0) <- NA
  set.values(grid0, as.numeric(as.character(pd.res$Cell)), as.numeric(pd.res$PD)) 
  
  ##############################################################
  #                         VILMA OBJECT                       #
  ##############################################################
  
  structure(
    list(
      distribution = dist$distribution,
      grid = dist$grid,
      pd.table = pd.res,
      rasters = list(ab.raster = dist$ab.raster,
                     r.raster = dist$r.raster,
                     pd.raster = grid0
      ),
      calculation.method = method,
      abundance = abundance,
      index = "mntd.calc"
    ),
    class = "vilma.pd"
  )
  
  
}
