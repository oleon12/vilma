#' @title Phylogenetic Endemism (PE)
#' @description
#' Compute Phylogenetic Endemism (PE) for each grid cell from a rooted
#' phylogeny and species distributions. Optionally returns Relative PE (RPE)
#' as PE / Faith's PD per cell.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (see \code{\link{points_to_raster}}).
#' @param RPE Logical; if \code{TRUE} also compute Relative PE (RPE = PE / PD). Default \code{TRUE}.
#' @param faith.method Character string specifying how to handle single-species cells.
#'   Options are:
#'   \itemize{
#'     \item \code{"root"} – For cells with single species, PD is calculated considering the root path to the species.
#'     \item \code{"node"} – For cells with single species, PD considers the closest ancestral node length.
#'     \item \code{"exclude"} – Single-species cells are excluded from the calculation (default).
#'   }
#'
#' @details
#' PE measures how much evolutionary history in a cell is geographically
#' restricted, summing branch lengths weighted by (i) the fraction of
#' descendants present in the cell and (ii) the inverse of the number of
#' cells those descendants occupy.
#'
#' @return
#' A \code{vilma.pd} object with:
#' \itemize{
#'   \item \code{distribution} – original species distribution (input).
#'   \item \code{grid} – reference raster grid (from \code{dist}).
#'   \item \code{pd.table} – a data frame with:
#'     \describe{
#'       \item{\code{Cell}}{Cell identifier.}
#'       \item{\code{SR}}{Species richness per cell.}
#'       \item{\code{PE}}{Phylogenetic Endemism value per cell.}
#'       \item{\code{RPE}}{Relative PE (if \code{RPE = TRUE}).}
#'     }
#'   \item \code{rasters} – list with abundance, richness, and PE rasters.
#'   \item \code{index} – the string \code{"pe.calc"}.
#' }
#'
#' @references
#' Rosauer D, Laffan SW, Crisp MD, Donnellan SC, Cook LG (2009) Phylogenetic endemism:
#' a new approach for identifying geographical concentrations of evolutionary history.
#' \emph{Molecular Ecology} 18(19):4061–4072. \doi{10.1111/j.1365-294X.2009.04311.x}
#'
#' Mishler BD, Knerr N, González-Orozco CE, Thornhill AH, Laffan SW, Miller JT (2014)
#' Phylogenetic measures of biodiversity and neo- and paleo-endemism in Australian \emph{Acacia}.
#' \emph{Nature Communications} 5:4473. \doi{10.1038/ncomms5473}
#'
#' Faith DP (1992) Conservation evaluation and phylogenetic diversity.
#' \emph{Biological Conservation} 61(1):1–10. \doi{10.1016/0006-3207(92)91201-3}
#'
#' @seealso \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- example_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' pe <- pe.calc(tree, dist, RPE = TRUE)
#' print(pe)
#' plot(pe)
#' view.vilma(pe)
#' }
#'
#' @export

pe.calc <- function(tree, dist, RPE = c(TRUE,FALSE), faith.method = c("node","root","exclude")){
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
  
  if(length(RPE) > 1){
    RPE <- TRUE
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
            " species from 'dist' not found in 'tree': ",
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
  
  tree.subset <- keep.tip(phy = tree, tip = as.character(unique(distM$Sp)))
  
  PD <- c()
  
  for(i in 1:length(rownames(pd.res))){
    cell.id <- pd.res$Cell[i]
    sp.list <- unique(distM$Sp[distM$Cell == cell.id])
    
    PE <- 0  # reset for each cell
    
    for(j in 1:nrow(tree.subset$edge)){
      
      child <- tree.subset$edge[j,2]
      branch <- tree.subset$edge.length[j]
      
      #Find all descendant tips of this branch
      if(child <= length(tree.subset$tip.label)){
        #Child is tip
        tips.names <- tree.subset$tip.label[child]
      }else{
        #Child is node
        desc_nodes <- getDescendants(tree.subset, child)
        tip_indices <- desc_nodes[desc_nodes <= length(tree.subset$tip.label)]
        tips.names <- tree.subset$tip.label[tip_indices]
      }
      
      # Fraction of descendant species present in the cell
      f <- sum(tips.names %in% sp.list) / length(tips.names)
      
      #Checl for null species
      if(length(f) == 0){
        f <- 0
      }
      
      # Geographic restriction: number of cells containing any descendant of this branch
      cells_with_tips <- unique(distM$Cell[distM$Sp %in% tips.names])
      
      #Check for null species
      if (length(cells_with_tips) == 0) {
        w <- 0
      } else {
        w <- 1 / length(cells_with_tips)
      }
      
      # Contribution to PE
      PE <- PE + branch * f * w
    }
    
    PD <- c(PD, PE)
  }
  
  pd.res$PD <- PD
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  # If multiple methods provided, default to "node"
  if(length(faith.method) > 1) {
    faith.method <- "exclude"
  }
  
  if(RPE == TRUE) {
    if(faith.method == "exclude") {
      # For exclude, calculate PD with "root" (or "node") as baseline
      fPD <- suppressMessages(faith.pd(tree, dist, method = "root"))
    }else{
      fPD <- suppressMessages(faith.pd(tree, dist, method = faith.method))
    }
  
    # Compute RPE
    pd.res$RPE <- PD / fPD$pd.table$PD
  
    # Mask singleton cells if exclude was requested
    if(faith.method == "exclude") {
      pd.res[pd.res$SR == 1, "RPE"] <- NA
    }
  }
  
  
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
      index = "pe.calc"
    ),
    class = "vilma.pd"
  )
  
  
}

