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
  
  # Keep your default behaviour: if user doesn't pass a single value -> "exclude"
  if (length(method) > 1L) {
    method <- "exclude"
    message("Using method 'exclude'\n")
  } else {
    method <- match.arg(method, c("root", "node", "exclude"))
  }
  
  ##############################################################
  #                     DATA PREPARATION                       #
  ##############################################################
  
  distM <- dist$distribution
  
  tree.species <- tree$tip.label
  dist.species <- unique(distM$Sp)
  
  missing.in.tree <- setdiff(dist.species, tree.species)
  missing.in.dist <- setdiff(tree.species, dist.species)
  
  if (length(missing.in.tree) > 0) {
    message(
      "Note: ", length(missing.in.tree), 
      " species from 'dist' not found in 'tree': ",
      paste(head(missing.in.tree, 5), collapse = ", "),
      ifelse(length(missing.in.tree) > 5, "...", "")
    )
  }
  
  if (length(missing.in.dist) > 0) {
    message(
      "Note: ", length(missing.in.dist),
      " species from 'tree' not found in 'dist': ",
      paste(head(missing.in.dist, 5), collapse = ", "),
      ifelse(length(missing.in.dist) > 5, "...", "")
    )
  }
  
  # Only common species
  common.species <- intersect(tree.species, dist.species)
  if (length(common.species) == 0) {
    stop("No species in common between the tree and distribution data.")
  }
  
  message("Using ", length(common.species), " species in common between tree and distribution.")
  
  # 1) Drop unused tips ONCE
  tips.to.drop <- setdiff(tree$tip.label, common.species)
  if (length(tips.to.drop) > 0) {
    tree <- drop.tip(tree, tips.to.drop)
  }
  tree.species <- tree$tip.label
  
  # 2) Filter distribution to common species
  distM <- distM[distM$Sp %in% tree.species, , drop = FALSE]
  
  ##############################################################
  #        COMMUNITY MATRIX (cells x species) + RICHNESS       #
  ##############################################################
  
  comm <- table(distM$Cell, distM$Sp)
  comm[comm > 0] <- 1L
  comm <- as.matrix(comm)
  
  # Reorder columns to match tree tip order
  comm <- comm[, tree.species, drop = FALSE]
  
  Cell <- rownames(comm)
  SR   <- rowSums(comm)
  
  pd.res <- data.frame(
    Cell = Cell,
    SR   = SR,
    stringsAsFactors = FALSE
  )
  
  ##############################################################
  #         PRECOMPUTE TIP–EDGE RELATIONSHIPS (BIG SPEEDUP)    #
  ##############################################################
  
  .build_descendant_matrix <- function(tree) {
    n_tip  <- length(tree$tip.label)
    n_edge <- nrow(tree$edge)
    
    max_node <- max(tree$edge)
    parent_edge <- integer(max_node)
    parent_edge[] <- NA_integer_
    
    # parent_edge[child_node] = edge_index
    for (e in seq_len(n_edge)) {
      child <- tree$edge[e, 2]
      parent_edge[child] <- e
    }
    
    # Root node is the one that never appears as a child
    root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
    
    # Logical matrix: edges x tips
    D <- matrix(FALSE, nrow = n_edge, ncol = n_tip)
    
    for (tip in seq_len(n_tip)) {
      cur <- tip
      # mark all edges on path tip -> root
      while (!is.na(cur) && cur != root_node) {
        e   <- parent_edge[cur]
        D[e, tip] <- TRUE
        cur <- tree$edge[e, 1]
      }
    }
    
    D
  }
  
  D <- .build_descendant_matrix(tree)
  edge_lengths <- tree$edge.length
  
  ##############################################################
  #                      PD CALCULATION                        #
  ##############################################################
  
  n_cells <- nrow(pd.res)
  PD      <- numeric(n_cells)
  names(PD) <- pd.res$Cell
  
  single_cells <- pd.res$SR == 1L
  multi_cells  <- pd.res$SR >  1L
  
  # If method == "exclude", we will throw singletons away at the end
  if (method == "exclude") {
    n_excl <- sum(single_cells)
    if (n_excl > 0) {
      message("Excluding ", n_excl, " single-species cells.")
    }
  }
  
  # Depths only needed if method == "root" for multi-species cells
  if (method == "root") {
    depths <- node.depth.edgelength(tree)
  }
  
  # 1) Multi-species cells: use descendant matrix instead of drop.tip()
  if (any(multi_cells)) {
    idx_multi <- which(multi_cells)
    
    for (i in idx_multi) {
      cell.id <- pd.res$Cell[i]
      
      # species present in this cell
      row_i   <- comm[cell.id, , drop = FALSE]
      sp_idx  <- which(row_i[1, ] == 1L)
      sp.list <- tree.species[sp_idx]
      
      # branches used by any of the species in this cell
      used_edges <- rowSums(D[, sp_idx, drop = FALSE]) > 0
      PDcore     <- sum(edge_lengths[used_edges])
      
      if (method == "root") {
        mrca_node    <- getMRCA(tree, sp.list)
        root_to_mrca <- depths[mrca_node]
        PD[i]        <- PDcore + root_to_mrca   # keep your original behaviour
      } else {
        PD[i] <- PDcore
      }
    }
  }
  
  # 2) Single-species cells: use your pd.taxon() helper
  if (any(single_cells) && method != "exclude") {
    idx_single <- which(single_cells)
    
    for (i in idx_single) {
      cell.id <- pd.res$Cell[i]
      row_i   <- comm[cell.id, , drop = FALSE]
      sp_idx  <- which(row_i[1, ] == 1L)
      sp.list <- tree.species[sp_idx]
      
      # one species per cell by construction
      PD[i] <- pd.taxon(tree, sp.list, method = method)
    }
  }
  
  # Attach PD to pd.res, dropping singletons if method == "exclude"
  if (method == "exclude") {
    keep <- !single_cells
    pd.res <- pd.res[keep, , drop = FALSE]
    pd.res$PD <- PD[pd.res$Cell]
  } else {
    pd.res$PD <- PD
  }
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  grid0 <- dist$grid
  values(grid0) <- NA
  
  set.values(
    grid0,
    as.numeric(as.character(pd.res$Cell)),
    as.numeric(pd.res$PD)
  )
  
  ##############################################################
  #                        VILMA OBJECT                        #
  ##############################################################
  
  structure(
    list(
      distribution = dist$distribution,
      grid         = dist$grid,
      pd.table     = pd.res,
      rasters      = list(
        ab.raster = dist$ab.raster,
        r.raster  = dist$r.raster,
        pd.raster = grid0
      ),
      calculation.method = method,
      index              = "faith.pd"
    ),
    class = "vilma.pd"
  )
}
