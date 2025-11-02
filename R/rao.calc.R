#' Rao's Q (within-community) phylogenetic diversity
#' @description
#' Computes Rao’s quadratic entropy (RaoQ) per spatial cell from a rooted
#' phylogenetic tree and a \code{vilma.dist} distribution. RaoQ is the expected
#' phylogenetic dissimilarity between two individuals randomly drawn from the
#' same community (cell), using patristic distances from the tree.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (e.g., from \code{points.to.raster()}),
#'   whose \code{$distribution} contains at least \code{Cell} and \code{Sp}, and whose
#'   \code{$grid} is used to rasterize outputs.
#' @param abundance Logical. If \code{FALSE} (default), a presence–absence Cell×Sp
#'   table is used (each species weighted equally within a cell). If \code{TRUE},
#'   the Cell×Sp counts from \code{table(Cell, Sp)} are used to compute within-cell
#'   relative abundances.
#' @param scale01 Logical. If \code{TRUE} (default), the patristic distance matrix
#'   is divided by its maximum off-diagonal so that \eqn{D \in [0,1]} (diagonal set
#'   to 0). Under this scaling, RaoQ lies in \eqn{[0,1]} up to numerical tolerance.
#'
#' @details
#' Species are intersected between \code{tree$tip.label} and \code{dist$distribution$Sp};
#' only the common species are used. The patristic matrix \eqn{D} is computed with
#' \code{ape::cophenetic.phylo()} and then **reindexed to match the column order of
#' the community matrix** to ensure correct behavior under taxa-label null models.
#'
#' For each cell, let \eqn{p} be the vector of within-cell relative abundances
#' (or presence–absence proportions). RaoQ is computed as \eqn{Q = p^\top D p}.
#' Empty cells (no species) are returned as \code{NA}.
#'
#' @return An object of class \code{vilma.pd}:
#' \itemize{
#'   \item \code{distribution} – Original \code{dist$distribution}.
#'   \item \code{grid} – Original \code{dist$grid}.
#'   \item \code{pd.table} – Data frame with columns:
#'         \code{Cell} (cell id), \code{SR} (species richness), and \code{PD} (RaoQ).
#'   \item \code{rasters} – List with \code{ab.raster}, \code{r.raster}, and
#'         \code{pd.raster} (RaoQ per cell).
#'   \item \code{calculation.method} – Character summary of the settings.
#'   \item \code{index} – The string \code{"rao"}.
#' }
#'
#' @section Interpretation:
#' Larger RaoQ indicates greater within-community phylogenetic dispersion.
#' With \code{scale01 = TRUE}, values are comparable across trees (range approx. \code{[0,1]}).
#'
#' @references
#' Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified approach.
#' \emph{Theoretical Population Biology}, 21, 24–43. \cr
#' Botta-Dukát, Z. (2005). Rao’s quadratic entropy as a measure of functional diversity.
#' \emph{Community Ecology}, 6, 283–290.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{rao.beta}} for between-community dissimilarity; 
#' \code{\link{rao.calc.null}} for null-model SES; 
#' \code{\link{points_to_raster}} to build \code{vilma.dist}.
#'
#' @examples
#' \dontrun{
#' 
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' rao <- rao.calc(tree = tree, dist = dist, abundance = FALSE, scale01 = TRUE)
#' print(rao)
#' plot(rao)
#' view.vilma(rao)
#' }
#'
#' @export


rao.calc <- function(tree, dist, abundance = FALSE, scale01 = TRUE){
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
    stop("Rao (phylogenetic) requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
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
  
  
  ###############################################################
  
  D <- as.matrix(cophenetic.phylo(tree))
  D <- D[common.species, common.species, drop =FALSE]
  
  if(scale01 == TRUE){
    diag(D) <- NA
    D <- D/max(D, na.rm = T)
    diag(D) <- 0
  }else{
    diag(D) <- 0
  }
  
  comm <- as.matrix(table(distM$Cell, distM$Sp))
  
  if(abundance == F){
    comm <- (comm > 0) * 1L
  }
  
  
  ###############################################################
  
  ### Avoid issues with taxa.label for null models
  
  sp_comm <- colnames(comm)
  
  D <- D[sp_comm, sp_comm, drop = FALSE]
  
  sum_comm <- apply(comm, 1, sum)
  sum_comm
  
  if (any(sum_comm == 0)) {
    message("Note: ", sum(sum_comm == 0), " empty cell(s) with total = 0. Returning NA for those.")
  }
  
  P <- comm/sum_comm
  
  M <- P %*% D
  alpha <- rowSums(M * P)
  
  P <- comm / sum_comm
  P[!is.finite(P)] <- 0
  
  M <- P %*% D
  alpha <- rowSums(M * P)
  names(alpha) <- rownames(P)
  
  alpha[sum_comm == 0] <- NA_real_
  
  ##############################################################
  
  pd.res$PD <- alpha
  
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
      calculation.method = paste0("Abundance = ", abundance),
      index = "rao.calc"
    ),
    class = "vilma.pd"
  )
  
  
}

