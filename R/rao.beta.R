#' Rao Beta (phylogenetic) dissimilarity index
#' @description
#' Computes **between-community** phylogenetic dissimilarity (Rao β) among spatial
#' cells using a rooted phylogeny and a \code{vilma.dist} object. Rao β is the
#' expected phylogenetic dissimilarity between two individuals drawn from
#' different communities:
#' \deqn{\Delta_{\mathrm{Rao}}(i,j) = -\tfrac{1}{2}\left(Q_i + Q_j - 2\,p_i^\top D\,p_j\right),}
#' where \eqn{Q = p^\top D p} is Rao’s α within a community, \eqn{p} are within-cell
#' relative abundances (presence–absence if \code{abundance = FALSE}), and
#' \eqn{D} is the patristic distance matrix among species.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (e.g., from \code{points.to.raster()}),
#'   containing at least columns \code{Cell} and \code{Sp} in \code{$distribution}, and a
#'   raster \code{$grid} for mapping outputs.
#' @param abundance Logical. If \code{FALSE} (default), the Cell×Sp table is binarized
#'   (presence–absence) before proportions are computed. If \code{TRUE}, raw Cell×Sp
#'   counts from \code{table(Cell, Sp)} define relative abundances.
#' @param scale01 Logical. If \code{TRUE} (default), \eqn{D} is divided by its maximum
#'   off-diagonal so \eqn{D \in [0,1]} (diagonal set to 0). Under this scaling, Rao β
#'   lies in \eqn{[0,1]} up to numerical tolerance.
#'
#' @details
#' Species are intersected between \code{tree$tip.label} and \code{dist$distribution$Sp};
#' only the common set is used. To ensure correct behavior under taxa-label null
#' models, the patristic matrix \eqn{D} is **reindexed to match** the community
#' matrix column order before matrix multiplications. Empty cells (no species after
#' optional binarization) have their corresponding rows/columns of the output
#' dissimilarity set to \code{NA}.
#'
#' In addition to the full pairwise Rao β matrix, the function returns:
#' \itemize{
#'   \item the mean Rao β per cell (raster),
#'   \item a PCoA of the dissimilarity (first two axes rasterized; Cailliez correction),
#'   \item a 2D NDMS (first two axes rasterized) and the stress value.
#' }
#'
#' @return An object of class \code{vilma.beta}:
#' \itemize{
#'   \item \code{distribution} – Original \code{dist$distribution}.
#'   \item \code{Rao.Beta} – Pairwise Rao β dissimilarity as a \code{stats::dist} object.
#'   \item \code{rasters} – List with:
#'         \code{mean.dissimilarity}, \code{pcoa.1}, \code{pcoa.2},
#'         \code{ndms.1}, \code{ndms.2}.
#'   \item \code{pcoa.eig} – Relative corrected eigenvalues for axes 1–2.
#'   \item \code{ndms.stress} – NDMS stress value.
#'   \item \code{calculation.method} – Summary string of \code{abundance} and \code{scale01}.
#'   \item \code{algorithm} – The string \code{"RaoBeta"}.
#' }
#'
#' @section Interpretation:
#' With \code{scale01 = TRUE}, Rao β ranges from \code{0} (identical phylogenetic
#' composition) to \code{1} (maximal turnover). A similarity matrix can be obtained
#' as \code{1 - as.matrix(Rao.Beta)} if needed.
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
#' \code{\link{rao.calc}} (Rao α), \code{\link{rao.calc.null}} (null models),
#' \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' out <- rao.beta(tree = tree, dist = dist, abundance = FALSE, scale01 = TRUE)
#' print(out)
#' plot(out)
#' view.vilma(out)
#' }
#'
#' @export

rao.beta <- function(tree, dist, abundance = FALSE, scale01 = TRUE){
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
    stop("Rao-Beta (phylogenetic) requires a rooted tree.")
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
  P[!is.finite(P)] <- 0        ### NEW: move NaN cleanup before multiplying
  
  M <- P %*% D
  alpha <- rowSums(M * P)
  names(alpha) <- rownames(P)  ### NEW: keep names here
  
  ###############################################################
  #                            BETA                             #
  ###############################################################
  
  C  <- M %*% t(P)
  DR <- -0.5 * (outer(alpha, alpha, "+") - 2 * C)  # <-- sign & 1/2 factor
  DR[DR < 0 & DR > -1e-12] <- 0                    # numerical cleanup
  diag(DR) <- 0
  rownames(DR) <- rownames(P)     ### NEW: set dimnames
  colnames(DR) <- rownames(P)
  
  ### NEW: set rows/cols to NA for empty cells
  if (any(sum_comm == 0)) {
    idx0 <- which(sum_comm == 0)
    DR[idx0, ] <- NA_real_
    DR[, idx0] <- NA_real_
  }
  
  #############################################################
  #                      Mean Rao Beta                        #
  #############################################################
  
  mean.sym <- DR
  diag(mean.sym) <- NA
  
  mean.sym <- apply(mean.sym, 1, function(x){mean(x, na.rm = T)})
  
  #############################################################
  #          Multidimensional reductions PCoA + NDMS          #
  #############################################################
  
  pcoa <- pcoa(as.dist(DR), correction = "cailliez")
  
  eig.dis <- pcoa$values$Rel_corr_eig[1:2]
  
  pcoa.dis.1 <-pcoa$vectors[, 1]
  pcoa.dis.2 <-pcoa$vectors[, 2]
 
  ndms <- monoMDS(as.dist(DR), k = 2)
  ndms.stress.dis <- ndms$stress
  
  ndms.dis.1 <- ndms$points[, 1]
  ndms.dis.2 <- ndms$points[, 2]
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  #### Mean Rao ####
  
  mean.rao <- dist$grid
  
  values(mean.rao) <- NA
  
  set.values(mean.rao, as.numeric(names(mean.sym)),  as.numeric(mean.sym))
  
  #### Raster PCoA ####
  
  pcoa.r1 <- dist$grid
  pcoa.r2 <- dist$grid
  
  values(pcoa.r1) <- NA
  values(pcoa.r2) <- NA
  
  set.values(pcoa.r1, as.numeric(names(pcoa.dis.1)),  as.numeric(pcoa.dis.1))
  set.values(pcoa.r2, as.numeric(names(pcoa.dis.2)),  as.numeric(pcoa.dis.2))
  
  #### Raster NDMS ####
  
  ndms.r1 <- dist$grid
  ndms.r2 <- dist$grid
  
  values(ndms.r1) <- NA
  values(ndms.r2) <- NA
  
  set.values(ndms.r1, as.numeric(names(ndms.dis.1)),  as.numeric(ndms.dis.1))
  set.values(ndms.r2, as.numeric(names(ndms.dis.2)),  as.numeric(ndms.dis.2))
  
  ##############################################################
  #                          Output                            #
  ##############################################################
  
  structure(
    list(
      distribution = dist$distribution,
      Rao.Beta = as.dist(DR),
      rasters = list(mean.dissimilarity = mean.rao,
                     pcoa.1 = pcoa.r1,
                     pcoa.2 = pcoa.r2,
                     ndms.1 = ndms.r1,
                     ndms.2 = ndms.r2),
      pcoa.eig = eig.dis,
      ndms.stress =ndms.stress.dis,
      calculation.method = paste0("Abudance: ", abundance, " Scale: ", scale01),
      algorithm = "RaoBeta"           
  
    ),
  class = "vilma.beta"
  )
  
  
  
}
