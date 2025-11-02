#' UniFrac Phylogenetic Dissimilarity between communities (cells)
#' @description
#' Computes pairwise UniFrac distances among grid cells and returns:
#' (i) a dissimilarity matrix, and (ii) convenient rasters summarizing
#' mean dissimilarity per cell as well as ordination (PCoA and NMDS)
#' axes mapped to space.
#'
#' The computation follows the classic UniFrac definitions:
#' when \code{method = "unweighted"}, distances depend only on presence/absence;
#' when \code{method = "weighted"}, distances are normalized and depend on
#' the relative abundances of taxa in each cell.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (e.g., from \code{points.to.raster()}),
#'   containing species occurrences per grid cell and a template raster.
#' @param method Character; UniFrac variant to compute:
#'   \itemize{
#'     \item \code{"unweighted"} – presence/absence UniFrac (Default).
#'     \item \code{"weighted"} – normalized weighted UniFrac using relative abundances.
#'   }
#'
#' @details
#' Let \eqn{E} be the set of tree edges with lengths \eqn{w_e}.
#' For two communities (cells) \eqn{A} and \eqn{B}, define \eqn{U_e(A)} and \eqn{U_e(B)}
#' as indicators that edge \eqn{e} has at least one descendant tip present in
#' community \eqn{A} (or \eqn{B}), respectively.
#'
#' \strong{Unweighted UniFrac} (presence/absence) is:
#' \deqn{\mathrm{UniFrac}_{\mathrm{unw}}(A,B) =
#'       \frac{\sum_{e \in E} w_e \cdot \mathbf{1}\{ U_e(A) + U_e(B) = 1 \}}
#'            {\sum_{e \in E} w_e \cdot \mathbf{1}\{ U_e(A) + U_e(B) \ge 1 \}}.}
#'
#' \strong{Weighted UniFrac} (normalized; relative abundances) defines the fraction of
#' abundance below each edge as \eqn{D_A(e)} and \eqn{D_B(e)} (each in \eqn{[0,1]}), and uses:
#' \deqn{\mathrm{UniFrac}_{\mathrm{w}}(A,B) =
#'       \frac{\sum_{e \in E} w_e \, |D_A(e) - D_B(e)|}
#'            {\sum_{e \in E} w_e \, [D_A(e) + D_B(e)]}.}
#'
#' Only species shared between \code{tree} and \code{dist} are used; informative
#' messages are printed for mismatches. The returned distance matrix has diagonal 0.
#'
#' After computing the matrix, the function summarizes each cell by its mean
#' dissimilarity to all other cells and performs ordinations on the dissimilarity
#' matrix using PCoA (with Cailliez correction) and NMDS, returning the first
#' two axes rasterized to the study grid.
#'
#' @return An object of class \code{vilma.beta}, a list with:
#' \itemize{
#'   \item \code{distribution} – Original species distribution (from \code{dist}).
#'   \item \code{UniFrac} – UniFrac dissimilarity as a \code{dist} object (labels are cell IDs).
#'   \item \code{rasters} – A named list of \code{SpatRaster} layers:
#'     \itemize{
#'       \item \code{mean.unifrac} – Raster of mean UniFrac per cell (diagonal excluded).
#'       \item \code{pcoa.1}, \code{pcoa.2} – PCoA axes 1 and 2 mapped to cells (from dissimilarity).
#'       \item \code{ndms.1}, \code{ndms.2} – NMDS axes 1 and 2 mapped to cells (from dissimilarity).
#'     }
#'   \item \code{pcoa.eig} – PCoA eigen values of the two first axes.
#'   \item \code{ndms.stress} – NMDS stress value.
#'   \item \code{calculation.method} – The resolved \code{method} used.
#'   \item \code{algorithm} – The string \code{"UniFrac"}.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item The tree must be rooted and have branch lengths.
#'   \item Only species present in both \code{tree} and \code{dist} are used.
#'   \item Ordinations are performed on the dissimilarity matrix. PCoA uses
#'         \code{ape::pcoa(..., correction = "cailliez")} to handle non-Euclidean cases.
#'   \item \code{method = "weighted"} implements the \emph{normalized} weighted UniFrac
#'         (values in \eqn{[0,1]}), based on relative abundances.
#' }
#'
#' @references
#' Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for
#' comparing microbial communities. \emph{Applied and Environmental Microbiology}, 71(12), 8228–8235.
#'
#' Lozupone, C., Hamady, M., Kelley, S.T., & Knight, R. (2007). Quantitative and
#' qualitative beta diversity measures lead to different insights into factors that
#' structure microbial communities. \emph{Applied and Environmental Microbiology}, 73(5), 1576–1585.
#'
#' Chen, J., Bittinger, K., Charlson, E.S., et al. (2012). Associating microbiome
#' composition with environmental covariates using generalized UniFrac distances.
#' \emph{Bioinformatics}, 28(16), 2106–2113.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{phylosor.calc}}, \code{\link{faith.pd}}, \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' # Unweighted UniFrac (presence/absence)
#' uf_unw <- unifrac.calc(tree, dist, method = "unweighted")
#' print(uf_unw)
#' plot(uf_unw)
#' view.vilma(uf_unw)
#'
#' # Weighted (normalized) UniFrac (relative abundances)
#' uf_w <- unifrac.calc(tree, dist, method = "weighted")
#' }
#'
#'
#' @export


unifrac.calc <- function(tree, dist, method = c("unweighted", "weighted")){
  
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
    stop("UniFrac requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(method)>1){
    method <- "unweighted"
    message("Using method 'unweighted' \n")
  }
  
  ##############################################################
  #                     DATA PREPARATION                       #
  ##############################################################
  
  ##### Dist Processing #####
  distM <- dist$distribution
  
  # Find and use only species that match between tree and dist
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
  
  ###############################################################
  
  D <- descendant.matrix(tree)
  
  cells <- unique(distM$Cell) 
  w <- tree$edge.length 
  
  UniFrac.Dist <- matrix(NA, nrow = length(cells), ncol = length(cells)) 
  colnames(UniFrac.Dist) <- cells
  rownames(UniFrac.Dist) <- cells
  
  cat("\n")
  message("Calculating UniFrac ...")
  cat("\n")
  
  
  pb <- txtProgressBar(min = 0, max = length(cells) , style = 3)
  on.exit(close(pb), add = TRUE)
  
  if(method == "unweighted"){
    for(a in seq_along(cells)){ 
      
      commA <- unique(distM$Sp[which(distM$Cell %in% cells[a])]) 
      Da <- D[ , which(colnames(D)%in%commA)] 
      
      if(!inherits(Da, "matrix")){ 
        Da <- as.matrix(Da) 
      } 
      
      Da <- apply(Da, 1, sum) 
      Da[Da > 0] <- 1 
      
      for(b in seq_along(cells)){ 
        
        commB <- unique(distM$Sp[which(distM$Cell %in% cells[b])]) 
        Db <- D[ , which(colnames(D)%in%commB)] 
        
        if(!inherits(Db, "matrix")){ 
          Db <- as.matrix(Db) 
        }
        
        Db <- apply(Db, 1, sum) 
        Db[Db > 0] <- 1 
        
        union.m <- as.data.frame(cbind(W = w, A = Da, B = Db, diffAB = Da+Db)) 
        union.m <- union.m[-which(union.m$diffAB == 0), ] 
        union.m$diffAB[union.m$diffAB == 2] <- 0 
        
        num <- sum((union.m$W * union.m$diffAB))
        den <- sum(union.m$W) 
        
        if(den > 0){
          UniFrac.Dist[a,b] <- num/den
        }else{
          UniFrac.Dist[a,b] <- 0
        }
      }
      setTxtProgressBar(pb, a) 
    }
  }
  
  if (method == "weighted") {
  for (a in seq_along(cells)) {

    commA <- distM$Sp[distM$Cell %in% cells[a]]
    A.counts <- table(commA)
    A.vec <- rep(0, ncol(D)); names(A.vec) <- colnames(D)
    A.vec[names(A.counts)] <- as.numeric(A.counts)

    sA <- sum(A.vec)
    if (sA == 0) {
      UniFrac.Dist[a, ] <- 0   # fill row a with zeros
      setTxtProgressBar(pb, a)
      next
    }
    A.vec <- A.vec / sA

    Da <- as.numeric(D %*% A.vec)

    for (b in seq_along(cells)) {

      commB <- distM$Sp[distM$Cell %in% cells[b]]
      B.counts <- table(commB)
      B.vec <- rep(0, ncol(D)); names(B.vec) <- colnames(D)
      B.vec[names(B.counts)] <- as.numeric(B.counts)

      sB <- sum(B.vec)
      if (sB == 0) { 
        UniFrac.Dist[a,b] <- 0
        next
      }
      B.vec <- B.vec / sB

      Db <- as.numeric(D %*% B.vec)

      num <- sum(w * abs(Da - Db))
      den <- sum(w * (Da + Db))
      UniFrac.Dist[a,b] <- if (den > 0) num/den else 0
    }
    setTxtProgressBar(pb, a)
  }
  }
  
  
  ################################################################
  #                        Mean PhyloSor                         #
  ################################################################
  
  mean.unifrac <- UniFrac.Dist
  diag(mean.unifrac) <- NA
  
  mean.unifrac <- apply(mean.unifrac, 1, function(x){mean(x, na.rm=T)})
  
  ################################################################
  #            Multidimensional reductions PCoA + NDMS           #
  ################################################################
  
  pcoa <- pcoa(as.dist(UniFrac.Dist), correction = "cailliez")
  eig.dis <- pcoa$values$Rel_corr_eig[1:2]
  
  pcoa.dis.1 <-pcoa$vectors[, 1]
  pcoa.dis.2 <-pcoa$vectors[, 2]
  
  ndms <- monoMDS(as.dist(UniFrac.Dist), k = 2)
  ndms.stress.dis <- ndms$stress
  
  ndms.dis.1 <- ndms$points[, 1]
  ndms.dis.2 <- ndms$points[, 2]
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  #### Raster Mean PhyloSor####
  
  ### Sym
  mean.unifrac.r <- dist$grid
  
  values(mean.unifrac.r) <- NA
  
  set.values(mean.unifrac.r, as.numeric(names(mean.unifrac)),  as.numeric(mean.unifrac))
  
  #### Raster PCoA ####
  
  pcoa.r1 <- dist$grid
  pcoa.r2 <- dist$grid
  
  values(pcoa.r1) <- NA
  values(pcoa.r2) <- NA
  
  set.values(pcoa.r1, as.numeric(names(pcoa.dis.1)),  as.numeric(pcoa.dis.1))
  set.values(pcoa.r2, as.numeric(names(pcoa.dis.2)),  as.numeric(pcoa.dis.2))
  
  #### Raster PCoA ####
  
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
      UniFrac = as.dist(UniFrac.Dist),
      rasters = list(
        mean.unifrac = mean.unifrac.r,
        pcoa.1 = pcoa.r1,
        pcoa.2 = pcoa.r2,
        ndms.1 = ndms.r1,
        ndms.2 = ndms.r2
      ),
      pcoa.eig = eig.dis,
      ndms.stress = ndms.stress.dis,
      calculation.method = method,
      algorithm = "UniFrac"
    ),
    class = "vilma.beta"
  )
  
}
