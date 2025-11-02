#' Between-community Mean Pairwise Distance (betaMPD) between communities (cells)
#'
#' @description
#' Computes pairwise between-community MPD among grid cells and returns:
#' (i) a betaMPD matrix (larger = more phylogenetically distant), and
#' (ii) rasters summarizing mean betaMPD per cell as well as
#' ordination (PCoA and NMDS) axes mapped to space.
#'
#' Two variants are provided:
#' \itemize{
#'   \item Unweighted (presence/absence) - betaMPD is the mean patristic
#'         distance over all cross-cell species pairs.
#'   \item Weighted (abundances) - betaMPD is the expectation of the
#'         patristic distance under within-cell relative abundance weights.
#' }
#' Optionally, the result can be normalized by the mean within-cell MPDs and/or
#' rescaled to the unit interval for visualization.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (e.g., from \code{points.to.raster()}),
#'   containing species occurrences per grid cell and a template raster.
#' @param mpd.method Character; how within-cell MPD is handled when
#'   \code{normalize = TRUE}. One of:
#'   \itemize{
#'     \item \code{"exclude"} - single-species cells are excluded (Default).
#'     \item \code{"root"} - singletons use root-to-tip distance.
#'     \item \code{"node"} - singletons use tip-to-nearest-node branch.
#'   }
#'   Ignored if \code{normalize = FALSE}.
#' @param abundance Logical; if \code{TRUE}, compute abundance-weighted betaMPD using
#'   relative abundances within each cell. Default \code{FALSE}.
#' @param exclude.conspecific Logical; if \code{TRUE}, conspecific matches across
#'   cells are excluded from betaMPD (and abundance weights are renormalized in the
#'   weighted case). Default \code{FALSE}.
#' @param normalize Logical; if \code{TRUE}, return the relative (dimensionless)
#'   index \eqn{betaMPD_rel = betaMPD / ((MPD_A + MPD_B)/2)},
#'   where within-cell MPDs follow \code{mpd.method} and \code{abundance}. Default \code{FALSE}.
#' @param scale01 Logical; if \code{TRUE}, divide the final betaMPD matrix by its
#'   maximum (over all finite entries) to rescale values to \eqn{[0,1]} for
#'   visualization. This is a display transform; it does not change the underlying
#'   metric. Default \code{FALSE}.
#'
#' @details
#' Let \eqn{D_{sp}} be the patristic (cophenetic) distance matrix among species.
#' For two cells with species sets \eqn{S_A} and \eqn{S_B}:
#'
#' Unweighted betaMPD:
#' \deqn{betaMPD(A,B) = mean\{ d(i,j) : i \in S_A,\ j \in S_B \}.}
#'
#' Weighted betaMPD (relative abundances \eqn{p_A}, \eqn{p_B}):
#' \deqn{betaMPD_w(A,B) = \sum_{i \in S} \sum_{j \in S} p_A(i)\, p_B(j)\, d(i,j)
#'      = \mathbf{p}_A^\top D_{sp}\, \mathbf{p}_B.}
#'
#' If \code{exclude.conspecific = TRUE}, conspecific pairs are removed; in the weighted
#' case, weights are renormalized within each cell. Comparisons involving an empty cell
#' return \code{NA}. The diagonal of the returned betaMPD matrix is set to 0 (purely
#' between-community distance).
#'
#' When \code{normalize = TRUE}, each pairwise value is divided by the mean of the
#' corresponding within-cell MPDs, computed via \code{mpd.calc} with the provided
#' \code{mpd.method} and \code{abundance} so that numerator and denominator are
#' consistent. Note that \code{normalize} produces a dimensionless ratio (not bounded
#' to \eqn{[0,1]}). If \code{scale01 = TRUE}, the final matrix (after any normalization)
#' is divided by its maximum for plotting.
#'
#' @return An object of class \code{vilma.beta}, a list with:
#' \itemize{
#'   \item \code{distribution} - Original species distribution (from \code{dist}).
#'   \item \code{bMPD} - Pairwise betaMPD as a \code{dist} object (labels are cell IDs).
#'   \item \code{rasters} - A named list of \code{SpatRaster} layers:
#'     \itemize{
#'       \item \code{mean.bMPD} - Raster of mean betaMPD per cell (diagonal excluded).
#'       \item \code{pcoa.1}, \code{pcoa.2} - PCoA axes 1 and 2 mapped to cells (from \code{bMPD}).
#'       \item \code{nmds.1}, \code{nmds.2} - NMDS axes 1 and 2 mapped to cells (from \code{bMPD}).
#'     }
#'   \item \code{pcoa.eig} - PCoA eigen values of the first two axes.
#'   \item \code{nmds.stress} - NMDS stress value.
#'   \item \code{calculation.method} - The resolved method used (see Notes).
#'   \item \code{algorithm} - The string \code{"beta.mpd"}.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item The tree must be rooted and have branch lengths.
#'   \item Only species present in both \code{tree} and \code{dist} are used; informative
#'         messages are printed for mismatches.
#'   \item \code{normalize = TRUE} uses \code{mpd.calc} internally with \code{mpd.method}
#'         and \code{abundance} to compute within-cell MPDs for the denominator.
#'   \item \code{scale01 = TRUE} rescales the final matrix by its maximum for visualization;
#'         values then lie in \eqn{[0,1]} but are not comparable across datasets/trees.
#' }
#'
#' @references
#' Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002).
#' Phylogenies and community ecology. Annual Review of Ecology and Systematics, 33, 475-505.
#'
#' Miller ET, Farine DR, Trisos CH (2017).
#' Phylogenetic community structure metrics and null models: a review with new methods.
#' Ecology Letters, 20(7), 807-818.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{mpd.calc}}, \code{\link{beta.mntd}}, \code{\link{phylosor.calc}}, \code{\link{unifrac.calc}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' # Unweighted betaMPD (raw, in tree units)
#' bm <- beta.mpd(tree, dist, abundance = FALSE)
#' print(bm)
#' plot(bm)
#' view.vilma(bm)
#'
#' # Weighted betaMPD (relative abundances)
#' bm_w <- beta.mpd(tree, dist, abundance = TRUE)
#'
#' # Relative (normalized) betaMPD
#' bm_rel <- beta.mpd(tree, dist, abundance = FALSE,
#'                    normalize = TRUE, mpd.method = "exclude")
#'
#' # Rescale to 0-1 for visualization
#' bm_scaled <- beta.mpd(tree, dist, abundance = TRUE,
#'                       normalize = FALSE, scale01 = TRUE)
#' }
#'
#' @export


beta.mpd <- function(tree, dist, mpd.method = c("exclude","root","node"), 
                     abundance = FALSE, exclude.conspecific = FALSE, normalize = FALSE, scale01=FALSE){
  
  
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
    stop("betaMPD requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(mpd.method)>1){
    mpd.method <- "exclude"
    message("Using MPD method 'exclude' \n")
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
  
  message("Using ", length(common.species), " species in common between tree and distribution. \n")
  
  
  message("Parameters used: \n", "  Abundance: ", abundance, "\n",
          "  Exclude conspecific: ", exclude.conspecific, "\n",
          "  Normalize: ", normalize, "\n",
          "  Scale to 0-1: ", scale01, "\n")
  
  distM <- distM[distM$Sp %in% common.species, ]
  
  ###############################################################
  #                    Data preparation                         # 
  ###############################################################
  
  # Conphenetic distances matrix
  Dsp <- cophenetic.phylo(tree)
  # align to common species
  Dsp <- Dsp[common.species, common.species, drop = FALSE]
  
  #Cells for the pair-wise matrix
  cells <- unique(distM$Cell)
  
  #Out put matrix
  beta.out <- matrix(NA, nrow = length(cells), ncol = length(cells),
                     dimnames = list(cells, cells))
  
  
  
  ############## Abundance Loop #################
  
  if(abundance == TRUE){
    relAb <- list() # Realtive abundance per species per cell
    for(i in seq_along(cells)){
      SpCount <- table(distM$Sp[which(distM$Cell %in% cells[i])])
      SpAb <- SpCount / sum(SpCount)
      relAb[[i]] <- SpAb
    }
    names(relAb) <- cells
  }
  
  ############## Normalization with MPD #####################
  
  if(normalize == TRUE){
    mpd_res <- suppressMessages(mpd.calc(tree = tree, dist = dist, method = mpd.method, abundance = abundance))
    MPD.within <- setNames(mpd_res$pd.table$PD, mpd_res$pd.table$Cell)
  }
  
  #############################################################
  
  message("Calculating betaMPD ...", "\n")
  
  pb <- txtProgressBar(min = 0, max = length(cells) , style = 3)
  
  for (a in seq_along(cells)) {
    
    ComA <- distM[which(distM$Cell %in% cells[a]), ]
    
    for(b in seq_along(cells)){
      # Species sets for each cell
      SpA <- unique(ComA$Sp)
      ComB <- distM[which(distM$Cell %in% cells[b]), ]
      SpB <- unique(ComB$Sp)
      
      # Empty-cell policy
      if(length(SpA) == 0 || length(SpB) == 0){
        beta.out[a, b] <- NA
        next
      }
      
      # Calculate bMPD (unweihgted and weighted)
      
      if(abundance == FALSE){
        
        uA <- SpA
        uB <- SpB
        
        # Excluding conspecific, keep only unique species between community
        if(exclude.conspecific == TRUE){
          shared <- intersect(uA, uB)
          if(length(shared) > 0){
            uA <- setdiff(uA, shared)
            uB <- setdiff(uB, shared)
          }
        }
        
        # If both cells share same species two options
        # When conspecific is TRUE set NA, else set 0
        
        if(length(uA) == 0 || length(uB) == 0){
          beta.out[a, b] <- if(exclude.conspecific) NA else 0
        }else{
          
          # Otherwise, calculate the mean distance between species uA and uB
          
          beta.out[a, b] <- mean(Dsp[uA, uB, drop = FALSE])
        }
      }else{
       
         ## Build probability vectors aligned to Dsp species order (common.species)
        pA <- rep(0, length(common.species)) 
        names(pA) <- common.species
        pB <- rep(0, length(common.species)) 
        names(pB) <- common.species
        
        ## Relative abundances already computed in relAb list
        pA[names(relAb[[a]])] <- as.numeric(relAb[[a]])
        pB[names(relAb[[b]])] <- as.numeric(relAb[[b]])
        
        ## Conspecific exclusion: zero shared species and renormalize each side
        if(exclude.conspecific == TRUE){
          both <- names(which(pA > 0 & pB > 0))
          if(length(both) > 0){
            pA[both] <- 0
            pB[both] <- 0
          }
          sA <- sum(pA) # Sum of proabilities
          sB <- sum(pB)
          if(sA > 0) pA <- pA / sA # How much probability remains after remove common species
          if(sB > 0) pB <- pB / sB # How much probability remains after remove common species
          if(sA == 0 || sB == 0){
            beta.out[a, b] <- NA
            next
          }
        }
        
        ## Expected cross-community distance: pA' * Dsp * pB
        val <- as.numeric(t(pA) %*% Dsp %*% pB)
        
        ## Normalization to betaMPD_rel (if requested)
        if(normalize == TRUE){
          denom <- mean(c(MPD.within[cells[a]], MPD.within[cells[b]]), na.rm = TRUE)
          if(!is.na(denom) && denom > 0){
            val <- val / denom
          }else{
            val <- NA
          }
        }
        
        beta.out[a, b] <- val
      }
      
      
    }
    setTxtProgressBar(pb, a)
  }
  
  diag(beta.out) <- 0
  
  # Optional 0-1 scaling of the beta matrix
  if (scale01 == TRUE) {
    maxD <- max(beta.out, na.rm = TRUE)
    if (maxD > 0) beta.out <- beta.out / maxD
  }
  
  bMPD <- as.dist(beta.out)
  
  ################################################################
  #                        Mean PhyloSor                         #
  ################################################################
  
  mean.bMPD <- beta.out
  diag(mean.bMPD) <- NA
  
  mean.bMPD <- apply(mean.bMPD, 1, function(x){mean(x, na.rm=T)})
  
  ################################################################
  #            Multidimensional reductions PCoA + NDMS           #
  ################################################################
  
  pcoa <- pcoa(bMPD, correction = "cailliez")
  eig.dis <- pcoa$values$Rel_corr_eig[1:2]
  
  pcoa.dis.1 <-pcoa$vectors[, 1]
  pcoa.dis.2 <-pcoa$vectors[, 2]
  
  ndms <- monoMDS(bMPD, k = 2)
  ndms.stress.dis <- ndms$stress
  
  ndms.dis.1 <- ndms$points[, 1]
  ndms.dis.2 <- ndms$points[, 2]
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  #### Raster Mean PhyloSor####
  
  ### Sym
  mean.bMPD.r <- dist$grid
  
  values(mean.bMPD.r) <- NA
  
  set.values(mean.bMPD.r, as.numeric(colnames(beta.out)),  as.numeric(mean.bMPD))
  
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
      bMPD = bMPD,
      rasters = list(
        mean.bMPD = mean.bMPD.r,
        pcoa.1 = pcoa.r1,
        pcoa.2 = pcoa.r2,
        ndms.1 = ndms.r1,
        ndms.2 = ndms.r2
      ),
      pcoa.eig = eig.dis,
      ndms.stress = ndms.stress.dis,
      calculation.method = mpd.method,
      algorithm = "beta.MPD"
    ),
    class = "vilma.beta"
  )
  

}


