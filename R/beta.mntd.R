#' Between-community Mean Nearest Taxon Distance (betaMNTD) between communities (cells)
#'
#' @description
#' Computes pairwise **betaMNTD** among grid cells and returns:
#' (i) a betaMNTD distance matrix (larger = more phylogenetically distant at the tip scale), and
#' (ii) rasters summarizing mean betaMNTD per cell as well as ordination
#' (PCoA and NMDS) axes mapped to space.
#'
#' Two variants are provided:
#' \itemize{
#'   \item **Unweighted (presence/absence)** - betaMNTD is the mean of nearest-taxon
#'         distances in both directions (A to B and B to A).
#'   \item **Weighted (abundances)** - betaMNTD is the expected nearest-taxon distance
#'         using within-cell relative abundances as sampling probabilities.
#' }
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} (e.g., from \code{points.to.raster()}),
#'   containing species occurrences per grid cell and a raster grid.
#' @param mntd.method Character; how within-cell MNTD is handled when
#'   \code{normalize = TRUE}. Options:
#'   \itemize{
#'     \item \code{"exclude"} - single-species cells excluded (Default).
#'     \item \code{"root"} - singletons use root-to-tip distance.
#'     \item \code{"node"} - singletons use tip-to-nearest-node branch.
#'   }
#' @param abundance Logical; if TRUE, compute abundance-weighted betaMNTD.
#' @param exclude.conspecific Logical; if TRUE, conspecific matches across
#'   cells are excluded. Default FALSE.
#' @param normalize Logical; if TRUE, return
#'   \eqn{betaMNTD_rel = betaMNTD / ((MNTD_A + MNTD_B)/2)}.
#' @param scale01 Logical; if TRUE, divide betaMNTD matrix by its maximum to
#'   rescale values to \eqn{[0,1]}. Default FALSE.
#'
#' @details
#' Let \eqn{D_sp} be the patristic distance matrix among species. For two cells:
#'
#' **Unweighted betaMNTD** (symmetrized):
#' \deqn{
#' betaMNTD(A,B) =
#' 0.5 * [
#'   mean_{i in S_A} min_{j in S_B} d(i,j)
#'   +
#'   mean_{j in S_B} min_{i in S_A} d(j,i)
#' ]
#' }
#'
#' **Weighted betaMNTD**:
#' \deqn{
#' betaMNTD_w(A,B) =
#' 0.5 * [
#'   sum p_A(i) * min_{j in S_B} d(i,j)
#'   +
#'   sum p_B(j) * min_{i in S_A} d(j,i)
#' ]
#' }
#'
#' Ordinations (PCoA, NMDS) are computed from the resulting distance matrix.
#'
#' @return A \code{vilma.beta} list containing:
#' \itemize{
#'   \item \code{distribution} - species distribution
#'   \item \code{bMNTD} - distance object (pairwise betaMNTD)
#'   \item \code{rasters} - SpatRaster layers:
#'     \itemize{
#'       \item \code{mean.bMNTD}
#'       \item \code{pcoa.1}, \code{pcoa.2}
#'       \item \code{nmds.1}, \code{nmds.2}
#'     }
#'   \item \code{pcoa.eig} - eigenvalues
#'   \item \code{nmds.stress} - NMDS stress
#'   \item \code{calculation.method}
#'   \item \code{algorithm} - "beta.mntd"
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item Tree must be rooted and have branch lengths.
#'   \item Only shared species are used.
#'   \item \code{normalize = TRUE} uses \code{mntd.calc}.
#'   \item \code{scale01 = TRUE} rescales values for visualization.
#' }
#'
#' @references
#' Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002)
#' Phylogenies and community ecology. Annual Review of Ecology and Systematics 33:475-505.
#'
#' Miller ET, Farine DR, Trisos CH (2017)
#' Phylogenetic community structure metrics and null models. Ecology Letters 20:807-818.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{mntd.calc}}, \code{\link{beta.mpd}}, \code{\link{phylosor.calc}}, \code{\link{unifrac.calc}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#' # Unweighted betaMNTD
#' bn <- beta.mntd(tree, dist, abundance = FALSE)
#' print(bn)
#' plot(bn)
#'
#' # Weighted betaMNTD
#' bn_w <- beta.mntd(tree, dist, abundance = TRUE)
#'
#' # Relative betaMNTD
#' bn_rel <- beta.mntd(tree, dist, abundance = FALSE, normalize = TRUE)
#'
#' # Scaled to 0-1
#' bn_scaled <- beta.mntd(tree, dist, abundance = TRUE, scale01 = TRUE)
#' }
#'
#' @export


beta.mntd <- function(tree, dist, mntd.method = c("exclude","root","node"),
                      abundance = FALSE, exclude.conspecific = FALSE, normalize = FALSE, scale01 = FALSE){
                      
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
    stop("betaMNTD requires a rooted tree.")
  }
  
  if (!inherits(dist, "vilma.dist")) {
    stop("Input 'dist' must be an object of class 'vilma.dist'. See 'points_to_raster()' function.")
  }
  
  if(length(mntd.method)>1){
    mntd.method <- "exclude"
    message("Using MNTD method 'exclude' \n")
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
    mntd_res <- suppressMessages(mntd.calc(tree = tree, dist = dist, method = mntd.method, abundance = abundance))
    MNTD.within <- setNames(mntd_res$pd.table$PD, mntd_res$pd.table$Cell)
  }                    
                      
  #############################################################
  
  message("Calculating betaMNTD ...", "\n")
  
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
          
          # Unweighted betaMNTD:
          dmin.AB <- apply(Dsp[uA, uB, drop = FALSE], 1, min)
          dmin.BA <- apply(Dsp[uB, uA, drop = FALSE], 1, min)
          beta.out[a,b] <- mean(c(mean(dmin.AB), mean(dmin.BA)))
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
        
        uA <- names(pA)[pA > 0]
        uB <- names(pB)[pB > 0]

        if(length(uA) == 0 || length(uB) == 0){
          beta.out[a, b] <- NA
          next
        }

        # Nearest-neighbor distances in each direction
        dmin.AB <- apply(Dsp[uA, uB, drop = FALSE], 1, min)  # for each i in A, nearest j in B
        dmin.BA <- apply(Dsp[uB, uA, drop = FALSE], 1, min)  # for each j in B, nearest i in A

        # Abundance-weighted means of those minima
        val.AB <- sum(pA[uA] * dmin.AB)   # weights = pA over A
        val.BA <- sum(pB[uB] * dmin.BA)   # weights = pB over B

        # Symmetrize
        val <- 0.5 * (val.AB + val.BA)

        # Normalization (unchanged)
        if(normalize == TRUE){
          denom <- mean(c(MNTD.within[cells[a]], MNTD.within[cells[b]]), na.rm = TRUE)
          if(!is.na(denom) && denom > 0) val <- val / denom else val <- NA
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
  
  bMNTD <- as.dist(beta.out)
  
  ################################################################
  #                        Mean PhyloSor                         #
  ################################################################
  
  mean.bMNTD <- beta.out
  diag(mean.bMNTD) <- NA
  
  mean.bMNTD <- apply(mean.bMNTD, 1, function(x){mean(x, na.rm=T)})
  
  ################################################################
  #            Multidimensional reductions PCoA + NDMS           #
  ################################################################
  
  pcoa <- pcoa(bMNTD, correction = "cailliez")
  eig.dis <- pcoa$values$Rel_corr_eig[1:2]
  
  pcoa.dis.1 <-pcoa$vectors[, 1]
  pcoa.dis.2 <-pcoa$vectors[, 2]
  
  ndms <- monoMDS(bMNTD, k = 2)
  ndms.stress.dis <- ndms$stress
  
  ndms.dis.1 <- ndms$points[, 1]
  ndms.dis.2 <- ndms$points[, 2]
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  #### Raster Mean PhyloSor####
  
  ### Sym
  mean.bMNTD.r <- dist$grid
  
  values(mean.bMNTD.r) <- NA
  
  set.values(mean.bMNTD.r, as.numeric(colnames(beta.out)),  as.numeric(mean.bMNTD))
  
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
      bMNTD = bMNTD,
      rasters = list(
        mean.bMNTD = mean.bMNTD.r,
        pcoa.1 = pcoa.r1,
        pcoa.2 = pcoa.r2,
        ndms.1 = ndms.r1,
        ndms.2 = ndms.r2
      ),
      pcoa.eig = eig.dis,
      ndms.stress = ndms.stress.dis,
      calculation.method = mntd.method,
      algorithm = "beta.MNTD"
    ),
    class = "vilma.beta"
  )                    
                      
                      
}
