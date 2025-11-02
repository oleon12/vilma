#' Phylogenetic Beta Diversity Partitioning (β_sor, β_sim, β_sne)
#' @description
#' Computes pairwise **phylogenetic** beta diversity between communities (cells)
#' by partitioning total dissimilarity into **turnover** and **nestedness**
#' components, using branch-length–based analogs of Baselga's decomposition.
#'
#' Two modes are available:
#' \itemize{
#'   \item \code{method = "unweighted"} – presence/absence on the tree.
#'   \item \code{method = "weighted"} – relative-abundance weighting along edges.
#' }
#' Optional \code{normalize = TRUE} (default) rescales per-community abundances to
#' sum to 1 before edge aggregation, so indices lie in \eqn{[0,1]} and are comparable
#' across communities with different total abundances.
#'
#' @param tree A rooted phylogenetic tree of class \code{phylo} with branch lengths.
#' @param dist An object of class \code{vilma.dist} containing species-by-cell
#'   occurrences (columns typically: \code{Sp}, \code{Lon}, \code{Lat}, \code{Cell})
#'   and a reference grid.
#' @param method Character, either \code{"unweighted"} (default) or \code{"weighted"}.
#'   See \emph{Mathematics} for definitions.
#' @param normalize Logical; if \code{TRUE} (default) and \code{method = "weighted"},
#'   converts within-cell abundances to relative abundances (sum to 1) prior to
#'   edge aggregation. Ignored for \code{method = "unweighted"}.
#'
#' @section Mathematics:
#' Let \eqn{E} be the set of tree edges with lengths \eqn{L_e}.
#' For two communities (cells) \eqn{A} and \eqn{B}, define for each edge \eqn{e}
#' the fraction of the community that lies \emph{below} that edge:
#' \itemize{
#'   \item \emph{Unweighted (presence/absence):}
#'     \eqn{p_A(e) = \mathbf{1}\{\text{any descendant tip of } e \text{ in } A\}},
#'     and similarly \eqn{p_B(e)}.
#'   \item \emph{Weighted (relative abundance):}
#'     \eqn{p_A(e) = \sum_{i \in \text{tips}(e)} \tilde{a}_i}, where
#'     \eqn{\tilde{a}_i = a_i / \sum_j a_j} are per-tip relative abundances in
#'     community \eqn{A} (and analogously \eqn{p_B(e)} for \eqn{B}).
#' }
#' From these, define branch-length partitions:
#' \deqn{a = \sum_{e \in E} L_e \, \min\{p_A(e), p_B(e)\},}
#' \deqn{b = \sum_{e \in E} L_e \, \max\{p_A(e) - p_B(e),\, 0\}, \qquad
#'       c = \sum_{e \in E} L_e \, \max\{p_B(e) - p_A(e),\, 0\}.}
#' The phylogenetic beta components are then
#' \deqn{\beta_{\mathrm{sor}} = \frac{b + c}{\,2a + b + c\,}, \qquad
#'       \beta_{\mathrm{sim}} = \frac{\min(b,c)}{\,a + \min(b,c)\,}, \qquad
#'       \beta_{\mathrm{sne}} = \beta_{\mathrm{sor}} - \beta_{\mathrm{sim}}.}
#' These reduce to Baselga's (2010) PA formulas when \eqn{p_A(e), p_B(e) \in \{0,1\}}.
#' By construction \eqn{0 \le \beta_{\mathrm{sim}}, \beta_{\mathrm{sne}}, \beta_{\mathrm{sor}} \le 1}
#' and \eqn{\beta_{\mathrm{sor}} = \beta_{\mathrm{sim}} + \beta_{\mathrm{sne}}}.
#'
#' @details
#' Only species present in both \code{tree} and \code{dist} are used; informative
#' messages report mismatches. Pairwise dissimilarities are computed for all
#' cell combinations. Ordinations (PCoA with Cailliez correction; NMDS) are run
#' on each dissimilarity to obtain spatial axes, and per-cell means are computed
#' by averaging off-diagonal dissimilarities.
#'
#' @return An object of class \code{vilma.beta}:
#' \itemize{
#'   \item \code{total.dissimilarity} – \code{dist} object of \eqn{\beta_{\mathrm{sor}}}.
#'   \item \code{turnover} – \code{dist} object of \eqn{\beta_{\mathrm{sim}}}.
#'   \item \code{nestedness} – \code{dist} object of \eqn{\beta_{\mathrm{sne}}}.
#'   \item \code{rasters} – named list with:
#'     \itemize{
#'       \item \code{mean.total}, \code{mean.turnover}, \code{mean.nestedness}.
#'       \item PCoA rasters: \code{total.pcoa.1}, \code{total.pcoa.2},
#'         \code{turnover.pcoa.1}, \code{turnover.pcoa.2},
#'         \code{nestedness.pcoa.1}, \code{nestedness.pcoa.2}.
#'       \item NMDS rasters: \code{total.ndms.1}, \code{total.ndms.2},
#'         \code{turnover.ndms.1}, \code{turnover.ndms.2},
#'         \code{nestedness.ndms.1}, \code{nestedness.ndms.2}.
#'     }
#'   \item \code{stats} – data frame summarizing eigenvalues (PCoA) and stress (NMDS).
#'   \item \code{calculation.method} – resolved \code{method} and whether normalization was applied.
#'   \item \code{algorithm} – the string \code{"PhyloBeta"}.
#' }
#'
#' @section Implementation Notes:
#' \itemize{
#'   \item Presence/absence mode is obtained by setting \eqn{p_A(e), p_B(e)} to
#'         edge-level indicators of occupancy.
#'   \item In weighted mode, setting \code{normalize = TRUE} yields relative-abundance
#'         fractions (sum to 1 per cell) before edge aggregation; \code{FALSE} uses
#'         raw counts and can emphasize differences in total sampling intensity.
#'   \item The distance matrices are symmetric with zero diagonals; internal checks
#'         confirm bounds and \eqn{\beta_{\mathrm{sor}} = \beta_{\mathrm{sim}} + \beta_{\mathrm{sne}}}.
#'   \item A progress bar can be used during pairwise loops; when present, it should
#'         be closed with \code{on.exit(close(pb), add = TRUE)}.
#' }
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity.
#' \emph{Global Ecology and Biogeography}, 19, 134–143.  
#' Leprieur, F., Albouy, C., de Bortoli, J., Cowman, P.F., Bellwood, D.R., & Mouillot, D. (2012).
#' Quantifying phylogenetic beta diversity: distinguishing turnover of lineages from PD gradients.
#' \emph{PLoS ONE}, 7(8), e42760.  
#' Nipperess, D.A., & Matsen IV, F.A. (2013). The mean and variance of phylogenetic diversity under rarefaction.
#' \emph{Methods in Ecology and Evolution}, 4, 566–572.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso \code{\link{unifrac.calc}}, \code{\link{phylosor.calc}}, \code{\link{faith.pd}},
#'   \code{\link{points_to_raster}}
#'
#' @examples
#' \dontrun{
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' # Presence/absence (Baselga-style) phylogenetic β
#' pb_unw <- phylo.beta(tree, dist, method = "unweighted")
#' print(pb_unw)
#' plot(pb_unw)
#' view.vilma(pb_unw)
#'
#' # Abundance-weighted phylogenetic β (relative abundances by default)
#' pb_w <- phylo.beta(tree, dist, method = "weighted", normalize = TRUE)
#' }
#'
#' @export

phylo.beta <- function(tree, dist, method = c("unweighted", "weighted"), normalize = TRUE){
  
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
    stop("Phylo Beta Diversity requires a rooted tree.")
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
  
  message("Using ", length(common.species), " species in common between tree and distribution. \n")
  
  distM <- distM[distM$Sp %in% common.species, ]
  
  ###############################################################
  #                    Data preparation                         # 
  ###############################################################
  
  cells <- unique(distM$Cell)
  
  #Out put matrices
  beta.sor <- matrix(NA, nrow = length(cells), ncol = length(cells), dimnames = list(cells, cells))
  beta.sim <- matrix(NA, nrow = length(cells), ncol = length(cells), dimnames = list(cells, cells))
  beta.sne <- matrix(NA, nrow = length(cells), ncol = length(cells), dimnames = list(cells, cells))
  
  # Just in case denominator is 0
  safe_div <- function(x, y){
    if(y == 0){
      0
    }else{
      x / y
    }
  }
  
  desc_tips <- .desc_tips_by_edge(tree)
  
  cat("\n")
  message("Calculating Phylo Beta Diversity ...")
  cat("\n")
  
  pb <- txtProgressBar(min = 0, max = length(cells) , style = 3)
  on.exit(close(pb), add = TRUE)
  
  
  for(a in seq_along(cells)){
    
    # 1st community
    commA <- unique(distM$Sp[which(distM$Cell %in% cells[a])])
    
    for(b in seq_along(cells)){
      
      if(a == b){
        beta.sor[a,b] <- 0
        beta.sim[a,b] <- 0
        beta.sne[a,b] <- 0
        next
      }
      
      
      # 2nd community
      commB <- unique(distM$Sp[which(distM$Cell %in% cells[b])])
      
      # Shared and unique species
      if(method == "weighted"){
        
        # Convert character species names to numeric abundances
        commA <- setNames(rep(1, length(commA)), commA)
        commB <- setNames(rep(1, length(commB)), commB)
        
        bl.sums <- branch.partition.weighted(tree, desc_tips, commA, commB, normalize = normalize)
      }else{
        bl.sums <- branch.partition(tree = tree, commA = commA, commB = commB)
      }
    
      blShared <- as.numeric(bl.sums["blShare"])
      blA <- as.numeric(bl.sums["blA"])
      blB <- as.numeric(bl.sums["blB"])
      
      beta.sor[a,b] <- safe_div(x =(blA + blB) , y = (2*blShared + blA + blB))
      beta.sim[a,b] <- safe_div(x = pmin(blA, blB), y = (blShared + pmin(blA, blB)))
      beta.sne[a,b] <- safe_div(x =(blA + blB) , y = (2*blShared + blA + blB)) - safe_div(x = pmin(blA, blB), y = (blShared + pmin(blA, blB)))
      
      
    }
    setTxtProgressBar(pb, a)
  }
  
  stopifnot(all(beta.sor >= 0 & beta.sor <= 1, na.rm=TRUE))
  stopifnot(all(beta.sim >= 0 & beta.sim <= 1, na.rm=TRUE))
  stopifnot(all(beta.sne >= 0 & beta.sne <= 1, na.rm=TRUE))
  stopifnot(max(abs(beta.sor - (beta.sim + beta.sne)), na.rm=TRUE) < 1e-10)
  
  
  ################################################################
  #                        Mean PhyloSor                         #
  ################################################################
  
  mean.sim <- beta.sim
  diag(mean.sim) <- NA
  
  mean.sim <- apply(mean.sim, 1, function(x){mean(x, na.rm=T)})
  
  mean.sor <- beta.sor
  diag(mean.sor) <- NA
  
  mean.sor <- apply(mean.sor, 1, function(x){mean(x, na.rm=T)})
  
  mean.sne <- beta.sne
  diag(mean.sne) <- NA
  
  mean.sne <- apply(mean.sne, 1, function(x){mean(x, na.rm=T)})
  
  
  ################################################################
  #            Multidimensional reductions PCoA + NDMS           #
  ################################################################
  
  ### Sim
  
  pcoa.sim <- pcoa(as.dist(beta.sim), correction = "cailliez")
  eig.sim <- pcoa.sim$values$Rel_corr_eig[1:2]
  
  pcoa.sim.1 <-pcoa.sim$vectors[, 1]
  pcoa.sim.2 <-pcoa.sim$vectors[, 2]
  
  ndms.sim <- monoMDS(as.dist(beta.sim), k = 2)
  ndms.stress.sim <- ndms.sim$stress
  
  ndms.sim.1 <- ndms.sim$points[, 1]
  ndms.sim.2 <- ndms.sim$points[, 2]
  
  ### Sor
  
  pcoa.sor <- pcoa(as.dist(beta.sor), correction = "cailliez")
  eig.sor <- pcoa.sor$values$Rel_corr_eig[1:2]
  
  pcoa.sor.1 <-pcoa.sor$vectors[, 1]
  pcoa.sor.2 <-pcoa.sor$vectors[, 2]
  
  ndms.sor <- monoMDS(as.dist(beta.sor), k = 2)
  ndms.stress.sor <- ndms.sor$stress
  
  ndms.sor.1 <- ndms.sor$points[, 1]
  ndms.sor.2 <- ndms.sor$points[, 2]
  
  ### Sne
  
  pcoa.sne <- pcoa(as.dist(beta.sne), correction = "cailliez")
  eig.sne <- pcoa.sne$values$Rel_corr_eig[1:2]
  
  pcoa.sne.1 <-pcoa.sne$vectors[, 1]
  pcoa.sne.2 <-pcoa.sne$vectors[, 2]
  
  ndms.sne <- monoMDS(as.dist(beta.sne), k = 2)
  ndms.stress.sne <- ndms.sne$stress
  
  ndms.sne.1 <- ndms.sne$points[, 1]
  ndms.sne.2 <- ndms.sne$points[, 2]
  
  
  ##############################################################
  #                       RASTER CREATION                      #
  ##############################################################
  
  #### Raster Mean Beta Phylo####
  
  ### Sym
  mean.sim.r <- dist$grid
  
  values(mean.sim.r) <- NA
  
  set.values(mean.sim.r, as.numeric(colnames(beta.sim)),  as.numeric(mean.sim))
  
  ### Sor
  mean.sor.r <- dist$grid
  
  values(mean.sor.r) <- NA
  
  set.values(mean.sor.r, as.numeric(colnames(beta.sor)),  as.numeric(mean.sor))
  
  ### Sne
  mean.sne.r <- dist$grid
  
  values(mean.sne.r) <- NA
  
  set.values(mean.sne.r, as.numeric(colnames(beta.sne)),  as.numeric(mean.sne))
  
  #### Raster PCoA ####
  
  ### Sim
  
  pcoa.sim.r1 <- dist$grid
  pcoa.sim.r2 <- dist$grid
  
  values(pcoa.sim.r1) <- NA
  values(pcoa.sim.r2) <- NA
  
  set.values(pcoa.sim.r1, as.numeric(names(pcoa.sim.1)),  as.numeric(pcoa.sim.1))
  set.values(pcoa.sim.r2, as.numeric(names(pcoa.sim.2)),  as.numeric(pcoa.sim.2))
  
  ### Sor
  
  pcoa.sor.r1 <- dist$grid
  pcoa.sor.r2 <- dist$grid
  
  values(pcoa.sor.r1) <- NA
  values(pcoa.sor.r2) <- NA
  
  set.values(pcoa.sor.r1, as.numeric(names(pcoa.sor.1)),  as.numeric(pcoa.sor.1))
  set.values(pcoa.sor.r2, as.numeric(names(pcoa.sor.2)),  as.numeric(pcoa.sor.2))
  
  ### Sne
  
  pcoa.sne.r1 <- dist$grid
  pcoa.sne.r2 <- dist$grid
  
  values(pcoa.sne.r1) <- NA
  values(pcoa.sne.r2) <- NA
  
  set.values(pcoa.sne.r1, as.numeric(names(pcoa.sne.1)),  as.numeric(pcoa.sne.1))
  set.values(pcoa.sne.r2, as.numeric(names(pcoa.sne.2)),  as.numeric(pcoa.sne.2))
  
  #### Raster PCoA ####
  
  ### Sim
  
  ndms.sim.r1 <- dist$grid
  ndms.sim.r2 <- dist$grid
  
  values(ndms.sim.r1) <- NA
  values(ndms.sim.r2) <- NA
  
  set.values(ndms.sim.r1, as.numeric(names(ndms.sim.1)),  as.numeric(ndms.sim.1))
  set.values(ndms.sim.r2, as.numeric(names(ndms.sim.2)),  as.numeric(ndms.sim.2))
  
  ### Sor
  
  ndms.sor.r1 <- dist$grid
  ndms.sor.r2 <- dist$grid
  
  values(ndms.sor.r1) <- NA
  values(ndms.sor.r2) <- NA
  
  set.values(ndms.sor.r1, as.numeric(names(ndms.sor.1)),  as.numeric(ndms.sor.1))
  set.values(ndms.sor.r2, as.numeric(names(ndms.sor.2)),  as.numeric(ndms.sor.2))
  
  ### Sne
  
  ndms.sne.r1 <- dist$grid
  ndms.sne.r2 <- dist$grid
  
  values(ndms.sne.r1) <- NA
  values(ndms.sne.r2) <- NA
  
  set.values(ndms.sne.r1, as.numeric(names(ndms.sne.1)),  as.numeric(ndms.sne.1))
  set.values(ndms.sne.r2, as.numeric(names(ndms.sne.2)),  as.numeric(ndms.sne.2))
  
  stats.phylo.beta <- data.frame(Eig.PCoA.1 = c(eig.sor[1], eig.sim[1], eig.sne[1]),
                                 Eig.PCoa.2 =c(eig.sor[2], eig.sim[2], eig.sne[2]),
                                 Stress.NDMS = c(ndms.stress.sor,ndms.stress.sim,ndms.stress.sne) )
  
  rownames(stats.phylo.beta) <- c("Total","Turnover","Nestedness")
  
  structure(
    list(distribution = dist$distribution,
         total.dissimilarity = as.dist(beta.sor),
         turnover = as.dist(beta.sim),
         nestedness = as.dist(beta.sne),
         rasters = list(mean.total = mean.sor.r,
                        mean.turnover = mean.sim.r,
                        mean.nestedness = mean.sne.r,
                        total.pcoa.1 = pcoa.sor.r1,
                        total.pcoa.2 = pcoa.sor.r2,
                        turnover.pcoa.1 = pcoa.sim.r1,
                        turnover.pcoa.2 = pcoa.sim.r2,
                        nestedness.pcoa.1 = pcoa.sne.r1,
                        nestedness.pcoa.2 = pcoa.sne.r2,
                        total.ndms.1 = ndms.sor.r1,
                        total.ndms.2 = ndms.sor.r2,
                        turnover.ndms.1 = ndms.sim.r1,
                        turnover.ndms.2 = ndms.sim.r2,
                        nestedness.ndms.1 = ndms.sne.r1,
                        nestedness.ndms.2 = ndms.sne.r2),
         stats = stats.phylo.beta,
         calculation.method = paste0(method, ifelse(normalize == T, ", with normalization", ", without normalization")),
         algorithm = "PhyloBeta"
      
    ),
    class = "vilma.beta"
  )
  
  
}
