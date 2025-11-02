#' Rao's Q (α) null models for vilma objects
#' @description
#' Generates null expectations for **within-community** Rao’s quadratic entropy
#' (RaoQ; phylogenetic α-diversity) using several randomization schemes, and
#' computes standardized effect sizes (SES) either **globally** or **per cell**.
#'
#' @param pd An object of class \code{vilma.pd} returned by \code{rao.calc()},
#'   containing the observed per-cell RaoQ table (\code{$pd.table}).
#' @param tree A rooted phylogeny of class \code{phylo} with branch lengths.
#' @param dist A \code{vilma.dist} object (e.g., from \code{points.to.raster()})
#'   that includes \code{$distribution} with columns \code{Cell}, \code{Sp}, and
#'   \code{$grid} used when rasterizing per-cell results.
#' @param iterations Integer; number of randomizations (default \code{999}).
#' @param method Character; whether to compute null statistics at the
#'   \code{"global"} level (summing RaoQ across cells) or per \code{"cell"}
#'   (default behavior resolves to \code{"global"} if multiple are supplied).
#' @param sampling Character; null model type:
#'   \itemize{
#'     \item \code{"taxa.label"} — Shuffle tip labels on the tree; community matrix fixed.
#'     \item \code{"range"} — Swap occurrences in a presence–absence matrix (preserves row/column sums).
#'     \item \code{"neighbor"} — Local swaps between adjacent cells on \code{dist$grid}.
#'     \item \code{"regional"} — Sample species for each cell from the regional pool (see \code{regional.weight}).
#'   }
#'   If multiple options are passed, defaults to \code{"taxa.label"}.
#' @param n.directions Neighborhood definition for the \code{"neighbor"} model:
#'   \code{"rook"}, \code{"bishop"}, or \code{"queen"} (default resolves to \code{"queen"}).
#' @param regional.weight Weighting of the regional pool for the \code{"regional"} model:
#'   \code{"uniform"} (equal), \code{"frequency"} (by number of records), or
#'   \code{"range"} (by number of occupied cells). Defaults to \code{"uniform"}.
#' @param abundance Logical; if \code{FALSE} (default) each Cell×Sp table is binarized
#'   prior to computing proportions; if \code{TRUE} the counts from \code{table(Cell, Sp)}
#'   are used as abundances.
#' @param scale01 Logical; if \code{TRUE} (default) patristic distances are divided
#'   by their maximum off-diagonal so \eqn{D \in [0,1]} (diagonal set to 0).
#'
#' @details
#' For each iteration, a randomized dataset (or tree, for \code{"taxa.label"}) is
#' produced, and \code{rao.calc()} is evaluated to obtain a null RaoQ per cell.
#' The function returns either:
#' \itemize{
#'   \item \strong{Global}: observed sum of RaoQ across cells, the vector of null sums,
#'     and \code{SES = (obs - mean(null)) / sd(null)} with a one-sided
#'     Monte Carlo p-value \eqn{( \#(null \ge obs) + 1 ) / (N + 1)}.
#'   \item \strong{Cell}: per-cell observed RaoQ, null mean and SD, SES per cell,
#'     and per-cell p-values computed the same way, plus a raster of SES.
#' }
#'
#' Species are restricted to the intersection between \code{tree$tip.label} and
#' \code{dist$distribution$Sp}. For \code{"neighbor"}, adjacency is derived from
#' \code{dist$grid} and the chosen \code{n.directions}. For \code{"regional"},
#' the number of species sampled in each cell equals its observed richness
#' (\code{pd$pd.table$SR}), with probabilities defined by \code{regional.weight}.
#'
#' @return An object of class \code{vilma.null}:
#' \describe{
#'   \item{\code{pd.obs}}{Observed total RaoQ (global mode) or \code{CellValues$PD} (cell mode).}
#'   \item{\code{null.pd}}{Vector of null totals (global mode) or omitted (cell mode).}
#'   \item{\code{SES}}{Standardized effect size (scalar in global mode; per-cell in cell mode).}
#'   \item{\code{Pvalue}}{Monte Carlo p-value(s).}
#'   \item{\code{CellValues}}{Data frame with per-cell results (cell mode).}
#'   \item{\code{Raster}}{Raster of SES (cell mode).}
#'   \item{\code{Iterations}}{Number of randomizations used.}
#'   \item{\code{Iter.table}}{Matrix of null RaoQ values (rows = cells, cols = iterations).}
#'   \item{\code{Method}}{The string \code{"global"} or \code{"cell"}.}
#' }
#'
#' @section Interpretation:
#' Negative SES indicates **phylogenetic clustering** (lower RaoQ than expected),
#' positive SES indicates **overdispersion** (higher RaoQ than expected), and
#' SES near 0 suggests random assembly under the chosen null model.
#'
#' @references
#' Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified approach.
#' \emph{Theoretical Population Biology}, 21, 24–43. \cr
#' Botta-Dukát, Z. (2005). Rao’s quadratic entropy as a measure of functional diversity.
#' \emph{Community Ecology}, 6, 283–290. \cr
#' Hardy, O. J. (2008). Testing the spatial phylogenetic structure of local communities:
#' statistical performances of different null models and test statistics on a locally neutral community.
#' \emph{Journal of Ecology}, 96(5), 914–926.
#'
#' @author
#' Omar Daniel Leon-Alvarado \url{https://leon-alvarado.weebly.com/}
#' J. Angel Soto-Centeno \url{https://www.mormoops.com}
#'
#' @seealso
#' \code{\link{rao.calc}} for observed α-diversity (RaoQ);
#' \code{\link{rao.beta}} for Rao β dissimilarity;
#' \code{\link{points_to_raster}} for building \code{vilma.dist}.
#'
#' @examples
#' \dontrun{
#' 
#' tree <- examplae_tree()
#' dist <- example_dist()
#' dist <- points_to_raster(dist)
#'
#' rao <- rao.calc(tree, dist)
#'
#' out.g <- rao.calc.null(rao, tree = tree, dist = dist,
#'                        iterations = 999, method = "global",
#'                        sampling = "taxa.label")
#' print(out.g)
#' plot(out.g)
#' view.vilma(out.g)
#'
#' # Per-cell SES under regional pool (frequency weights)
#' out.c <- rao.calc.null(pd = pd_obs, tree = tree, dist = dist,
#'                        iterations = 499, method = "cell",
#'                        sampling = "regional", regional.weight = "frequency")
#' }
#'
#' @export


rao.calc.null <- function(pd, tree, dist, iterations = 999, 
                          method = c("global","cell"), 
                          sampling = c("taxa.label","range","neighbor","regional"), 
                          n.directions = c("rook","bishop","queen"),
                          regional.weight = c("uniform","frequency","range"),
                          abundance = FALSE, scale01 = TRUE){
  
  ##############################################################
  #                        VERIFICATION                        #
  ##############################################################
  
  if (!inherits(pd, "vilma.pd")) {
    stop("Input 'pd' must be an object of class 'vilma.pd'.")
  }
  
  if(!inherits(tree, "phylo")){
    stop("Input 'tree' must be an objecto of class 'phylo'.")
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
  
  # Defaults when user passes c("a","b") or multiple choices
  if (length(method) != 1) method <- "global"
  if (length(sampling) != 1) sampling <- "taxa.label"
  if (length(n.directions) != 1) n.directions <- "queen"
  if (length(regional.weight) != 1) regional.weight <- "uniform"
  
  ##############################################################
  #                      END VERIFICATION                      #
  ##############################################################
  
  # Messaging (use scale01, not 'normalization')
  if (sampling == "neighbor"){
    message("Running null-model with: ", iterations, " iterations.\n",
            "Method: ", method, "\n",
            "Sampling: ", sampling, "\n",
            "Neighbor direction: ", n.directions, "\n",
            "Abundance: ", abundance, "\n",
            "Scale 0-1: ", scale01, "\n")
  } else if (sampling == "regional") {
    message("Running null-model with: ", iterations, " iterations.\n",
            "Method: ", method, "\n",
            "Sampling: ", sampling, "\n",
            "Weight: ", regional.weight, "\n",
            "Abundance: ", abundance, "\n",
            "Scale 0-1: ", scale01, "\n")
  } else {
    message("Running null-model with: ", iterations, " iterations.\n",
            "Method: ", method, "\n",
            "Sampling: ", sampling, "\n",
            "Abundance: ", abundance, "\n",
            "Scale 0-1: ", scale01, "\n")
  }
  
  ###############################
  #         Taxa Method         #
  ###############################
  
  if (sampling == "taxa.label") {
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    dist1 <- dist
    taxa.label <- tree$tip.label
    taxa.label <- intersect(taxa.label, dist$distribution$Sp)
    taxa.tree.pos <- which(tree$tip.label %in% taxa.label)
    
    set.seed(123)
    iter.tree <- tree
    iter.tree.out <- vector("list", iterations)
    
    message("\n","Iterating trees ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for (i in 1:iterations) {
      iter.tree$tip.label[taxa.tree.pos] <- sample(taxa.label, length(taxa.label))
      iter.tree.out[[i]] <- iter.tree
      setTxtProgressBar(pb, i)
    }
    
    message("\n","Calculating PD ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    for (i in 1:iterations) {
      pd.iter <- suppressMessages(
        rao.calc(iter.tree.out[[i]], dist = dist1, abundance = abundance, scale01 = scale01)
      )
      samples[, i] <- pd.iter$pd.table$PD
      setTxtProgressBar(pb, i)
    }
    colnames(samples) <- paste0("It.", seq_len(ncol(samples)))
  }
  
  ###############################
  #        Range Method         #
  ###############################
  
  if (sampling == "range") {
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    dist1 <- dist
    
    dist0 <- dist$distribution[, c(1,4)]
    dist01 <- dist0[which(dist0$Sp %in% intersect(dist0$Sp, tree$tip.label)), ]
    dist01 <- table(dist01$Cell, dist01$Sp)
    dist01 <- (dist01 > 0) * 1
    dist01 <- as.matrix(dist01)
    
    iter.dist <- vector("list", iterations)
    
    message("\n","Swaping distributions ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    for (i in 1:iterations) {
      iter.d <- swap.null(dist01)
      dist1$distribution <- return.vilma.dist(iter.d)
      iter.dist[[i]] <- dist1
      dist1$distribution <- NA
      setTxtProgressBar(pb, i)
    }
    
    message("\n","Calculating PD ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    for (i in 1:iterations) {
      pd.iter <- suppressMessages(
        rao.calc(tree, iter.dist[[i]], abundance = abundance, scale01 = scale01)
      )  # <-- closed the parenthesis here
      samples[, i] <- pd.iter$pd.table$PD
      setTxtProgressBar(pb, i)
    }
    colnames(samples) <- paste0("It.", seq_len(ncol(samples)))
  }
  
  ###############################
  #       Neighbor Method       #
  ###############################
  if (sampling == "neighbor") {
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    dist0 <- dist
    dist0$distribution <- dist0$distribution[
      which(dist0$distribution$Sp %in% intersect(dist0$distribution$Sp, tree$tip.label)), ]
    dist0$distribution <- dist0$distribution[
      which(dist0$distribution$Cell %in% as.integer(pd$pd.table$Cell)), ]
    
    dist1 <- dist0
    dist2 <- dist0
    iter.dist <- vector("list", iterations)
    
    iter <- (ncell(dist1$grid) * (log(ncell(dist1$grid)) + 0.577))
    iter <- iter * 2
    
    message("\n","Swaping neighbors ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    for (i in 1:iterations) {
      dist1$distribution$Sp <- as.factor(dist1$distribution$Sp)
      dist1$distribution$Cell <- as.integer(dist1$distribution$Cell)
      
      n.cells <- adjacent(x = dist1$grid,
                                 cells = 1:ncell(dist1$grid),
                                 directions = n.directions,
                                 pairs = TRUE)
      n.cells <- split(n.cells[,2], n.cells[,1])
      spp.cells <- split(dist1$distribution$Sp, dist1$distribution$Cell)
      
      for (j in 1:iter) {
        cell.a <- sample(names(spp.cells), 1)
        sp.a <- spp.cells[[cell.a]]
        if (length(sp.a) == 0) next
        
        neighs <- n.cells[[cell.a]]
        if (is.null(neighs) || length(neighs) == 0) next
        cell.b <- as.character(sample(neighs, 1))
        sp.b <- spp.cells[[cell.b]]
        if (length(sp.b) == 0) next
        
        sp.a.pick <- sample(sp.a, 1)
        sp.b.pick <- sample(sp.b, 1)
        
        if (!(sp.a.pick %in% sp.b) && !(sp.b.pick %in% sp.a)) {
          spp.cells[[cell.a]][spp.cells[[cell.a]] == sp.a.pick][1] <- sp.b.pick
          spp.cells[[cell.b]][spp.cells[[cell.b]] == sp.b.pick][1] <- sp.a.pick
        }
      }
      
      Sp <- character(0)
      Cell <- integer(0)
      for (x in seq_along(spp.cells)) {
        Sp <- c(Sp, as.character(spp.cells[[x]]))
        Cell <- c(Cell, as.integer(rep(names(spp.cells)[x], length(spp.cells[[x]]))))
      }
      
      dist2$distribution <- data.frame(Sp = Sp,
                                       Lon = rep(0, length(Sp)),
                                       Lat = rep(0, length(Sp)),
                                       Cell = Cell)
      iter.dist[[i]] <- dist2
      setTxtProgressBar(pb, i)
    }
    
    message("\n","Calculating PD ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    for (i in 1:iterations) {
      pd.iter <- suppressMessages(
        rao.calc(tree, iter.dist[[i]], abundance = abundance, scale01 = scale01)
      )  # <-- fixed 'abundace' -> 'abundance'
      samples[, i] <- pd.iter$pd.table$PD
      setTxtProgressBar(pb, i)
    }
    colnames(samples) <- paste0("It.", seq_len(ncol(samples)))
  }
  
  ###############################
  #        Regional Pool        #
  ###############################
  
  if (sampling == "regional") {
    samples <- matrix(NA, nrow = nrow(pd$pd.table), ncol = iterations)
    
    dist0 <- dist
    dist0$distribution <- dist0$distribution[
      which(dist0$distribution$Sp %in% intersect(dist0$distribution$Sp, tree$tip.label)), ]
    dist0$distribution <- dist0$distribution[
      which(dist0$distribution$Cell %in% as.integer(pd$pd.table$Cell)), ]
    
    dist1 <- dist0
    
    distM <- dist0$distribution
    tree.sp <- tree$tip.label
    obs.sp <- unique(distM$Sp)
    common.sp <- intersect(tree.sp, obs.sp)
    if (length(common.sp) == 0) stop("No overlapping species for regional null.")
    
    weight <- rep(1, length(common.sp))
    names(weight) <- common.sp
    
    if (regional.weight == "frequency") {
      species.freq <- table(distM$Sp)
      weight[names(species.freq)] <- as.numeric(species.freq[names(species.freq)])
    }
    if (regional.weight == "range") {
      occ.tab <- table(distM$Sp, distM$Cell)
      species.cells <- rowSums(occ.tab > 0)
      weight[names(species.cells)] <- as.numeric(species.cells[names(species.cells)])
    }
    prob_vec <- weight / sum(weight)
    
    pd_table <- pd$pd.table
    cell_ids <- pd_table$Cell
    richness_vec <- pd_table$SR
    
    iter.dist <- vector("list", iterations)
    
    message("\n","Getting distributions ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    
    for (it in seq_len(iterations)) {
      null_list <- vector("list", length = length(cell_ids))
      names(null_list) <- cell_ids
      for (i in seq_along(cell_ids)) {
        k <- richness_vec[i]
        if (k > 0) {
          null_list[[i]] <- sample(common.sp, size = k, replace = FALSE, prob = prob_vec)
        } else {
          null_list[[i]] <- character(0)
        }
      }
      # Build a vilma.dist-like distribution from the null_list
      null_dist <- return.vilma.dist2(null_list)
      dist1$distribution <- null_dist
      iter.dist[[it]] <- dist1
      setTxtProgressBar(pb, it)
    }
    
    message("\n","Calculating PD ...","\n")
    pb <- txtProgressBar(min = 0, max = iterations , style = 3)
    for (i in seq_len(iterations)) {
      pd_null <- suppressMessages(rao.calc(tree, iter.dist[[i]], abundance = abundance, scale01 = scale01)
      )
      samples[, i] <- pd_null$pd.table$PD
      setTxtProgressBar(pb, i)
    }
  }
  
  cat("\n")
  
  ##############################################################
  #                     NULL MODEL STATS                       #
  ##############################################################

  if (method == "global") {
    pd_obs <- sum(pd$pd.table$PD, na.rm = TRUE)
    null_pd <- as.vector(apply(samples, 2, sum, na.rm = TRUE))
    null_mean <- mean(null_pd)
    null_sd <- sd(null_pd)
  
    ## SAFE SES (avoid Inf when sd = 0)
  if (is.finite(null_sd) && null_sd == 0) {
    ses <- if (isTRUE(all.equal(pd_obs, null_mean))) 0 else NA_real_
  } else {
    ses <- (pd_obs - null_mean) / null_sd
  }
  
  ## Use only finite nulls in denominator
  n_eff <- sum(is.finite(null_pd))                                    # <<< ADDED
  if (n_eff == 0) {                                                   # <<< ADDED
    p_val <- NA_real_                                                 # <<< ADDED
  } else {                                                            # <<< ADDED
    p_val <- (sum(null_pd >= pd_obs, na.rm = TRUE) + 1) / (n_eff + 1) # <<< ADDED
  }                                                                   # <<< ADDED
  
  return(
    structure(
      list(pd.obs = pd_obs,
           null.pd = null_pd,
           SES = ses,
           Pvalue = p_val,
           Iterations = iterations,
           Iter.table = samples,
           Method = method),
      class = "vilma.null"
      )
    )
  }

  if (method == "cell") {
    pd_obs <- pd$pd.table$PD
    null_mean <- as.vector(apply(samples, 1, mean, na.rm = TRUE))
    null_sd <- as.vector(apply(samples, 1, sd,   na.rm = TRUE))
    ses <- (pd_obs - null_mean) / null_sd
  
    ## SAFE SES per-cell (avoid Inf when sd = 0)
    zero_sd <- is.finite(null_sd) & (null_sd == 0)
    same_as_mean <- abs(pd_obs - null_mean) < .Machine$double.eps^0.5
    ses[zero_sd &  same_as_mean] <- 0
    ses[zero_sd & !same_as_mean] <- NA_real_
  
    p_val <- numeric(nrow(samples))
  
  for (i in 1:nrow(samples)) {
    pd.obs1 <- pd_obs[i]
    null.pd <- samples[i, ]
    
    denom_i <- sum(is.finite(null.pd))                                # <<< ADDED
    if (denom_i == 0) {                                               # <<< ADDED
      p_val[i] <- NA_real_                                            # <<< ADDED
    } else {                                                          # <<< ADDED
      p_val[i] <- (sum(null.pd >= pd.obs1, na.rm = TRUE) + 1) / (denom_i + 1)  # <<< ADDED
    }                                                                 # <<< ADDED
  }
  
  cell.values <- data.frame(Cell = pd$pd.table$Cell,
                            PD = pd_obs,
                            NullPD = null_mean,
                            NullSD = null_sd,
                            SES = ses,
                            Pvalue = p_val)
  
  grid0 <- dist$grid
  values(grid0) <- NA
  set.values(grid0,
             as.numeric(as.character(cell.values$Cell)),
             as.numeric(ses))
  
  return(
    structure(
      list(CellValues = cell.values,
           Raster = grid0,
           Iterations = iterations,
           Iter.table = samples,
           Method = method),
      class = "vilma.null"
      )
    )
  }
}





